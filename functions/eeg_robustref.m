function [EEG, info] = eeg_robustref(EEG, varargin)
% EEG_ROBUSTREF  Robust average reference of the PREP pipeline: the average
% is computed without the bad channels, which are interpolated first, and
% detection and referencing are repeated until the set of bad channels stops
% changing.
%
%   [EEG, info] = eeg_robustref(EEG, 'Channels', chans, ...)
%
% An ordinary average reference is pulled towards the bad channels: a single
% noisy electrode ends up in every other channel with the opposite sign. PREP
% therefore estimates the reference from an interpolated, bad-channel-free
% version of the data.
%
%   1. A working copy is high-passed at 1 Hz (detection only) and referenced
%      to the median over channels, as a first robust estimate.
%   2. Bad channels are detected on that copy (see below) and added to the
%      bad set.
%   3. They are interpolated (spherical splines) in the high-passed data, the
%      mean over channels of that interpolated data is the new reference, and
%      the working copy becomes the high-passed data minus it.
%   4. Steps 2-3 repeat until the bad set is unchanged or MaxIterations is
%      reached.
%   5. The reference is then computed the same way from the ORIGINAL (not
%      high-passed) data and subtracted from every channel. With
%      'InterpolateBad' the bad channels are left interpolated in the output,
%      as PREP does ('post-reference' interpolation order).
%
% The tests, with PREP's defaults:
%   deviation    robust z (0.7413*iqr) of the channel SD above 'Deviation' (5)
%   correlation  a channel whose largest correlation with any other channel is
%                below 'Correlation' (0.4) in more than 'BadTime' (0.01) of
%                the 1 s windows
%   HF noise     robust z of the ratio (noise above 45 Hz)/(signal below) above
%                'HighFrequency' (5). Only with a sampling rate above 100 Hz.
%   RANSAC       eeg_ransac: correlation with the prediction from random
%                channel subsets below 0.75 in more than 40% of 5 s windows
%
% OPTIONS (name/value)
%   'Channels'       channels referenced and tested (default: located channels)
%   'MaxIterations'  maximum detect/reference rounds (default 4, as PREP)
%   'Deviation'      robust z for the deviation test, 0 = off (default 5)
%   'Correlation'    minimum correlation with any other channel (default 0.4)
%   'BadTime'        fraction of windows a channel may fail it (default 0.01)
%   'HighFrequency'  robust z for the high-frequency noise test (default 5)
%   'Ransac'         run the RANSAC test (default true)
%   'InterpolateBad' leave the bad channels interpolated in the output
%                    (default true, as PREP); false restores them unreferenced
%                    values minus the reference
%   'Verbose'        print a line per iteration (default false)
%
% OUTPUT
%   EEG    referenced dataset; EEG.etc.robustref holds info
%   info   .bad          bad channel indices (into EEG)
%          .labels       their labels
%          .iterations   rounds used
%          .reference    the reference signal (1 x pnts)
%          .interpolated true when the bad channels were interpolated
%          .tests        struct with the per-test channel lists
%
% Reference: Bigdely-Shamlo N, Mullen T, Kothe C, Su KM, Robbins KA (2015).
% The PREP pipeline: standardized preprocessing for large-scale EEG analysis.
% Front Neuroinform 9:16. Reimplemented after robustReference.m and
% performReference.m of PREP 0.57.0; the detection tests are close to
% findNoisyChannels but not identical (see eeg_ransac for the RANSAC part).
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

p = inputParser;
p.addParameter('Channels', []);
p.addParameter('MaxIterations', 4);
p.addParameter('Deviation', 5);
p.addParameter('Correlation', 0.4);
p.addParameter('BadTime', 0.01);
p.addParameter('HighFrequency', 5);
p.addParameter('Ransac', true);
p.addParameter('InterpolateBad', true);
p.addParameter('Verbose', false);
p.parse(varargin{:});
o = p.Results;
if EEG.trials > 1
    error('eeg_robustref:Epoched', 'eeg_robustref needs continuous data.');
end
chans = o.Channels;
if isempty(chans)
    chans = find(~cellfun(@isempty, {EEG.chanlocs.X}));
end
chans = chans(:)';
noLoc = chans(cellfun(@isempty, {EEG.chanlocs(chans).X}));
if ~isempty(noLoc)
    error('eeg_robustref:NoLocation', ...
        ['%d of the channels have no location (%s): they can neither be interpolated nor ' ...
         'predict the others. Look their locations up, or leave them out of Channels.'], ...
        numel(noLoc), strjoin({EEG.chanlocs(noLoc).labels}, ' '));
end
if numel(chans) < 4
    error('eeg_robustref:TooFewChannels', 'At least 4 located channels are needed.');
end

% high-passed copy for the detection only
H = EEG;
H.data = highpass1Hz(double(EEG.data(chans, :)), EEG.srate);
H.chanlocs = EEG.chanlocs(chans);
H.nbchan = numel(chans);
H.icaweights = []; H.icasphere = []; H.icawinv = []; H.icaact = []; H.icachansind = [];

% first estimate: median over channels
W = H.data - median(H.data, 1);

bad = false(1, numel(chans));
tests = struct('deviation', [], 'correlation', [], 'highFrequency', [], 'ransac', []);
it = 0;
while it <= o.MaxIterations
    [hit, tests] = noisyChannels(W, H, o, tests);
    newBad = bad | hit;
    if it > 0 && isequal(newBad, bad)
        break
    end
    bad = newBad;
    if all(bad)
        error('eeg_robustref:AllBad', 'Every channel failed the tests; no robust reference.');
    end
    W = H.data - referenceOf(H.data, bad, H.chanlocs);
    it = it + 1;
    if o.Verbose
        fprintf('eeg_robustref: iteration %d, %d bad channels\n', it, nnz(bad));
    end
end

% the reference itself, from the original data
X   = double(EEG.data(chans, :));
ref = referenceOf(X, bad, EEG.chanlocs(chans));
if o.InterpolateBad && any(bad)
    X = interpolateBad(X, bad, EEG.chanlocs(chans));
end
EEG.data(chans, :) = cast(X - ref, 'like', EEG.data);
rest = setdiff(1:EEG.nbchan, chans);
if ~isempty(rest)
    EEG.data(rest, :) = cast(double(EEG.data(rest, :)) - ref, 'like', EEG.data);
end
EEG.ref   = 'averageref';
EEG.icaact = [];
[EEG.chanlocs.ref] = deal('average');

info = struct('bad', chans(bad), 'labels', {{EEG.chanlocs(chans(bad)).labels}}, ...
    'iterations', it, 'reference', ref, 'interpolated', o.InterpolateBad && any(bad), ...
    'tests', tests);
info.tests = structfun(@(v) chans(v), info.tests, 'UniformOutput', false);
if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc), EEG.etc = struct(); end
EEG.etc.robustref = info;
EEG.history = [EEG.history newline sprintf( ...
    'EEG = eeg_robustref(EEG); %% PREP robust average reference, %d iterations, %d bad channels', it, nnz(bad))];
end

% -----------------------------------------------------------------------------
function ref = referenceOf(X, bad, chanlocs)
% Mean over all channels after interpolating the bad ones.
if any(bad)
    X = interpolateBad(X, bad, chanlocs);
end
ref = mean(X, 1);
end

% -----------------------------------------------------------------------------
function X = interpolateBad(X, bad, chanlocs)
% Spherical-spline interpolation (as eeg_interp) of the bad channels from the
% good ones, as one matrix operation.
xyz = [[chanlocs.X]; [chanlocs.Y]; [chanlocs.Z]]';
xyz = xyz ./ sqrt(sum(xyz.^2, 2));
G    = splineG(xyz);
good = find(~bad);
bd   = find(bad);
ng   = numel(good);
M    = [G(good, good); ones(1, ng)];
W    = G(bd, good) * pinv(M);
W    = W(:, 1:end-1);
W    = W - (sum(W, 2) - 1) / ng;          % data centred on the mean of the good channels
X(bd, :) = W * X(good, :);
end

% -----------------------------------------------------------------------------
function g = splineG(xyz)
% computeg of eeg_interp (m = 4, 7 Legendre terms) for all pairs.
EI = 1 - sqrt(max(0, sum(xyz.^2, 2) + sum(xyz.^2, 2)' - 2*(xyz*xyz')));
g    = zeros(size(EI));
Pnm1 = ones(size(EI));
Pn   = EI;
for nn = 1:7
    g = g + ((2*nn+1)/(nn^4*(nn+1)^4)) * Pn;
    Pnp1 = ((2*nn+1).*EI.*Pn - nn.*Pnm1) / (nn+1);
    Pnm1 = Pn;
    Pn   = Pnp1;
end
g = g/(4*pi);
end

% -----------------------------------------------------------------------------
function [bad, tests] = noisyChannels(W, H, o, tests)
% PREP's channel tests on the referenced, high-passed data W.
fs = H.srate;
n  = size(W, 1);
bad = false(1, n);

% deviation: robust z of the channel SD
if o.Deviation > 0
    sd = 0.7413 * iqr(W, 2)';
    z  = (sd - median(sd)) / (0.7413 * iqr(sd));
    hit = abs(z) > o.Deviation | isnan(z);
    tests.deviation = union(tests.deviation, find(hit));
    bad = bad | hit;
end

% correlation: largest correlation with any other channel, per 1 s window
if o.Correlation > 0
    L = round(fs);
    nw = floor(size(W, 2) / L);
    if nw > 1
        low = 0;
        for w = 1:nw
            C = corrcoef(W(:, (w-1)*L + (1:L))');
            C(1:n+1:end) = 0;
            low = low + (max(abs(C), [], 2)' < o.Correlation);
        end
        hit = low / nw > o.BadTime;
        tests.correlation = union(tests.correlation, find(hit));
        bad = bad | hit;
    end
end

% high-frequency noise: power above 45 Hz relative to the rest
if o.HighFrequency > 0 && fs > 100
    lo = lowpass45Hz(W, fs);
    r  = mad(W - lo, 1, 2)' ./ max(mad(lo, 1, 2)', eps);
    z  = (r - median(r)) / (1.4826 * mad(r, 1));
    hit = z > o.HighFrequency | isnan(z);
    tests.highFrequency = union(tests.highFrequency, find(hit));
    bad = bad | hit;
end

% RANSAC predictability
if o.Ransac
    try
        R = eeg_ransac(setData(H, W), 'Channels', 1:n);
        frac = mean(R < 0.75, 2)';
        hit  = frac > 0.4;
        tests.ransac = union(tests.ransac, find(hit));
        bad = bad | hit;
    catch
    end
end
end

% -----------------------------------------------------------------------------
function E = setData(E, X)
E.data = X;
end

% -----------------------------------------------------------------------------
function Y = highpass1Hz(X, fs)
% Zero-phase 1 Hz high-pass (PREP detrends at 1 Hz before detection).
Y = X - movmean(X, max(3, round(fs)), 2);     % boxcar high-pass, ~1 Hz
end

% -----------------------------------------------------------------------------
function Y = lowpass45Hz(X, fs)
Y = movmean(X, max(3, round(fs/45)), 2);
end
