function [R, info] = eeg_ransac(EEG, varargin)
% EEG_RANSAC  Correlation of every channel with its RANSAC prediction, per
% window, as in the PREP pipeline's bad-channel test (findNoisyChannels).
%
%   R          = eeg_ransac(EEG)
%   [R, info]  = eeg_ransac(EEG, 'Name', value, ...)
%
% Each channel is predicted by spherical-spline interpolation from Draws
% random subsets of Fraction of the channels; the prediction is the median
% over the subsets. Because a few bad channels are rarely in most subsets,
% the prediction stays sensible even when neighbours are bad, which a
% leave-one-out prediction from all other channels does not. Per window, the
% recorded signal is compared with that prediction.
%
% OPTIONS
%   'Channels'      channels evaluated and used as predictors (default: all
%                   channels with a location). Unlocated channels are an error.
%   'Windows'       [start stop] sample pairs, one row per window. Default:
%                   PREP's non-overlapping windows of WindowSeconds.
%   'WindowSeconds' window length when Windows is not given (default 5)
%   'Draws'         number of random subsets (default 50)
%   'Fraction'      subset size as a fraction of Channels (default 0.25)
%   'Seed'          seed of the subset draw (default 435656, as PREP, so the
%                   same subsets are drawn every run)
%   'Style'         'prep'    (default) PREP's test: the data are low-passed
%                             at 45-50 Hz when srate > 100 and the correlation
%                             is uncentred (sum(x.*y)/sqrt(sum(x.^2)*sum(y.^2)))
%                   'detrend' the test of the Interpolation Clean button:
%                             unfiltered data, recorded and predicted signal
%                             linearly detrended per window, Pearson r
%
% OUTPUT
%   R      numel(Channels) x nWindows correlations
%   info   .chans      the channels, .windows [start stop] per window,
%          .sd         numel(Channels) x nWindows SD of predicted minus
%                      recorded signal (after the same filtering/detrending)
%          .subsetSize channels per random subset, .draws, .seed
%
% PREP flags a channel when R < 0.75 in more than 40% of the windows
% (ransacCorrelationThreshold, ransacUnbrokenTime); that decision is left to
% the caller.
%
% The algorithm follows findNoisyChannels.m, spherical_interpolate.m,
% design_fir.m and filtfilt_fast.m of the PREP pipeline (Bigdely-Shamlo et al.
% 2015, Front Neuroinform 9:16; VisLab/EEG-Clean-Tools; filter and spline code
% by C. Kothe, SCCN). Reimplemented, not copied: subset draw, spline
% interpolation, filter and median are matched to PREP 0.57.0.
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

p = inputParser;
p.addParameter('Channels', []);
p.addParameter('Windows', []);
p.addParameter('WindowSeconds', 5);
p.addParameter('Draws', 50);
p.addParameter('Fraction', 0.25);
p.addParameter('Seed', 435656);
p.addParameter('Style', 'prep');
p.parse(varargin{:});
o = p.Results;
style = lower(o.Style);
if ~ismember(style, {'prep', 'detrend'})
    error('eeg_ransac:Style', 'Style must be ''prep'' or ''detrend''.');
end
if EEG.trials > 1
    error('eeg_ransac:Epoched', 'eeg_ransac needs continuous data.');
end

chans = o.Channels;
if isempty(chans)
    chans = find(~cellfun(@isempty, {EEG.chanlocs.X}));
end
chans = chans(:)';
if any(cellfun(@isempty, {EEG.chanlocs(chans).X}))
    error('eeg_ransac:NoLocation', 'Every channel in Channels needs a location.');
end
m = numel(chans);
subsetSize = round(o.Fraction * m);
if m < subsetSize + 1 || m < 3 || subsetSize < 2
    error('eeg_ransac:TooFewChannels', ...
        '%d channels are too few for RANSAC with a %.2f subset fraction.', m, o.Fraction);
end

% samples x channels, as PREP
X = double(EEG.data(chans, :))';
if strcmp(style, 'prep') && EEG.srate > 100
    B = designFir(100, [2*[0 45 50]/EEG.srate 1], [1 1 0 0]);
    X = filtfiltFast(B, X);
end

% windows
if isempty(o.Windows)
    frames  = o.WindowSeconds * EEG.srate;
    starts  = 1:frames:(size(X, 1) - frames);         % PREP drops the tail
    windows = [starts(:), starts(:) + frames - 1];
else
    windows = round(o.Windows);
end
nWin = size(windows, 1);

% one reconstruction matrix per random subset: rows of the subset map to
% every channel (m x m), all subsets side by side (m x m*Draws)
locs = [[EEG.chanlocs(chans).X]; [EEG.chanlocs(chans).Y]; [EEG.chanlocs(chans).Z]];
subsets = randomSubsets(m, subsetSize, o.Draws, o.Seed);
P = zeros(m, m * o.Draws);
for k = 1:o.Draws
    s   = subsets(k, :);
    tmp = zeros(m);
    tmp(s, :) = sphericalInterpolate(locs(:, s), locs)';
    P(:, (k-1)*m + (1:m)) = tmp;
end

R  = ones(m, nWin);
SD = zeros(m, nWin);
for w = 1:nWin
    XX = X(windows(w,1):windows(w,2), :);
    n  = size(XX, 1);
    YY = sort(reshape(XX * P, n, m, o.Draws), 3);
    YY = YY(:, :, round(end/2));                       % PREP's median
    if strcmp(style, 'prep')
        R(:, w) = (sum(XX.*YY) ./ (sqrt(sum(XX.^2)) .* sqrt(sum(YY.^2))))';
        SD(:, w) = std(YY - XX)';
    else
        real = detrend(XX);
        surr = detrend(YY);
        R(:, w)  = diag(corr(surr, real));
        SD(:, w) = std(surr - real)';
    end
end

info = struct('chans', chans, 'windows', windows, 'sd', SD, ...
    'subsetSize', subsetSize, 'draws', o.Draws, 'seed', o.Seed);
end

% -----------------------------------------------------------------------------
function subsets = randomSubsets(m, subsetSize, draws, seed)
% PREP's draw: subsetSize channels without replacement per subset, from a
% Mersenne twister with a fixed seed.
stream  = RandStream('mt19937ar', 'Seed', seed);
subsets = zeros(draws, subsetSize);
for k = 1:draws
    pool = 1:m;
    for j = 1:subsetSize
        pick = round(1 + (numel(pool) - 1) .* rand(stream));
        subsets(k, j) = pool(pick);
        pool(pick) = [];
    end
end
end

% -----------------------------------------------------------------------------
function W = sphericalInterpolate(src, dest)
% Spherical-spline weights from positions src (3 x nSrc) to dest (3 x nDest):
% lambda 1e-5, order 4, Legendre series summed until converged (<= 500 terms).
lambda = 1e-5;
src  = src  ./ sqrt(sum(src.^2));
dest = dest ./ sqrt(sum(dest.^2));
Gss = splineG(src' * src);
Gds = splineG(dest' * src);
Gss = Gss + lambda * eye(size(Gss));
C   = [Gss ones(size(Gss, 1), 1); ones(1, size(Gss, 2)) 0];
iC  = pinv(C);
W   = [Gds ones(size(Gds, 1), 1)] * iC(:, 1:end-1);
end

% -----------------------------------------------------------------------------
function G = splineG(x)
% sum_n (2n+1) P_n(x) / (n(n+1))^4 / (4 pi), element by element stopped when
% the running means of the absolute change of this sum and of the matching
% surface-Laplacian sum (terms times n(n+1)) both drop below eps, as PREP.
order = 4;
tol   = eps;
Pns1 = ones(size(x));
Pn   = x;
G    = 3 * Pn / (2^order);
dG   = abs(G);
dH   = abs(2 * G);
active = true(size(x));
for n = 2:500
    Pns2 = Pns1;
    Pns1 = Pn;
    Pn   = ((2*n - 1) * x .* Pns1 - (n - 1) * Pns2) / n;
    step = ((2*n + 1) * Pn) / ((n*n + n)^order);
    step(~active) = 0;
    G  = G + step;
    dG(active) = (abs(step(active)) + dG(active)) / 2;
    dH(active) = (abs((n*n + n) * step(active)) + dH(active)) / 2;
    active = active & ~(dG < tol & dH < tol);
    if ~any(active(:))
        break
    end
end
G = G / (4*pi);
end

% -----------------------------------------------------------------------------
function B = designFir(N, F, A)
% Frequency-sampled linear-phase FIR with a Hamming window (PREP design_fir).
nfft = max(512, 2^ceil(log(N)/log(2)));
W    = 0.54 - 0.46*cos(2*pi*(0:N)/N);
F    = interp1(round(F*nfft), A, (0:nfft), 'pchip');
F    = F .* exp(-(0.5*N)*sqrt(-1)*pi*(0:nfft)./nfft);
B    = real(ifft([F conj(F(end-1:-1:2))]));
B    = B(1:N+1) .* W;
end

% -----------------------------------------------------------------------------
function X = filtfiltFast(B, X)
% Zero-phase filtering with reflected padding (PREP filtfilt_fast, A = 1).
w = numel(B);
t = size(X, 1);
X = [2*X(1,:) - X(1 + mod(((w+1):-1:2) - 1, t), :); X; ...
     2*X(t,:) - X(1 + mod(((t-1):-1:(t-w)) - 1, t), :)];
X = filter(B, 1, X);  X = X(end:-1:1, :);
X = filter(B, 1, X);  X = X(end:-1:1, :);
X([1:w, t + w + (1:w)], :) = [];
end
