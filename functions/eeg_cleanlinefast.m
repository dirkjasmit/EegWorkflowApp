function [EEG, info] = eeg_cleanlinefast(EEG, varargin)
% EEG_CLEANLINEFAST  CleanLine's line-noise removal, every channel at once.
%
%   [EEG, info] = eeg_cleanlinefast(EEG, 'LineFrequencies', 60*(1:4), ...)
%
% Same algorithm and same numbers as the CleanLine plugin in its 'newversion'
% form (cleanLineNoise / removeLinesMovingWindow / fitSignificantFrequencies
% of Cleanline 2.1, which PREP also uses): in sliding windows, a multitaper
% Thomson F-test looks for a significant sinusoid near each line frequency,
% and the fitted sinusoid is subtracted, with a sigmoidal cross-fade over the
% overlap between windows; this repeats until a pass no longer lowers the
% power at that frequency, or MaximumIterations is reached.
%
% The plugin walks through the channels one at a time (a parfor over
% channels), and within each channel it redoes the tapered FFT per window and
% per iteration. The maths in testSignificantFrequencies is written for many
% channels at once but is never called that way. Here every window is done
% for ALL channels in one batched FFT, which is where the time goes: on a
% 127-channel, 6-minute recording that is one 1024 x 7 x 127 transform per
% window instead of 127 separate ones.
%
% OPTIONS (name/value; the defaults are CleanLine's own)
%   'LineFrequencies'   frequencies to remove (default 60*(1:4)); a single
%                       value is expanded to its first four harmonics, as
%                       pop_cleanline does
%   'Channels'          channels to clean (default all)
%   'p'                 significance level of the F-test (default 0.01)
%   'fScanBandWidth'    the exact line is looked for within +/- half of this
%                       around each frequency (default 2)
%   'taperBandWidth'    multitaper bandwidth in Hz (default 2)
%   'taperWindowSize'   sliding window in seconds (default 4)
%   'taperWindowStep'   step between windows in seconds (default 1)
%   'tau'               steepness of the cross-fade (default 100)
%   'pad'               FFT padding exponent (default 0)
%   'fPassBand'         frequency band kept (default [0 Fs/2])
%   'maximumIterations' passes over the data (default 10)
%   'Verbose'           print a line per iteration (default false)
%
% OUTPUT
%   EEG    cleaned dataset; EEG.etc.cleanline holds info
%   info   .frequencies  the line frequencies asked for
%          .iterations   passes actually used
%          .removed      per channel, the proportion of its variance removed
%          .active       channels x frequencies, still being fitted at the end
%
% Verified against the plugin: same data in, the two outputs agree to
% rounding (see the header of the test in the repository).
%
% Reference: Mitra & Bokil, Observed Brain Dynamics (2008), chapter 7; the
% implementation follows Cleanline 2.1 by T. Mullen (GPL), reorganised to run
% over channels at once.
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

p = inputParser;
p.addParameter('LineFrequencies', 60*(1:4));
p.addParameter('Channels', []);
p.addParameter('p', 0.01);
p.addParameter('fScanBandWidth', 2);
p.addParameter('taperBandWidth', 2);
p.addParameter('taperWindowSize', 4);
p.addParameter('taperWindowStep', 1);
p.addParameter('tau', 100);
p.addParameter('pad', 0);
p.addParameter('fPassBand', []);
p.addParameter('maximumIterations', 10);
p.addParameter('Verbose', false);
p.parse(varargin{:});
o = p.Results;

Fs = EEG.srate;
if isempty(o.fPassBand), o.fPassBand = [0 Fs/2]; end
f0 = o.LineFrequencies(:)';
if numel(f0) == 1
    f0 = f0 * (1:4);                      % as pop_cleanline expands it
end
f0 = f0(f0 < Fs/2);
chans = o.Channels;
if isempty(chans), chans = 1:EEG.nbchan; end
chans = chans(:)';

X0 = double(EEG.data(chans, :))';          % samples x channels
X  = X0;
[N, C] = size(X);
F  = numel(f0);

% tapers, once (checkTapers of the plugin: TW = halfBandwidth * windowSize)
Nwin   = round(Fs * o.taperWindowSize);
TW     = (o.taperBandWidth/2) * o.taperWindowSize;
K      = floor(2*TW - 1);
tapers = dpss(Nwin, TW, K) * sqrt(Fs);
nfft   = max(2^(nextpow2(Nwin) + o.pad), Nwin);
fax    = (0:nfft-1) * Fs / nfft;
findx  = find(fax >= o.fPassBand(1) & fax <= o.fPassBand(end));
f      = fax(findx);
Kodd   = 1:2:K;
Keven  = 2:2:K;
H0     = sum(tapers(:, Kodd), 1);          % 1 x Kodd, the same for every channel
H0sq   = sum(H0.^2);
sig    = finv(1 - o.p, 2, 2*K - 2);

% windows and the sigmoidal cross-fade over their overlap
Nstep    = round(o.taperWindowStep * Fs);
Noverlap = Nwin - Nstep;
winstart = 1:Nstep:(N - Nwin + 1);
nw       = numel(winstart);
if nw < 1
    error('eeg_cleanlinefast:TooShort', 'The data are shorter than one window.');
end
x      = (1:Noverlap)';
smooth = 1 ./ (1 + exp(-o.tau .* (x - Noverlap/2) / Noverlap));
last   = winstart(nw) + Nwin - 1;
t      = (0:Nwin-1)' / Fs;

% which frequencies are still being fitted, per channel
active = true(C, F);
spec   = 10*log10(segmentSpectrum(X, Fs, o, tapers, nfft, findx));   % f x C
fidx   = zeros(1, F);
for k = 1:F
    [~, fidx(k)] = min(abs(f - f0(k)));
end

iter = 0;
for iteration = 1:o.maximumIterations
    iter = iteration;
    datafit = zeros(last, C);
    fitted  = false(C, F);
    prev    = zeros(Nwin, C);
    for n = 1:nw
        idx  = winstart(n):(winstart(n) + Nwin - 1);
        J    = taperedFft(X(idx, :), tapers, nfft, Fs);              % nfft x K x C
        Jp   = J(findx, Kodd, :);
        A    = reshape(sum(Jp .* H0, 2), [], C) / H0sq;              % f x C
        Jhat = reshape(A, [], 1, C) .* H0;
        num  = (K - 1) * abs(A).^2 * H0sq;
        den  = reshape(sum(abs(Jp - Jhat).^2, 2), [], C) + ...
               reshape(sum(abs(J(findx, Keven, :)).^2, 2), [], C);
        Fval = num ./ den;
        A    = A * Fs;

        fit = zeros(Nwin, C);
        for k = 1:F
            cols = find(active(:, k))';
            if isempty(cols), continue, end
            lo = findNearest(f, f0(k) - o.fScanBandWidth/2);
            hi = findNearest(f, f0(k) + o.fScanBandWidth/2);
            band = lo:hi;
            Fb = Fval(band, cols);
            Fb(Fb < sig) = 0;
            [best, where] = max(Fb, [], 1);
            take = best > 0;                                         % significant
            if ~any(take), continue, end
            cols  = cols(take);
            rows  = band(where(take));
            lin   = sub2ind(size(A), rows(:), cols(:));
            amps  = A(lin).';                                        % 1 x nTake
            freqs = f(rows);
            fit(:, cols) = fit(:, cols) + 2*real(exp(1i*2*pi*t*freqs) .* amps);
            fitted(cols, k) = true;
        end
        if n > 1
            fit(1:Noverlap, :) = smooth .* fit(1:Noverlap, :) + ...
                (1 - smooth) .* prev(Nwin - Noverlap + 1:Nwin, :);
        end
        datafit(idx, :) = fit;
        prev = fit;
    end

    X(1:last, :) = X(1:last, :) - datafit;

    % keep a frequency only while removing it still lowers the power there
    if any(fitted(:))
        newspec = 10*log10(segmentSpectrum(X, Fs, o, tapers, nfft, findx));
        drop = (spec(fidx, :) - newspec(fidx, :))' < 0;              % C x F
        active = active & fitted & ~drop;
        spec = newspec;
    else
        active(:) = false;
    end
    if o.Verbose
        fprintf('eeg_cleanlinefast: pass %d, %d channel-frequencies still active\n', ...
            iteration, nnz(active));
    end
    if ~any(active(:))
        break
    end
end

EEG.data(chans, :) = cast(X', 'like', EEG.data);
EEG.icaact = [];
removed = sum((X0 - X).^2, 1) ./ max(sum(X0.^2, 1), eps);
info = struct('frequencies', f0, 'iterations', iter, 'removed', removed, 'active', active);
if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc), EEG.etc = struct(); end
EEG.etc.cleanline = info;
EEG.history = [EEG.history newline sprintf( ...
    'EEG = eeg_cleanlinefast(EEG, ''LineFrequencies'', %s); %% %d passes', mat2str(f0), iter)];
end

% -----------------------------------------------------------------------------
function J = taperedFft(data, tapers, nfft, Fs)
% mtfftc of the plugin: samples x channels -> nfft x tapers x channels.
J = fft(reshape(data, size(data,1), 1, []) .* tapers, nfft, 1) / Fs;
end

% -----------------------------------------------------------------------------
function i = findNearest(f, value)
[~, i] = min(abs(f - value));
end

% -----------------------------------------------------------------------------
function S = segmentSpectrum(X, Fs, o, tapers, nfft, findx)
% calculateSegmentSpectrum of the plugin, for every channel at once: the data
% are cut into non-overlapping windows of taperWindowSize and the multitaper
% power spectra are averaged over those windows.
Nwin = size(tapers, 1);
nseg = floor(size(X, 1) / Nwin);
if nseg < 1
    S = ones(numel(findx), size(X, 2));
    return
end
S = zeros(numel(findx), size(X, 2));
for s = 1:nseg
    idx = (s-1)*Nwin + (1:Nwin);
    J   = taperedFft(X(idx, :), tapers, nfft, Fs);
    S   = S + reshape(mean(abs(J(findx, :, :)).^2, 2), [], size(X, 2));
end
S = S / nseg;
end
