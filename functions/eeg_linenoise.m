function [EEG, info] = eeg_linenoise(EEG, varargin)
% EEG_LINENOISE  Remove line noise by regressing out sinusoids of estimated
% frequency, fitted over whole stretches of data at once.
%
%   [EEG, info] = eeg_linenoise(EEG, 'LineFreq', 50, ...)
%
% Mains interference is a sinusoid. Within a stretch of continuous data its
% frequency is (nearly) constant, so it is modelled here over the whole
% stretch rather than in short sliding windows:
%
%   1. The FFT of the complete stretch (all channels at once) gives the peak
%      closest to the nominal line frequency. With T seconds of data the
%      frequency resolution is 1/T Hz, so the peak is already close.
%   2. The frequency is refined by a grid search around that peak, zooming in
%      a few times. For every candidate frequency f, a cosine and a sine of f
%      are regressed out of ALL channels in one least-squares step; the
%      candidate explaining most variance summed over channels wins.
%      Regressing both a sine and a cosine yields the exact best phase and
%      amplitude for every channel, so no search over phase is needed.
%   3. The best cosine/sine pair is regressed out. Every channel uses the same
%      two regressors (same frequency); only their weights differ per channel,
%      because line noise reaches each electrode with its own amplitude and
%      phase.
%
% Harmonics are handled the same way, each with its own refined frequency.
% A frequency is only removed when there is a clear peak: the FFT peak must
% stand 'MinPeakDb' above the median power of the surrounding bins (1.5 to 5 Hz
% away). Otherwise nothing is removed for that piece and it is reported as
% skipped; without a peak the search would just fit the largest noise bin.
% Nothing loops over channels: every step is a matrix operation on the full
% channels x samples data. Any tiny mismatch in the nominal sampling rate is
% absorbed into the estimated frequency.
%
% Data are split into stretches at boundary events (removed data breaks the
% phase) and, for epoched data, per epoch. Long stretches are split further
% into pieces of 'SegmentLength' seconds, because the mains frequency is not
% constant: on a 5.4 min recording the estimated frequency ranged from 59.99
% to 60.02 Hz between 20 s pieces. A single sinusoid over the whole recording
% then removed almost nothing (60 Hz peak 3.9 dB -> 3.9 dB), whereas 20 s
% pieces removed it to the noise floor (0.5 dB). Pieces must stay long enough
% for a fine frequency estimate (resolution 1/SegmentLength Hz before the grid
% search).
%
% INPUT
%   EEG             EEGLAB dataset
%
% OPTIONS (name/value)
%   'LineFreq'      nominal line frequency in Hz (default 50)
%   'Harmonics'     multiples of LineFreq to remove, e.g. 1 or [1 2 3]
%                   (default 1). Multiples at or above Nyquist are skipped.
%   'SearchHz'      the FFT peak is searched within +/- this many Hz of each
%                   nominal frequency (default 0.5)
%   'MinPeakDb'     minimum height of that peak above the median power 1.5-5 Hz
%                   around it, in dB, for the frequency to be removed
%                   (default 5). Measured: real 60 Hz mains 6-22 dB per 20 s
%                   piece; no line noise at most 1.8 dB with 32 or more
%                   channels, but up to 6 dB with only 4 channels, where the
%                   channel-averaged spectrum is noisier.
%   'GridPoints'    points per grid-search level (default 7)
%   'Levels'        number of zoom levels of the grid search (default 4)
%   'SegmentLength' maximum piece length in seconds; Inf = whole stretch
%                   between boundaries (default 20)
%   'Channels'      channel indices to clean (default all)
%   'Verbose'       print a summary (default true)
%
% OUTPUT
%   EEG             cleaned dataset; EEG.etc.linenoise holds info
%   info            struct array, one element per stretch x harmonic:
%                     .segment    [first last] sample
%                     .nominal    nominal frequency (Hz)
%                     .freq       estimated frequency (Hz); FFT peak when skipped
%                     .peakDb     peak height above its surroundings (dB)
%                     .skipped    true when no clear peak: nothing removed
%                     .removed    per channel, the proportion of that
%                                 channel's variance removed by this sinusoid
%                     .amplitude  per-channel amplitude of the removed sinusoid
%
% EXAMPLE
%   EEG = eeg_linenoise(EEG, 'LineFreq', 60, 'Harmonics', 1:3);
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

p = inputParser;
p.FunctionName = 'eeg_linenoise';
addRequired(p,  'EEG', @isstruct);
addParameter(p, 'LineFreq', 50, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Harmonics', 1, @(x) isnumeric(x) && all(x >= 1) && all(x == round(x)));
addParameter(p, 'SearchHz', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'MinPeakDb', 5, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'GridPoints', 7, @(x) isnumeric(x) && isscalar(x) && x >= 3);
addParameter(p, 'Levels', 4, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'SegmentLength', 20, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Channels', [], @isnumeric);
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));
parse(p, EEG, varargin{:});
o = p.Results;

fs    = EEG.srate;
chans = o.Channels;
if isempty(chans)
    chans = 1:EEG.nbchan;
end
freqs = o.LineFreq * unique(o.Harmonics(:)');
freqs = freqs(freqs < fs/2 - o.SearchHz);
info  = struct('segment', {}, 'nominal', {}, 'freq', {}, 'peakDb', {}, 'skipped', {}, 'removed', {}, 'amplitude', {});
if isempty(freqs)
    warning('eeg_linenoise:nothingToDo', 'No line frequencies below Nyquist.');
    return
end

isEpoched = EEG.trials > 1;
nEp       = size(EEG.data, 3);
for ep = 1:nEp
    % ---- stretches of continuous data ---------------------------------------
    if isEpoched
        cuts = [1, EEG.pnts + 1];
    else
        cuts = 1;
        if isfield(EEG, 'event') && ~isempty(EEG.event)
            isB  = strcmpi({EEG.event.type}, 'boundary');
            lat  = round([EEG.event(isB).latency]);
            cuts = [cuts, lat(lat > 1 & lat <= EEG.pnts)];
        end
        cuts = unique([cuts, EEG.pnts + 1]);
    end
    for c = 1:numel(cuts) - 1
        first = cuts(c);
        last  = cuts(c+1) - 1;
        % optional further split into pieces of at most SegmentLength
        if isinf(o.SegmentLength)
            starts = first;
            pieceLen = last - first + 1;
        else
            pieceLen = max(1, round(o.SegmentLength * fs));
            starts = first:pieceLen:last;
        end
        for s0 = starts
            s1 = min(last, s0 + pieceLen - 1);
            if s1 - s0 + 1 < 2*fs          % too short for a meaningful fit
                continue
            end
            Y  = double(EEG.data(chans, s0:s1, ep))';     % samples x channels
            mu = mean(Y, 1);
            [Y, segInfo] = removeSinusoids(Y - mu, fs, freqs, o);
            EEG.data(chans, s0:s1, ep) = cast((Y + mu)', 'like', EEG.data);
            for k = 1:numel(segInfo)
                segInfo(k).segment = [s0 s1];
            end
            info = [info, segInfo]; %#ok<AGROW>
        end
    end
end

EEG.icaact = [];
if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc), EEG.etc = struct(); end
EEG.etc.linenoise = info;
cmd = sprintf('EEG = eeg_linenoise(EEG, ''LineFreq'', %g, ''Harmonics'', %s, ''SegmentLength'', %g);', ...
    o.LineFreq, mat2str(o.Harmonics), o.SegmentLength);
EEG.history = [EEG.history newline cmd];

if o.Verbose && ~isempty(info)
    for f = freqs
        k    = [info.nominal] == f;
        done = k & ~[info.skipped];
        if ~any(done)
            fprintf('eeg_linenoise: %g Hz -> no clear peak in any of %d pieces, nothing removed\n', f, sum(k));
            continue
        end
        est = [info(done).freq];
        fprintf('eeg_linenoise: %g Hz -> estimated %.4f Hz (range %.4f-%.4f) in %d of %d pieces, median %.2f%% of channel variance removed\n', ...
            f, median(est), min(est), max(est), sum(done), sum(k), 100*median([info(done).removed]));
    end
end
end


% =============================================================================
function [Y, segInfo] = removeSinusoids(Y, fs, freqs, o)
% Y: samples x channels, mean removed. Removes one sinusoid per frequency.
n   = size(Y, 1);
t   = (0:n-1)' / fs;
tot = sum(Y.^2, 1);                                   % per-channel energy
segInfo = struct('segment', {}, 'nominal', {}, 'freq', {}, 'peakDb', {}, 'skipped', {}, 'removed', {}, 'amplitude', {});

% FFT of the whole stretch, all channels at once; mean power over channels
F   = fft(Y);
P   = mean(abs(F(1:floor(n/2)+1, :)).^2, 2);
fax = (0:floor(n/2))' * fs / n;
df  = fs / n;
clear F

for f0 = freqs
    inb = fax >= f0 - o.SearchHz & fax <= f0 + o.SearchHz;
    if ~any(inb)
        continue
    end
    idx   = find(inb);
    [~, m] = max(P(idx));
    fpk   = fax(idx(m));

    % is there a clear peak? compare with the median power 1.5-5 Hz away
    around = abs(fax - f0) >= 1.5 & abs(fax - f0) <= 5 & fax > 0;
    peakDb = 10*log10(P(idx(m)) / median(P(around)));
    k = numel(segInfo) + 1;
    segInfo(k).nominal = f0;
    segInfo(k).peakDb  = peakDb;
    if ~(peakDb >= o.MinPeakDb)
        segInfo(k).freq      = fpk;
        segInfo(k).skipped   = true;
        segInfo(k).removed   = zeros(1, size(Y, 2));
        segInfo(k).amplitude = zeros(1, size(Y, 2));
        continue
    end

    % grid search around the FFT peak, zooming in
    centre = fpk;
    half   = df;                                      % +/- one FFT bin
    best   = fpk;
    bestE  = -Inf;
    for lev = 1:o.Levels
        grid = centre + linspace(-half, half, o.GridPoints);
        for g = grid
            E = explained(Y, t, g);
            if E > bestE
                bestE = E;
                best  = g;
            end
        end
        centre = best;
        half   = 2*half / (o.GridPoints - 1);         % next level spans one step
    end

    % regress the best cosine/sine pair out of all channels
    X = [cos(2*pi*best*t), sin(2*pi*best*t)];
    B = (X'*X) \ (X'*Y);                              % 2 x channels
    fit = X*B;
    Y = Y - fit;

    segInfo(k).freq      = best;
    segInfo(k).skipped   = false;
    segInfo(k).removed   = sum(fit.^2, 1) ./ max(tot, eps);
    segInfo(k).amplitude = sqrt(sum(B.^2, 1));
end
end

% -----------------------------------------------------------------------------
function E = explained(Y, t, f)
% Variance explained, summed over channels, by a cosine/sine pair at f.
X   = [cos(2*pi*f*t), sin(2*pi*f*t)];
XtY = X' * Y;                                         % 2 x channels
E   = sum(sum(XtY .* ((X'*X) \ XtY), 1));
end
