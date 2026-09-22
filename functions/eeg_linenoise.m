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
%      because line noise reaches each electrode wit its own amplitude and
%      phase.
%
%   4. The mains frequency also drifts WITHIN a piece (HBN: typically 0.017 Hz
%      between neighbouring 5 s pieces). A constant-frequency sinusoid then
%      drifts out of phase towards the ends of the piece, which limits the
%      removal to roughly 30-35 dB. With 'Drift' on (default), the frequency
%      is modelled as changing linearly over the piece: the regressors are the
%      cosine and sine of 2*pi*(f*t + k*t^2/2), t centred on the middle of the
%      piece, so f is the frequency at the centre and k the drift in Hz/s.
%      f and k are found together by the same zooming grid search (now over a
%      2-D grid). k = 0 is on the grid, so the drift model never fits worse
%      than a constant frequency.
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
% phase) and, for epoched data, per epoch. Each stretch is split further into
% equal pieces of about 'SegmentLength' seconds (a 38 s stretch with 5 s
% pieces becomes 8 pieces of 4.75 s, not 7 x 5 s plus a 3 s rest), because
% the mains frequency is not constant: on a 5.4 min recording the estimated
% frequency ranged from 59.99 to 60.02 Hz between 20 s pieces, and a single
% sinusoid over the whole recording removed almost nothing (60 Hz peak
% 3.9 dB -> 3.9 dB). Tracking HBN mains over time (eeg_linedrift) showed slow
% trends of ~0.02 Hz per 20-60 s plus wiggles of ~0.002 Hz over seconds;
% 5 s pieces with the drift fit follow that best. Stretches shorter than 2 s
% are left unchanged (too short for a frequency estimate). Run line-noise
% removal before cutting out periods, so the stretches stay long.
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
%   'Drift'         also fit a linear frequency drift within each piece
%                   (default true); false gives the constant-frequency fit
%   'MaxDrift'      largest drift searched, in Hz per second (default 0.02,
%                   i.e. +/-0.1 Hz over a 5 s piece)
%   'SegmentLength' approximate piece length in seconds; each stretch between
%                   boundaries is split into equal pieces of about this
%                   length. Inf = whole stretch (default 5)
%   'Channels'      channel indices to clean (default all)
%   'Verbose'       print a summary (default true)
%
% OUTPUT
%   EEG             cleaned dataset; EEG.etc.linenoise holds info
%   info            struct array, one element per stretch x harmonic:
%                     .segment    [first last] sample
%                     .nominal    nominal frequency (Hz)
%                     .freq       estimated frequency (Hz) at the middle of the
%                                 piece; FFT peak when skipped
%                     .drift      estimated drift (Hz/s); 0 without 'Drift'
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
addParameter(p, 'Drift', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'MaxDrift', 0.02, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'SegmentLength', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);
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
info  = struct('segment', {}, 'nominal', {}, 'freq', {}, 'drift', {}, 'peakDb', {}, 'skipped', {}, 'removed', {}, 'amplitude', {});
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
        % split into equal pieces of about SegmentLength (no short rest piece)
        len = last - first + 1;
        if isinf(o.SegmentLength)
            nP = 1;
        else
            nP = max(1, round(len / (o.SegmentLength * fs)));
        end
        edges = round(linspace(first, last + 1, nP + 1));
        for ip = 1:nP
            s0 = edges(ip);
            s1 = edges(ip + 1) - 1;
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
cmd = sprintf('EEG = eeg_linenoise(EEG, ''LineFreq'', %g, ''Harmonics'', %s, ''SegmentLength'', %g, ''Drift'', %d);', ...
    o.LineFreq, mat2str(o.Harmonics), o.SegmentLength, logical(o.Drift));
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
t   = ((0:n-1)' - (n-1)/2) / fs;                      % centred: f is the mid-piece frequency
tot = sum(Y.^2, 1);                                   % per-channel energy
segInfo = struct('segment', {}, 'nominal', {}, 'freq', {}, 'drift', {}, 'peakDb', {}, 'skipped', {}, 'removed', {}, 'amplitude', {});

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
        segInfo(k).drift     = 0;
        segInfo(k).skipped   = true;
        segInfo(k).removed   = zeros(1, size(Y, 2));
        segInfo(k).amplitude = zeros(1, size(Y, 2));
        continue
    end

    % grid search around the FFT peak, zooming in: frequency, and with Drift
    % on also the drift rate (2-D grid; drift 0 is always a grid point)
    cF = fpk;  hF = df;                               % +/- one FFT bin
    cK = 0;    hK = o.MaxDrift * logical(o.Drift);
    best  = [fpk 0];
    bestE = -Inf;
    for lev = 1:o.Levels
        gF = cF + linspace(-hF, hF, o.GridPoints);
        if hK > 0
            gK = cK + linspace(-hK, hK, o.GridPoints);
        else
            gK = 0;
        end
        for kk = gK
            for g = gF
                E = explained(Y, t, g, kk);
                if E > bestE
                    bestE = E;
                    best  = [g kk];
                end
            end
        end
        cF = best(1);  hF = 2*hF / (o.GridPoints - 1);    % next level spans one step
        cK = best(2);  hK = 2*hK / (o.GridPoints - 1);
    end

    % regress the best cosine/sine pair out of all channels
    ph = 2*pi*(best(1)*t + 0.5*best(2)*t.^2);
    X = [cos(ph), sin(ph)];
    B = (X'*X) \ (X'*Y);                              % 2 x channels
    fit = X*B;
    Y = Y - fit;

    segInfo(k).freq      = best(1);
    segInfo(k).drift     = best(2);
    segInfo(k).skipped   = false;
    segInfo(k).removed   = sum(fit.^2, 1) ./ max(tot, eps);
    segInfo(k).amplitude = sqrt(sum(B.^2, 1));
end
end

% -----------------------------------------------------------------------------
function E = explained(Y, t, f, k)
% Variance explained, summed over channels, by a cosine/sine pair at
% frequency f (at t = 0) drifting by k Hz/s.
ph  = 2*pi*(f*t + 0.5*k*t.^2);
X   = [cos(ph), sin(ph)];
XtY = X' * Y;                                         % 2 x channels
E   = sum(sum(XtY .* ((X'*X) \ XtY), 1));
end
