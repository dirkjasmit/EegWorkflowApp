function track = eeg_linedrift(EEG, varargin)
% EEG_LINEDRIFT  Plot how the mains (line-noise) frequency and amplitude
% change over the recording.
%
%   eeg_linedrift(EEG)
%   track = eeg_linedrift(EEG, 'LineFreq', 60, ...)
%
% In short sliding windows, the frequency near LineFreq that explains most
% variance (summed over channels) is found with the same exact cosine/sine
% fit eeg_linenoise uses, and the amplitude of that sinusoid is taken per
% channel. The top panel shows the frequency over time, the bottom panel the
% amplitude (median and maximum over channels). Boundary events are drawn as
% vertical lines; windows that cross one are skipped.
%
% When EEG.etc.linenoise exists (the data went through eeg_linenoise), its
% fits are drawn over the top panel: per piece the fitted frequency at the
% centre, with the fitted drift as the slope of the line. Run it on the data
% BEFORE line-noise removal to see the mains itself; on cleaned data the
% sliding estimate follows whatever is left.
%
% OPTIONS (name/value)
%   'LineFreq'   nominal line frequency in Hz (default 60)
%   'Window'     window length in seconds (default 1)
%   'Step'       step between windows in seconds (default 0.25)
%   'SearchHz'   frequency searched within +/- this of LineFreq (default 0.3)
%   'Channels'   channels to use (default all)
%   'Title'      figure title (default EEG.setname or filename)
%   'Plot'       draw the figure (default true)
%
% OUTPUT
%   track        struct: .time (s, window centres), .freq (Hz),
%                .ampMedian and .ampMax (over channels, data units),
%                .amp (channels x windows), .peakDb (height of the fitted
%                sinusoid's power above the rest of the window, a rough
%                indication of how clear the line is)
%
% EXAMPLE
%   eeg_linedrift(EEG, 'LineFreq', 60);             % before Line noise
%   EEG2 = eeg_linenoise(EEG, 'LineFreq', 60);
%   eeg_linedrift(EEG2, 'LineFreq', 60);            % what is left, plus the fits
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

p = inputParser;
p.addParameter('LineFreq', 60);
p.addParameter('Window', 1);
p.addParameter('Step', 0.25);
p.addParameter('SearchHz', 0.3);
p.addParameter('Channels', []);
p.addParameter('Title', '');
p.addParameter('Plot', true);
p.parse(varargin{:});
o = p.Results;
if EEG.trials > 1
    error('eeg_linedrift:Epoched', 'eeg_linedrift needs continuous data.');
end
fs    = EEG.srate;
chans = o.Channels;
if isempty(chans)
    chans = 1:EEG.nbchan;
end
L    = round(o.Window * fs);
step = max(1, round(o.Step * fs));

% boundaries (windows may not cross them)
bnd = [];
if isfield(EEG, 'event') && ~isempty(EEG.event)
    isB = strcmpi({EEG.event.type}, 'boundary');
    bnd = [EEG.event(isB).latency];
    bnd = bnd(bnd > 1 & bnd < EEG.pnts);
end

starts = 1:step:(EEG.pnts - L + 1);
if ~isempty(bnd)
    crosses = any(bnd(:)' > starts(:) & bnd(:)' < starts(:) + L - 1, 2);
    starts  = starts(~crosses);
end
nW  = numel(starts);
tt  = ((0:L-1)' - (L-1)/2) / fs;
opt = optimset('TolX', 1e-6, 'Display', 'off');

freq   = nan(1, nW);
amp    = nan(numel(chans), nW);
peakDb = nan(1, nW);
for w = 1:nW
    Y = double(EEG.data(chans, starts(w):starts(w) + L - 1))';
    Y = Y - mean(Y, 1);
    f = fminbnd(@(f) -explained(Y, tt, f), o.LineFreq - o.SearchHz, o.LineFreq + o.SearchHz, opt);
    X = [cos(2*pi*f*tt), sin(2*pi*f*tt)];
    B = X \ Y;
    fitE = sum(sum((X*B).^2));
    freq(w)   = f;
    amp(:, w) = sqrt(sum(B.^2, 1))';
    peakDb(w) = 10*log10(fitE / max(sum(Y(:).^2) - fitE, eps) * (L/2));
end

track = struct('time', (starts + (L-1)/2 - 1) / fs, 'freq', freq, ...
    'ampMedian', median(amp, 1), 'ampMax', max(amp, [], 1), 'amp', amp, 'peakDb', peakDb);

if ~o.Plot
    return
end
ttl = o.Title;
if isempty(ttl)
    ttl = EEG.setname;
    if isempty(ttl) && isfield(EEG, 'filename'), ttl = EEG.filename; end
end
fig = figure('Name', sprintf('Line-noise drift %s', ttl), 'Color', 'w');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
title(tl, sprintf('%s  (%g Hz mains, %g s windows)', ttl, o.LineFreq, o.Window), 'Interpreter', 'none');

ax1 = nexttile(tl);
hold(ax1, 'on');
scatter(ax1, track.time, track.freq, 8, track.peakDb, 'filled');
cb = colorbar(ax1);  cb.Label.String = 'line clarity (dB)';
yline(ax1, o.LineFreq, ':', 'Color', [.5 .5 .5]);
if isfield(EEG, 'etc') && isfield(EEG.etc, 'linenoise') && ~isempty(EEG.etc.linenoise)
    Ln = EEG.etc.linenoise;
    Ln = Ln([Ln.nominal] == o.LineFreq & ~[Ln.skipped]);
    for k = 1:numel(Ln)
        seg = Ln(k).segment;
        ts  = [seg(1) seg(end)] / fs;
        tc  = mean(ts);
        dr  = 0;
        if isfield(Ln, 'drift') && ~isempty(Ln(k).drift), dr = Ln(k).drift; end
        h = plot(ax1, ts, Ln(k).freq + dr * (ts - tc), 'r-', 'LineWidth', 1.5);
    end
    if ~isempty(Ln)
        legend(ax1, h, 'eeg\_linenoise fit per piece', 'Location', 'best');
    end
end
ylabel(ax1, 'frequency (Hz)');
grid(ax1, 'on');

ax2 = nexttile(tl);
hold(ax2, 'on');
plot(ax2, track.time, track.ampMax, '-', 'Color', [.85 .45 .1]);
plot(ax2, track.time, track.ampMedian, '-', 'Color', [.1 .3 .7], 'LineWidth', 1.2);
set(ax2, 'YScale', 'log');
legend(ax2, {'max over channels', 'median over channels'}, 'Location', 'best');
ylabel(ax2, 'amplitude');
xlabel(ax2, 'time (s)');
grid(ax2, 'on');

for ax = [ax1 ax2]
    for b = bnd
        xline(ax, b / fs, '-', 'Color', [.6 .6 .6]);
    end
end
linkaxes([ax1 ax2], 'x');
xlim(ax1, [0 EEG.pnts / fs]);
end

% -----------------------------------------------------------------------------
function E = explained(Y, t, f)
X   = [cos(2*pi*f*t), sin(2*pi*f*t)];
XtY = X' * Y;
E   = sum(sum(XtY .* ((X'*X) \ XtY)));
end
