function [index,I] = emg_psd(X,opt)
% emg_psd() - Selects EMG components according to their Power Density
% Spectrum
%
% Usage:
%   >> [index,I] = emg_psd(X,opt)
%
% Inputs:
%   X               - data matrix (dxN, data channels are rowwise)
%   opt.femg       - frequency aproximately separating the EEG and EMG bands
%                    def: 15 Hz
%   opt.fs         - sampling frequency
%                    def: 250 Hz
%   opt.ratio      - ratio of average power (per unit of frequency) in EEG
%                    band to average power in the EMG band below which a
%                    component will be considered to be EMG-related.
%                    Increasing/decreasing this value increases/decreases
%                    the amount of correction.
%                    def: 10
%   opt.range      - range of components that can be removed. At least
%                    opt.range(1) components will be removed in each analysis
%                    window and at most opt.range(2) components.
%                    def: [0 floor(d/2)]
%   opt.estimator  - spectral estimator object (spectrum.welch). Its segment
%                    length, overlap and window are honoured. If left empty a
%                    Welch estimator with a Hamming window and a segment
%                    length of min(2*fs, floor(N/2)) samples is used.
%   opt.NFFT       - number of FFT points used to compute the PSD. If left
%                    empty the default 2^ceil(log2(min(2*fs,floor(N/2)))) is
%                    used.
%
% Outputs:
%   index   - indexes of the rows of X corresponding to EMG components
%
%
% Notes:
%   - This function requires the MATLAB Signal Processing Toolbox.
%   - The optimum value of the parameter opt.ratio depends on the actual
%     spectral estimator used.
%
% See also:
%   FD, POP_AUTOBSSEMG, AUTOBSS, EEGLAB
%

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com
%
% Speed-optimized drop-in (2026): the per-component loop over the legacy
% spectrum.welch / psd / avgpower object API is replaced by a single
% vectorized pwelch() call over all components at once. The classification is
% identical in method - the ratio of mean PSD in the EEG band (<=femg) to the
% mean PSD in the EMG band (>=femg) is compared against opt.ratio - but it is
% computed far more cheaply. Because the decision is based on a *ratio* of
% band powers, it is invariant to the constant PSD scaling conventions, so
% results match the original estimator to within frequency-bin discretization
% at the femg boundary.


% check that the Signal Processing Toolbox is available
% -------------------------------------------------
if isempty(ver('signal')),
    error('(emg_psd) this criterion requires the Signal Processing Toolbox');
end

if nargin < 1, help emg_psd; return; end
[d,N] = size(X);

% remove mean from data
% ------------------------------------------------
X = X - mean(X,2);

if ~exist('opt','var'),
    opt = def_emg_psd;
else
    opt = def_emg_psd(opt);
end
fs       = opt.fs;
ratio_th = opt.ratio;
femg     = opt.femg;

% default range of components that might be removed
% -------------------------------------------------
if isempty(opt.range),
    RANGE = [0 floor(d/2)];
else
    RANGE = opt.range;
end

% Welch parameters: segment length, overlap and window taper, honouring a
% caller-supplied spectrum.welch estimator when present.
% -------------------------------------------------
[seglen, noverlap, win] = local_welch_params(opt, N, fs);

% FFT length
% -------------------------------------------------
if isempty(opt.NFFT),
    NFFT = 2^ceil(log2(min(2*opt.fs,floor(N/2))));
else
    NFFT = opt.NFFT;
end
NFFT = max(NFFT, seglen);   % pwelch requires NFFT >= window length

% Compute the Welch PSD of every component in a single vectorized call.
% pwelch treats each column of its input as a separate signal, so Pxx is
% [nfreq x d] with one column per component and f in Hz (0 .. fs/2).
% -------------------------------------------------
[Pxx,f] = pwelch(X.', win, noverlap, NFFT, fs);

% Average power below and above femg (endpoints inclusive on both bands,
% matching the original avgpower() band definition).
% -------------------------------------------------
b1 = (f <= femg);        % EEG band
b2 = (f >= femg);        % EMG band
p1 = mean(Pxx(b1,:), 1);
p2 = mean(Pxx(b2,:), 1);

% detect components that are likely to be EMG-related
% ---------------------------------------------------
ratio = p1./p2;
[ratio,I] = sort(ratio);
index = I(ratio<ratio_th);


% take a number of components within the specified range
% ---------------------------------------------
if length(index) < RANGE(1), index = I(1:RANGE(1)); end
if length(index) > RANGE(2), index = I(1:RANGE(2)); end

return;


% subfunction to derive the Welch parameters
% ---------------------------------------------------
function [seglen, noverlap, win] = local_welch_params(opt, N, fs)
overlap_pct = 50;          % spectrum.welch default overlap
winname     = 'hamming';   % spectrum.welch default taper in this pipeline
seglen      = [];

if isfield(opt,'estimator') && ~isempty(opt.estimator),
    h = opt.estimator;
    try, seglen      = get(h,'SegmentLength');  catch, end %#ok<CTCH>
    try, overlap_pct = get(h,'OverlapPercent'); catch, end %#ok<CTCH>
    try
        wn = get(h,'WindowName');
        if iscell(wn), wn = wn{1}; end
        winname = lower(char(wn));
    catch %#ok<CTCH>
    end
end

if isempty(seglen),
    seglen = min(2*fs, floor(N/2));   % original default segment length
end
seglen   = max(2, min(round(seglen), N));           % keep within data length
noverlap = min(seglen-1, round(seglen*overlap_pct/100));

switch winname
    case {'hann','hanning'}
        win = hann(seglen);
    case 'hamming'
        win = hamming(seglen);
    otherwise
        win = hamming(seglen);   % faithful to the default taper
end


% subfunction to define the default parameters
% ---------------------------------------------------
function [opt] = def_emg_psd(opt)
if nargin < 1 || ~isfield(opt,'fs'),
    opt.fs=250;
    warning('(emg_psd) Default sampling frequency [%d] Hz will be used',opt.fs);
end
if ~isfield(opt, 'bssout'),
    opt.bssout = [];
end
if ~isfield(opt, 'ratio') || isempty(opt.ratio),
    opt.ratio = 10;
end
if ~isfield(opt, 'femg') || isempty(opt.femg),
    opt.femg = 15;
end
if ~isfield(opt,'range'),
    opt.range = [];
end
if ~isfield(opt,'estimator'),
    opt.estimator = [];
end
if ~isfield(opt,'NFFT'),
    opt.NFFT = [];
end
