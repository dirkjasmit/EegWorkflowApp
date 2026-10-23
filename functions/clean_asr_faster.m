function [signal, state] = clean_asr_faster(signal, cutoff, windowlen, stepsize, maxdims, ...
    ref_maxbadchannels, ref_tolerances, ref_wndlen, splits, maxmem)
% CLEAN_ASR_FASTER  ASR (clean_asr) without the channels x channels x samples
% covariance series, so it needs no memory blocks and gives the same answer
% every time.
%
%   signal = clean_asr_faster(signal, cutoff)
%   signal = clean_asr_faster(signal, cutoff, windowlen, stepsize, maxdims, ...
%                             ref_maxbadchannels, ref_tolerances, ref_wndlen, splits)
%
% WHAT IS DIFFERENT
% asr_process smooths the outer product x*x' over a window of N samples at
% EVERY sample (moving_average over C*C rows), then uses only every
% stepsize-th column of it. For 129 channels and 190 s that series is 12.6 GB,
% which is why clean_asr cuts the data into memory blocks - and because ASR
% restarts its update grid at every block, the cleaned data depend on how much
% memory it was given.
%
% The same covariance is obtained here directly where it is used:
%
%   Xcov(:,:,j) = X(:, u_j-N+1 : u_j) * X(:, u_j-N+1 : u_j)' / N
%
% one small matrix product per update position (BLAS, threaded) instead of a
% cumulative sum over C*C*S elements. Nothing else changes: the same IIR
% filter, the same eigendecomposition and keep rule, the same blending
% between updates, the same carry at the edges.
%
%   'maxmem'  memory in MB for the CALIBRATION only (default 64, as
%             clean_asr): asr_calibrate averages its covariance in blocks
%             whose size follows from it, so the same value must be used to
%             compare with clean_asr. The processing step needs no memory
%             setting at all.
%   'splits'  how the update grid is laid out (default 1)
%             1   : one grid over the whole recording - deterministic, and
%                   what this function is for
%             n>1 : restart the grid every 1/n-th of the data, which is what
%                   clean_asr does with n memory blocks. Use it to compare
%                   with the original: with the same n the two agree to
%                   rounding.
%
% The other arguments are those of clean_asr; calibration still runs through
% the plugin's asr_calibrate, which is a small part of the time.
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

if ~exist('cutoff', 'var') || isempty(cutoff), cutoff = 20; end
if ~exist('windowlen', 'var') || isempty(windowlen)
    windowlen = max(0.5, 1.5*signal.nbchan/signal.srate);
end
if ~exist('stepsize', 'var'), stepsize = []; end
if ~exist('maxdims', 'var') || isempty(maxdims), maxdims = 0.66; end
if ~exist('ref_maxbadchannels', 'var') || isempty(ref_maxbadchannels), ref_maxbadchannels = 0.075; end
if ~exist('ref_tolerances', 'var') || isempty(ref_tolerances), ref_tolerances = [-3.5 5.5]; end
if ~exist('ref_wndlen', 'var') || isempty(ref_wndlen), ref_wndlen = 1; end
if ~exist('splits', 'var') || isempty(splits), splits = 1; end
% asr_calibrate derives its block size from maxmem, so the calibration state
% itself depends on it: pass the same value clean_asr would (its default 64)
if ~exist('maxmem', 'var') || isempty(maxmem), maxmem = 64; end
if isempty(stepsize), stepsize = floor(signal.srate*windowlen/2); end

signal.data = double(signal.data);

% ---- calibration, as clean_asr does it --------------------------------------
if isnumeric(ref_maxbadchannels) && isnumeric(ref_tolerances) && isnumeric(ref_wndlen)
    ref_section = clean_windows(signal, ref_maxbadchannels, ref_tolerances, ref_wndlen);
else
    ref_section = signal;
end
% exactly as clean_asr calls it, so the state is the same (its disk cache
% calls asr_calibrate without maxmem, which is why that must match too)
if exist('hlp_diskcache', 'file')
    state = hlp_diskcache('filterdesign', @asr_calibrate, ref_section.data, ref_section.srate, cutoff);
else
    state = asr_calibrate(ref_section.data, ref_section.srate, cutoff, [], [], [], [], [], [], [], maxmem);
end

% ---- processing --------------------------------------------------------------
% clean_asr appends a mirrored tail so that the lookahead of asr_process does
% not shift the output; without it the result comes out P samples early
sig = [signal.data bsxfun(@minus, 2*signal.data(:,end), ...
    signal.data(:, (end-1):-1:end-round(windowlen/2*signal.srate)))];
out = asr_process_fast(sig, signal.srate, state, windowlen, ...
    windowlen/2, stepsize, maxdims, splits);
% clean_asr then drops the leading carry columns, which is what lines the
% result up with the input again
out(:, 1:round(windowlen/2*signal.srate)) = [];
signal.data = out;
signal.etc.clean_asr_faster = struct('cutoff', cutoff, 'windowlen', windowlen, ...
    'stepsize', stepsize, 'splits', splits);
end
