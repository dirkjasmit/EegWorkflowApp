function data = asr_process_fast(data, srate, state, windowlen, lookahead, stepsize, maxdims, splits)
% ASR_PROCESS_FAST  asr_process of clean_rawdata, without the channels x
% channels x samples covariance series.
%
%   data = asr_process_fast(data, srate, state, windowlen, lookahead, ...
%                           stepsize, maxdims, splits)
%
% asr_process smooths the outer product x*x' over N samples at EVERY sample
% (moving_average over C*C rows) and then uses only every stepsize-th column.
% For 129 channels and 190 s that series is 12.6 GB, which is why clean_asr
% cuts the data into memory blocks; because the update grid restarts at every
% block, its result depends on how much memory it was given.
%
% Here the same covariance is computed where it is used,
%
%   Xcov = X(:, u-N+1 : u) * X(:, u-N+1 : u)' / N
%
% one small matrix product per update position. Everything else is the
% original: the same IIR filter, eigendecomposition, keep rule, blending and
% edge carry.
%
%   splits  1 (default) puts one update grid over the whole recording, which
%           is deterministic; n > 1 restarts the grid every 1/n-th of the
%           data, which is what the original does with n memory blocks. Pass
%           the same n to compare the two.
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

if nargin < 8 || isempty(splits), splits = 1; end
[C, S] = size(data);
N = round(windowlen*srate);
P = round(lookahead*srate);
if maxdims < 1, maxdims = round(C*maxdims); end
T = state.T;
M = state.M;
A = state.A;
B = state.B;

% the same edge padding as the original
carry = repmat(2*data(:,1), 1, P) - data(:, 1 + mod(((P+1):-1:2) - 1, S));
data  = [carry data];
data(~isfinite(data(:))) = 0;

last_R = [];
last_trivial = true;
iir = [];
prev = zeros(C, N);          % the N filtered samples before this block
for i = 1:splits
    range = 1 + floor((i-1)*S/splits) : min(S, floor(i*S/splits));
    if isempty(range), continue, end

    [X, iir] = filter(B, A, double(data(:, range + P)), iir, 2);

    % update positions, exactly as the original lays them out per block
    n = size(X, 2);
    update_at = min(stepsize:stepsize:(n + stepsize - 1), n);
    if isempty(last_R)
        update_at = [1 update_at]; %#ok<AGROW>
        last_R = eye(C);
    end

    % the covariance over the N samples ending at each update position, from
    % the block's own data plus what came before it
    XX = [prev X];
    last_n = 0;
    for j = 1:numel(update_at)
        u    = update_at(j);
        stop = u + size(prev, 2);                        % position inside XX
        w    = XX(:, max(1, stop - N + 1) : stop);
        if size(w, 2) < N                                % at the very start the
            w = [zeros(C, N - size(w, 2)) w]; %#ok<AGROW> % original pads with zeros
        end
        Xcov = (w * w') / N;

        [V, D] = eig(Xcov);
        [D, order] = sort(reshape(diag(D), 1, C));
        V = V(:, order);
        keep = D < sum((T*V).^2) | (1:C) < (C - maxdims);
        trivial = all(keep);
        if ~trivial
            R = real(M * pinv(bsxfun(@times, keep', V'*M)) * V');
        else
            R = eye(C);
        end
        if ~trivial || ~last_trivial
            subrange = range((last_n+1):u);
            blend = (1 - cos(pi*(1:(u - last_n))/(u - last_n)))/2;
            data(:, subrange) = bsxfun(@times, blend, R*data(:, subrange)) + ...
                                bsxfun(@times, 1 - blend, last_R*data(:, subrange));
        end
        last_n = u;
        last_R = R;
        last_trivial = trivial;
    end
    prev = XX(:, end-N+1:end);
end
data = data(:, 1:(end-P));
end
