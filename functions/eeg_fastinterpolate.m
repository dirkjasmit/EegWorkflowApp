function [imp, W] = eeg_fastinterpolate(EEG, chans, method)
% EEG_FASTINTERPOLATE  Every channel predicted from all other channels by
% spherical spline interpolation, in one matrix operation.
%
%   imp      = eeg_fastinterpolate(EEG)
%   imp      = eeg_fastinterpolate(EEG, chans)
%   [imp, W] = eeg_fastinterpolate(EEG, chans, method)
%
% For each channel in chans, the signal is interpolated from all OTHER
% channels that have a location, exactly as
%
%   tmp = pop_interp(EEG, k, 'spherical');  imp(k,:) = tmp.data(k,:);
%
% would do, but without running pop_interp once per channel. Spherical spline
% interpolation (eeg_interp) is linear in the data, so the leave-one-out
% predictions of all channels are a single weight matrix W, computed once from
% the electrode positions, times the data. The Legendre series (the expensive
% part) is evaluated once for all electrode pairs, the per-channel spline
% systems are small (channels x channels), and the data are touched in one
% matrix multiply.
%
% INPUT
%   EEG      EEGLAB dataset (continuous or epoched)
%   chans    channel indices to predict (default: all channels). Channels
%            without a location are returned unchanged, as pop_interp does.
%   method   'spherical' (default: lambda 0, m 4, n 7, Perrin et al. 1989) or
%            'sphericalKang' (lambda 1e-8, m 3, n 50, Kang et al. 2015), as in
%            eeg_interp
%
% OUTPUT
%   imp      numel(chans) x pnts (x trials) interpolated data, class of EEG.data
%   W        numel(chans) x nbchan weights: imp = W * EEG.data (per sample)
%
% The spline math follows EEGLAB's eeg_interp (spheric_spline, computeg):
% unit-sphere positions, the data centred on the mean over the predicting
% channels at every sample, and pinv of the spline system with the constraint
% row. Results match the pop_interp loop to rounding error.
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

if nargin < 2 || isempty(chans)
    chans = 1:EEG.nbchan;
end
if nargin < 3 || isempty(method)
    method = 'spherical';
end
switch lower(method)
    case 'spherical',     params = [0 4 7];
    case 'sphericalkang', params = [1e-8 3 50];
    otherwise, error('eeg_fastinterpolate:method', 'Unknown method ''%s''.', method);
end
chans = chans(:)';

% channels with a location, on the unit sphere
hasLoc = ~cellfun(@isempty, {EEG.chanlocs.theta}) & ~cellfun(@isempty, {EEG.chanlocs.X});
L = find(hasLoc);
n = numel(L);
xyz = [[EEG.chanlocs(L).X]; [EEG.chanlocs(L).Y]; [EEG.chanlocs(L).Z]]';
xyz = xyz ./ sqrt(sum(xyz.^2, 2));

% Legendre series for all located pairs, once
G = splineG(xyz, params);

% leave-one-out weights. For target k, predicting from the other located
% channels g (mean over g removed per sample, then added back):
%   imp_k = w*(X_g - m) + m,  m = mean(X_g)  ->  imp_k = w*X_g - (sum(w)-1)*m
% and m = (sum over all located - x_k)/(n-1), which gives one row of W.
W = zeros(numel(chans), EEG.nbchan);
lambda = params(1);
for i = 1:numel(chans)
    k = chans(i);
    pos = find(L == k, 1);
    if isempty(pos) || n < 2
        W(i, k) = 1;                                  % no location: unchanged
        continue
    end
    g  = [1:pos-1, pos+1:n];
    M  = [G(g,g) + lambda*eye(n-1); ones(1, n-1)];
    w  = G(pos, g) * pinv(M);
    w  = w(1:end-1);                                  % constraint column
    c  = sum(w) - 1;
    row = zeros(1, n);
    row(g) = w - c/(n-1);
    W(i, L) = row;                                     % own column stays 0
end

% apply to all samples at once
sz  = size(EEG.data);
X   = reshape(double(EEG.data), sz(1), []);
imp = W * X;
imp = reshape(imp, [numel(chans), sz(2:end)]);
imp = cast(imp, 'like', EEG.data);
end

% -----------------------------------------------------------------------------
function g = splineG(xyz, params)
% computeg of eeg_interp for all pairs of unit vectors (rows of xyz).
EI = 1 - sqrt(max(0, sum(xyz.^2, 2) + sum(xyz.^2, 2)' - 2*(xyz*xyz')));
m    = params(2);
maxn = params(3);
g    = zeros(size(EI));
Pnm1 = ones(size(EI));
Pn   = EI;
for nn = 1:maxn
    g = g + ((2*nn+1)/(nn^m*(nn+1)^m)) * Pn;
    Pnp1 = ((2*nn+1).*EI.*Pn - nn.*Pnm1) / (nn+1);
    Pnm1 = Pn;
    Pn   = Pnp1;
end
g = g/(4*pi);
end
