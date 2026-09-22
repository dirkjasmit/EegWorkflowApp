function [EEG, out] = eeg_autoreject(EEG, varargin)
% EEG_AUTOREJECT  Local autoreject (Jas et al. 2017) on continuous data cut
% into non-overlapping epochs: per-channel peak-to-peak thresholds learned by
% cross-validation, bad channels interpolated within an epoch, and epochs
% with too many bad channels marked for removal.
%
%   [EEG, out] = eeg_autoreject(EEG, 'Epochs', [start stop], ...)
%
% 1. Per channel and epoch, the peak-to-peak (ptp) amplitude.
% 2. Per channel, a threshold: a cell (channel x epoch) is bad when its ptp
%    exceeds it. The threshold is chosen by K-fold cross-validation over all
%    observed ptp values of that channel: for each candidate, the mean of the
%    training epochs under the threshold is compared with the median of the
%    validation epochs (RMSE over time); the candidate with the lowest mean
%    loss over folds wins.
% 3. Per epoch, with n bad cells: when n > kappa * nChannels the epoch is
%    dropped; otherwise its min(n, rho) worst bad channels (largest ptp) are
%    interpolated from the other channels (spherical splines, as eeg_interp).
%    kappa and rho are chosen together by the same cross-validation, now on
%    all channels: mean of the repaired, kept training epochs against the
%    median of the raw validation epochs.
%
% Epochs are baseline corrected (epoch mean removed per channel) for the loss,
% as MNE epochs normally are; the ptp does not depend on it. Folds are
% contiguous blocks of epochs (scikit-learn KFold without shuffling, as
% autoreject), so the result is deterministic. The thresholds are searched
% exhaustively over the observed ptp values instead of by autoreject's
% Bayesian optimisation; with contiguous folds this finds the same optimum
% up to the spacing of those values.
%
% OPTIONS (name/value)
%   'Epochs'     [start stop] sample pairs, non-overlapping (required)
%   'Channels'   channels tested and repaired; they need a location
%                (default: all located channels). Other channels are left
%                as they are.
%   'Folds'      number of cross-validation folds (default 10)
%   'Kappa'      candidate kappa values (default 0:0.1:1)
%   'Rho'        candidate rho values (default [1 4 32], limited to fewer
%                than the number of channels)
%
% OUTPUT
%   EEG          data with the repaired cells replaced by their
%                interpolation. Dropped epochs are NOT removed here.
%   out          .epochs      [start stop] per epoch
%                .chans       the tested channels
%                .thresholds  per channel (data units)
%                .ptp         channels x epochs peak-to-peak
%                .bad         channels x epochs, ptp above threshold
%                .repaired    channels x epochs, interpolated cells
%                .dropped     1 x epochs, epochs to remove
%                .kappa, .rho the chosen values; .loss the kappa x rho
%                             cross-validation loss
%
% Reference: Jas M, Engemann DA, Bekhti Y, Raimondo F, Gramfort A (2017).
% Autoreject: Automated artifact rejection for MEG and EEG data. NeuroImage
% 159:417-429. doi:10.1016/j.neuroimage.2017.06.030. Reimplemented in MATLAB,
% not a port of the Python package.
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

p = inputParser;
p.addParameter('Epochs', []);
p.addParameter('Channels', []);
p.addParameter('Folds', 10);
p.addParameter('Kappa', 0:0.1:1);
p.addParameter('Rho', [1 4 32]);
p.parse(varargin{:});
o = p.Results;
if EEG.trials > 1
    error('eeg_autoreject:Epoched', 'eeg_autoreject needs continuous data.');
end
if isempty(o.Epochs)
    error('eeg_autoreject:NoEpochs', 'Give the epochs as [start stop] sample pairs.');
end
chans = o.Channels;
if isempty(chans)
    chans = find(~cellfun(@isempty, {EEG.chanlocs.X}));
end
chans = chans(:)';
if any(cellfun(@isempty, {EEG.chanlocs(chans).X}))
    error('eeg_autoreject:NoLocation', 'Every channel in Channels needs a location.');
end
ep    = round(o.Epochs);
nEp   = size(ep, 1);
L     = ep(1,2) - ep(1,1) + 1;
if any(ep(:,2) - ep(:,1) + 1 ~= L)
    error('eeg_autoreject:EpochLength', 'All epochs must have the same length.');
end
nCh   = numel(chans);
K     = min(o.Folds, nEp);
if K < 2 || nEp < 4
    error('eeg_autoreject:TooFewEpochs', 'Too few epochs (%d) for cross-validation.', nEp);
end
kappas = o.Kappa(:)';
rhos   = unique(o.Rho(o.Rho >= 1 & o.Rho < nCh));
if isempty(rhos)
    rhos = 1;
end

% epochs: channels x time x epochs, baseline corrected
X = zeros(nCh, L, nEp);
for e = 1:nEp
    X(:, :, e) = double(EEG.data(chans, ep(e,1):ep(e,2)));
end
ptp = squeeze(max(X, [], 2) - min(X, [], 2));          % channels x epochs
if nCh == 1, ptp = ptp(:)'; end
X = X - mean(X, 2);

% contiguous folds
fold = zeros(1, nEp);
edges = round(linspace(0, nEp, K + 1));
for k = 1:K
    fold(edges(k)+1:edges(k+1)) = k;
end
% validation medians (all channels), once per fold
med = zeros(nCh, L, K);
for k = 1:K
    med(:, :, k) = median(X(:, :, fold == k), 3);
end

% ---- 1. per-channel thresholds -------------------------------------------
thr = zeros(1, nCh);
for c = 1:nCh
    Xc   = squeeze(X(c, :, :));                           % time x epochs
    pc   = ptp(c, :);
    cand = unique(pc);                                    % candidate thresholds
    loss = zeros(K, numel(cand));
    for k = 1:K
        tr = find(fold ~= k);
        [ps, ord] = sort(pc(tr));
        cm = cumsum(Xc(:, tr(ord)), 2) ./ (1:numel(tr));  % mean of the i smallest
        nKeep = sum(ps(:) <= cand, 1);                    % per candidate
        err = sqrt(mean((cm - med(c, :, k)').^2, 1));     % RMSE per i
        l = inf(1, numel(cand));
        ok = nKeep > 0;
        l(ok) = err(nKeep(ok));
        loss(k, :) = l;
    end
    ml = mean(loss, 1);
    best = find(ml == min(ml));
    thr(c) = cand(best(end));                             % ties: keep more data
end
bad  = ptp > thr(:);
nBad = sum(bad, 1);

% ---- 2. kappa and rho ------------------------------------------------------
xyz = [[EEG.chanlocs(chans).X]; [EEG.chanlocs(chans).Y]; [EEG.chanlocs(chans).Z]]';
xyz = xyz ./ sqrt(sum(xyz.^2, 2));
G   = splineG(xyz);
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');

% epoch e is kept for kappa index j when nBad(e) <= kappas(j)*nCh
keepIdx = zeros(1, nEp);                                  % first kappa index keeping e
for e = 1:nEp
    j = find(nBad(e) <= kappas * nCh, 1);
    if isempty(j), j = numel(kappas) + 1; end
    keepIdx(e) = j;
end
nK = numel(kappas);
lossKR = zeros(nK, numel(rhos));
for r = 1:numel(rhos)
    % sums of repaired epochs per fold and per first-keeping kappa index
    S = zeros(nCh, L, K, nK);
    N = zeros(K, nK);
    for e = 1:nEp
        j = keepIdx(e);
        if j > nK, continue; end
        Y = repairEpoch(X(:, :, e), ptp(:, e), bad(:, e), rhos(r), G, cache);
        S(:, :, fold(e), j) = S(:, :, fold(e), j) + Y;
        N(fold(e), j) = N(fold(e), j) + 1;
    end
    S = cumsum(S, 4);  N = cumsum(N, 2);                  % kept for kappa <= j
    St = sum(S, 3);    Nt = sum(N, 1);
    for j = 1:nK
        l = zeros(1, K);
        for k = 1:K
            n = Nt(j) - N(k, j);
            if n < 1
                l(k) = Inf;
                continue
            end
            m = (St(:, :, 1, j) - S(:, :, k, j)) / n;
            l(k) = sqrt(mean((m - med(:, :, k)).^2, 'all'));
        end
        lossKR(j, r) = mean(l);
    end
end
[~, i] = min(lossKR(:));
[jBest, rBest] = ind2sub(size(lossKR), i);
kappa = kappas(jBest);
rho   = rhos(rBest);

% ---- 3. apply ------------------------------------------------------------------
dropped  = nBad > kappa * nCh;
repaired = false(nCh, nEp);
for e = find(~dropped & nBad > 0)
    raw = double(EEG.data(chans, ep(e,1):ep(e,2)));
    [Y, fix] = repairEpoch(raw, ptp(:, e), bad(:, e), rho, G, cache);
    EEG.data(chans, ep(e,1):ep(e,2)) = cast(Y, 'like', EEG.data);
    repaired(fix, e) = true;
end
if any(repaired(:))
    EEG.icaact = [];
end

out = struct('epochs', ep, 'chans', chans, 'thresholds', thr, 'ptp', ptp, ...
    'bad', bad, 'repaired', repaired, 'dropped', dropped, ...
    'kappa', kappa, 'rho', rho, 'kappas', kappas, 'rhos', rhos, 'loss', lossKR);
end

% -----------------------------------------------------------------------------
function [Y, fix] = repairEpoch(Y, ptpE, badE, rho, G, cache)
% Interpolate the rho worst bad channels (largest ptp) from all others.
idx = find(badE);
fix = [];
if isempty(idx)
    return
end
[~, o] = sort(ptpE(idx), 'descend');
fix  = sort(idx(o(1:min(rho, numel(idx)))));
key  = sprintf('%d,', fix);
if isKey(cache, key)
    W = cache(key);
else
    good = setdiff(1:size(G, 1), fix);
    W = splineWeights(G, good, fix);
    cache(key) = W;                               % containers.Map is a handle
end
good = setdiff(1:size(G, 1), fix);
Y(fix, :) = W * Y(good, :);
end

% -----------------------------------------------------------------------------
function W = splineWeights(G, good, bad)
% eeg_interp's spherical spline (lambda 0), data centred on the mean of the
% good channels, written as one linear map from good to bad channels.
ng = numel(good);
M  = [G(good, good); ones(1, ng)];
w  = G(bad, good) * pinv(M);
w  = w(:, 1:end-1);
c  = sum(w, 2) - 1;
W  = w - c / ng;
end

% -----------------------------------------------------------------------------
function g = splineG(xyz)
% computeg of eeg_interp (m = 4, 7 Legendre terms) for all pairs of unit vectors.
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
