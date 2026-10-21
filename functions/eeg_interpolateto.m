function [EEG, info] = eeg_interpolateto(EEG, target, varargin)
% EEG_INTERPOLATETO  Replace the montage by a fixed set of electrode positions,
% every target channel spherically interpolated from the recorded ones, so that
% datasets recorded with different caps become comparable.
%
%   [EEG, info] = eeg_interpolateto(EEG, target)
%   [EEG, info] = eeg_interpolateto(EEG, target, 'Align', 'labels', ...)
%
% target is either a chanlocs struct array, or a cell array of labels that are
% looked up in a template file ('Template', default resources/
% standard-10-5-cap385.elp: a sphere of radius 85 with 10-5 names).
%
% MATCHING TWO MONTAGES
% Spherical splines use only the DIRECTION of each electrode from the centre
% of the head: every position is normalised to the unit sphere first, so the
% radius of a montage (85 mm, 1, or a realistic head with radii of 83-117)
% makes no difference, and there is no need to project onto a sphere of any
% particular size first.
%
% Both montages do have to be in the same coordinate convention, and through
% the app they are: EEGLAB's lookup (pop_chanedit) puts every template in its
% own convention. Measured on the files this app ships: after a lookup with
% standard-10-5-cap385.elp the positions are identical to the template
% (0.00 degrees for every shared label), and after a lookup with EEGLAB's
% realistic head template they differ by 6.1 degrees median (10.7 worst),
% which is the real difference between a realistic head and a sphere.
%
% So no rotation is fitted by default. The source is only recentred on its
% best-fitting sphere (chancenter), which matters for a realistic montage
% whose electrode cloud is not centred on the origin, and is a no-op for a
% spherical one. The shared labels are compared anyway and reported in
% info.mismatch: a large angle means the montage is NOT in EEGLAB's
% convention (raw coordinates from a file, read without a lookup, can be off
% by 80 degrees), and then either run Lookup first or pass 'Align' 'labels',
% which fits a rotation on the shared labels (and mirrors when the montage
% turns out to be left-handed).
%
% OPTIONS (name/value)
%   'Template'   file the target labels are looked up in (default:
%                resources/standard-10-5-cap385.elp next to this function)
%   'Align'      'none' (default) use the coordinates as they are, which is
%                right when both montages went through EEGLAB's lookup;
%                'labels' fit a rotation on the shared labels, for montages
%                that are not in EEGLAB's convention
%   'Recenter'   recentre the source on its best-fitting sphere (default true)
%   'MinShared'  minimum shared labels for the fit (default 4)
%   'Method'     'spherical' (default) or 'sphericalKang', as eeg_interp
%   'Channels'   source channels to interpolate from (default: every located
%                channel whose label does not contain 'eog')
%
% OUTPUT
%   EEG    dataset with the target channels only, in the target order; the
%          target positions are kept as they are in the template
%   info   .labels     target labels
%          .source     source channels used
%          .shared     labels used for the rotation
%          .aligned    true when a rotation was fitted ('Align' 'labels')
%          .mismatch   angle per shared label between the source and the
%                      target positions, before any alignment (degrees):
%                      0 after a lookup with the same template
%          .mirrored   true when the source montage turned out to be
%                      left-handed (X/Y swapped, or left and right reversed)
%                      and had to be mirrored to fit: worth checking, since
%                      it means the recorded positions do not follow EEGLAB's
%                      convention
%          .residual   angle per shared label after alignment (degrees)
%          .exact      target labels that were also recorded (their signal is
%                      still interpolated, so that every channel is treated
%                      the same way)
%
% Old-style labels of the source (T3 T4 T5 T6 A1 A2) are matched to their
% 10-10 names (T7 T8 P7 P8 M1 M2) for the label comparison only.
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

p = inputParser;
p.addParameter('Template', '');
p.addParameter('Align', 'none');
p.addParameter('Recenter', true);
p.addParameter('MinShared', 4);
p.addParameter('Method', 'spherical');
p.addParameter('Channels', []);
p.parse(varargin{:});
o = p.Results;

% ---- target positions --------------------------------------------------------
if iscell(target) || isstring(target)
    tmpl = o.Template;
    if isempty(tmpl)
        tmpl = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'resources', ...
            'standard-10-5-cap385.elp');
    end
    if ~exist(tmpl, 'file')
        error('eeg_interpolateto:NoTemplate', 'Template file %s not found.', tmpl);
    end
    T = readlocs(tmpl);
    want = cellstr(target(:))';
    [found, where] = ismember(canonical(want), canonical({T.labels}));
    if ~all(found)
        error('eeg_interpolateto:UnknownLabels', 'Not in %s: %s', ...
            tmpl, strjoin(want(~found), ' '));
    end
    tgt = T(where);
    for k = 1:numel(tgt)
        tgt(k).labels = want{k};                  % keep the names as asked for
    end
else
    tgt = target(:)';
end
if any(cellfun(@isempty, {tgt.X}))
    error('eeg_interpolateto:TargetNoLocation', 'Every target channel needs a location.');
end

% ---- source channels ---------------------------------------------------------
src = o.Channels;
if isempty(src)
    located = find(~cellfun(@isempty, {EEG.chanlocs.X}));
    iseog   = contains(lower({EEG.chanlocs(located).labels}), 'eog');
    src     = located(~iseog);
end
src = src(:)';
if numel(src) < 4
    error('eeg_interpolateto:TooFewChannels', ...
        'At least 4 located source channels are needed (found %d).', numel(src));
end

% ---- positions on the unit sphere, source recentred --------------------------
S = [[EEG.chanlocs(src).X]; [EEG.chanlocs(src).Y]; [EEG.chanlocs(src).Z]]';
if o.Recenter
    try
        [cx, cy, cz] = chancenter(S(:,1), S(:,2), S(:,3), []);
        S = [cx cy cz];                            % best-fitting sphere centre
    catch
        S = S - mean(S, 1);
    end
end
S = S ./ sqrt(sum(S.^2, 2));
Tg = [[tgt.X]; [tgt.Y]; [tgt.Z]]';
Tg = Tg ./ sqrt(sum(Tg.^2, 2));

% ---- align the source to the target by the labels they share -----------------
info = struct('labels', {{tgt.labels}}, 'source', src, 'shared', {{}}, ...
    'aligned', false, 'mirrored', false, 'residual', [], 'mismatch', [], 'exact', {{}});
srcLab = canonical({EEG.chanlocs(src).labels});
tgtLab = canonical({tgt.labels});
[isShared, pos] = ismember(srcLab, tgtLab);
info.exact = {tgt(pos(isShared)).labels};
if any(isShared)
    % how well do the two montages agree where they share a label? (a check,
    % not a correction: after a lookup this is 0 for the same template)
    info.mismatch = acosd(max(-1, min(1, sum(S(isShared,:) .* Tg(pos(isShared),:), 2))))';
    info.shared   = {EEG.chanlocs(src(isShared)).labels};
end
if strcmpi(o.Align, 'labels') && nnz(isShared) >= o.MinShared
    A = S(isShared, :);
    B = Tg(pos(isShared), :);
    [U, ~, V] = svd(A' * B);
    Rrot = U * diag([1 1 sign(det(U*V'))]) * V';   % rotation only
    Rref = U * diag([1 1 -sign(det(U*V'))]) * V';  % with a mirroring
    fit = @(R) median(acosd(max(-1, min(1, sum((A*R) .* B, 2)))));
    R = Rrot;
    if fit(Rref) < 0.5 * fit(Rrot)
        % a left-handed coordinate file (X and Y swapped, or left/right
        % reversed): mirroring fits far better than any rotation can
        R = Rref;
        info.mirrored = true;
    end
    S = S * R;
    info.aligned  = true;
    info.shared   = {EEG.chanlocs(src(isShared)).labels};
    info.residual = acosd(max(-1, min(1, sum(S(isShared,:) .* B, 2))))';
end

% ---- interpolate: target channels predicted from the source channels ---------
% one dataset holding both montages, so eeg_fastinterpolate does the work
C = EEG;
C.data     = [double(EEG.data(src, :)); zeros(numel(tgt), size(EEG.data(:,:), 2))];
C.nbchan   = numel(src) + numel(tgt);
C.chanlocs = minimalLocs([{EEG.chanlocs(src).labels}, {tgt.labels}], [S; Tg]);
C.trials   = 1;
C.pnts     = size(C.data, 2);
C.icaweights = []; C.icasphere = []; C.icawinv = []; C.icaact = []; C.icachansind = [];
newIdx = numel(src) + (1:numel(tgt));
imp = eeg_fastinterpolate(C, newIdx, o.Method, 1:numel(src));

% ---- the new dataset ---------------------------------------------------------
sz = size(EEG.data);
EEG.data     = cast(reshape(imp, [numel(tgt), sz(2:end)]), 'like', EEG.data);
EEG.chanlocs = tgt;
EEG.nbchan   = numel(tgt);
EEG.icaweights = []; EEG.icasphere = []; EEG.icawinv = []; EEG.icaact = []; EEG.icachansind = [];
if isfield(EEG, 'chaninfo')
    EEG.chaninfo.removedchans = [];
    EEG.chaninfo.nodatchans   = [];
end
if ~isfield(EEG, 'etc') || ~isstruct(EEG.etc), EEG.etc = struct(); end
EEG.etc.interpolateto = info;
EEG.history = [EEG.history newline sprintf( ...
    'EEG = eeg_interpolateto(EEG, {%s}); %% %d channels interpolated from %d recorded ones', ...
    strjoin(cellfun(@(s) ['''' s ''''], {tgt.labels}, 'UniformOutput', false), ' '), ...
    numel(tgt), numel(src))];
end

% -----------------------------------------------------------------------------
function loc = minimalLocs(labels, xyz)
% chanlocs with just the labels and the unit-sphere positions: that is all the
% spline uses, and it avoids mixing the fields of two montages.
loc = struct('labels', labels(:)', 'X', num2cell(xyz(:,1))', ...
    'Y', num2cell(xyz(:,2))', 'Z', num2cell(xyz(:,3))');
loc = convertlocs(loc, 'cart2all');   % the angles too: eeg_fastinterpolate
end                                   % counts a channel as located only with both

% -----------------------------------------------------------------------------
function c = canonical(labels)
% Upper case, and the old 10-20 names mapped to their 10-10 equivalents, for
% comparing labels between montages.
old = {'T3','T7'; 'T4','T8'; 'T5','P7'; 'T6','P8'; 'A1','M1'; 'A2','M2'};
c = upper(strtrim(cellstr(labels)));
for k = 1:size(old, 1)
    c(strcmp(c, old{k,1})) = old(k,2);
end
c = c(:)';
end
