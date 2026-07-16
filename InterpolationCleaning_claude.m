function [EEGOut, Info] = InterpolationCleaning_claude(EEGIn, varargin)

% InterpolationCleaning_claude  Automated EEG cleaning by interpolation.
%
% For every channel a surrogate signal is estimated from ALL other
% channels (spherical interpolation, see InterpolationReplace). The
% original and the interpolated data are cut into (optionally overlapping)
% epochs and, per channel/epoch, the amount of difference in activity
% between the two is quantified. Where that difference exceeds a threshold
% the channel-by-epoch combination is marked bad.
%
% A channel with more than ChanRejPct(1) (e.g. 30%) bad epochs is marked
% for deletion. Any epoch that still contains a bad channel/epoch mark
% (ignoring the channels that will be deleted) is marked for deletion. If
% too little data would remain (less than MinRemain of the recording), the
% channel-rejection percentage is lowered to the next value in ChanRejPct
% (e.g. 20%) and the channels and epochs to reject are recomputed. Once the
% final set of bad channels and epochs is established they are removed:
% channels with pop_select, epochs (as continuous time ranges) with
% eeg_eegrej.
%
% Requires channel location information in the EEG struct and CONTINUOUS
% data (EEGIn.trials == 1).
%
% USAGE
%   [EEGOut, Info] = InterpolationCleaning_claude(EEG, ...
%                        'ChanRejPct',[.30 .20], 'EpochSec',2, ...
%                        'Overlap',0.5, 'Plot',true);
%
% INPUT
%   EEGIn       EEGLAB EEG struct (continuous). Positional, required.
%
% NAME-VALUE PARAMETERS
%   'ChanRejPct' vector of channel-rejection thresholds to try in order,
%                each the max fraction of bad epochs a channel may have
%                before it is deleted (default [.30 .20]).
%   'EpochSec'   epoch length in seconds (default 2).
%   'Overlap'    fraction of overlap between successive epochs, in [0 1)
%                (default 0.5, i.e. 50% overlap).
%   'SDcrit'     threshold on the (median-scaled) SD of the difference
%                between real and imputed data (default 10).
%   'Rcrit'      threshold on the correlation between real and imputed data;
%                below this a channel/epoch is bad (default 0).
%   'MinRemain'  minimum fraction of recording time that must remain before
%                the next (stricter) channel threshold is tried (default .5).
%   'doDetrend'  detrend each epoch before comparing (default true).
%   'UseParallel' use parfor inside the interpolation step (default false).
%   'Plot'       open a visualisation of what was removed (default false).
%   'Plotter'    which viewer to use when Plot is true (default
%                'vis_artifacts'):
%                  'vis_artifacts' - clean_rawdata viewer; overlays original
%                        (red) and cleaned (dark) data. Removed channels and
%                        removed time segments appear as red-only. Reads the
%                        clean_channel_mask / clean_sample_mask stored in
%                        EEGOut.etc (always set by this function).
%                  'eegplot' - eegplot with a winrej matrix: bad-channel
%                        traces redrawn in red, bad sections a green band.
%                  'none'    - do not plot.
%   'PlotChunkSec' (eegplot only) chunk length (s) used to tile the
%                bad-channel marks so their red traces are drawn on every
%                scroll page; must be <= the eegplot display window
%                (default 1).
%
% OUTPUT
%   EEGOut      cleaned EEG struct (bad channels and bad epochs removed).
%               EEGOut.etc.clean_channel_mask and .clean_sample_mask are
%               set (clean_rawdata convention), so you can visualise the
%               result at any time with:  vis_artifacts(EEGOut, EEGIn)
%   Info        struct with the details of the decision process:
%     .Mask         nchan x nepoch logical, bad channel/epoch combinations
%     .BadChans     logical vector, channels deleted
%     .BadEpochs    logical vector, epochs deleted (bad channels ignored)
%     .R, .SD       the per channel/epoch metrics
%     .EpochStart/.EpochEnd  epoch boundaries in samples
%     .UsedPct      channel-rejection percentage that was finally used
%     .RemainFrac   fraction of recording time retained
%     .RejRegions   merged [start end] sample ranges that were removed
%     .WinRej       the eegplot 'winrej' matrix that was (or would be) drawn
%     .CleanChannelMask  1 x nbchan_orig logical, true = channel kept
%     .CleanSampleMask   1 x pnts_orig  logical, true = sample kept

% ---- parse name-value parameters ---------------------------------------
p = inputParser;
p.FunctionName = 'InterpolationCleaning_claude';
addRequired(p,  'EEGIn',       @isstruct);
addParameter(p, 'ChanRejPct',  [.30 .20], @(x)isnumeric(x)&&isvector(x)&&all(x>0&x<1));
addParameter(p, 'EpochSec',    2,         @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p, 'Overlap',     0.5,       @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<1);
addParameter(p, 'SDcrit',      10,        @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'Rcrit',       0,         @(x)isnumeric(x)&&isscalar(x));
addParameter(p, 'MinRemain',   0.5,       @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p, 'doDetrend',   true,      @(x)islogical(x)||ismember(x,[0 1]));
addParameter(p, 'UseParallel', false,     @(x)islogical(x)||ismember(x,[0 1]));
addParameter(p, 'Plot',        false,     @(x)islogical(x)||ismember(x,[0 1]));
addParameter(p, 'Plotter',    'vis_artifacts', @(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p, 'PlotChunkSec', 1,        @(x)isnumeric(x)&&isscalar(x)&&x>0);
parse(p, EEGIn, varargin{:});

ChanRejPct  = p.Results.ChanRejPct;
EpochSec    = p.Results.EpochSec;
Overlap     = p.Results.Overlap;
SDcrit      = p.Results.SDcrit;
Rcrit       = p.Results.Rcrit;
MinRemain   = p.Results.MinRemain;
doDetrend   = logical(p.Results.doDetrend);
UseParallel = logical(p.Results.UseParallel);
doPlot      = logical(p.Results.Plot);
Plotter     = lower(char(p.Results.Plotter));
PlotChunkSec = p.Results.PlotChunkSec;

if EEGIn.trials ~= 1
    error('InterpolationCleaning_claude requires continuous data (trials==1).');
end

% ---- 1) interpolated surrogate for every channel -----------------------
RealData = EEGIn.data;
tmp = InterpolationReplace(EEGIn, UseParallel);
ImputedData = tmp.data;

% ---- 2) build (overlapping) epoch boundaries in samples ----------------
len  = round(EpochSec * EEGIn.srate);        % epoch length in samples
step = max(1, round(len * (1 - Overlap)));   % advance between epochs
if len > EEGIn.pnts
    error('Epoch length (%d samples) exceeds the data length (%d samples).', len, EEGIn.pnts);
end
EpochStart = (1:step:(EEGIn.pnts - len + 1))';
EpochEnd   = EpochStart + len - 1;
nepochs    = numel(EpochStart);

% ---- 3) quantify the real-vs-imputed difference per channel/epoch ------
nchan = EEGIn.nbchan;
R  = zeros(nchan, nepochs);
SD = zeros(nchan, nepochs);
for e = 1:nepochs
    idx  = EpochStart(e):EpochEnd(e);
    real = RealData(:, idx)';                 % samples x channels
    imp  = ImputedData(:, idx)';
    if doDetrend
        real = detrend(real);
        imp  = detrend(imp);
    end
    R(:, e)  = diag(corr(imp, real));         % similarity per channel
    SD(:, e) = std(imp - real)';              % difference in activity
end
SD = SD ./ median(SD(:));                     % scale to median SD

% ---- 4) channel/epoch artifact mask ------------------------------------
Mask = (R < Rcrit) | (SD > SDcrit);           % nchan x nepoch, true = bad

% ---- 5) pick channel threshold that keeps enough data ------------------
BadFracPerChan = mean(Mask, 2);               % fraction of bad epochs / chan
UsedPct    = ChanRejPct(end);
BadChans   = false(nchan, 1);
BadEpochs  = false(1, nepochs);
RemainFrac = 0;
RejRegions = zeros(0, 2);

for k = 1:numel(ChanRejPct)
    pct = ChanRejPct(k);

    % channels with too many bad epochs are deleted
    BadChans = BadFracPerChan > pct;

    % epochs are bad if any *surviving* channel is bad in them
    if all(BadChans)
        BadEpochs = false(1, nepochs);        % nothing left to judge epochs
    else
        BadEpochs = any(Mask(~BadChans, :), 1);
    end

    % how much recording time survives (epochs overlap -> merge in samples)
    RejRegions = mergeRegions([EpochStart(BadEpochs) EpochEnd(BadEpochs)], EEGIn.pnts);
    rejSamp    = sum(RejRegions(:,2) - RejRegions(:,1) + 1);
    RemainFrac = (EEGIn.pnts - rejSamp) / EEGIn.pnts;

    UsedPct = pct;
    if RemainFrac >= MinRemain
        break;                                % enough data kept, stop
    end
end
if RemainFrac < MinRemain
    warning(['Even at the strictest channel threshold (%.0f%%) only %.1f%% ' ...
             'of the data remains.'], UsedPct*100, RemainFrac*100);
end

% ---- 6) apply the decisions --------------------------------------------
EEGOut = EEGIn;

% remove bad channels
if any(BadChans)
    EEGOut = pop_select(EEGOut, 'nochannel', find(BadChans));
end

% remove bad epochs as continuous time ranges
if ~isempty(RejRegions)
    EEGOut = eeg_eegrej(EEGOut, RejRegions);
end
EEGOut = eeg_checkset(EEGOut);

% ---- 7) store clean_rawdata-style masks (used by vis_artifacts) --------
% Masks are relative to the ORIGINAL data dimensions: true = kept.
% vis_artifacts(EEGOut, EEGIn) reads these to reconstruct a full-size view
% in which removed channels and removed time segments appear as red-only
% (original trace with no cleaned trace over it).
CleanChanMask = ~BadChans(:)';                 % 1 x original nbchan
CleanSampMask = true(1, EEGIn.pnts);           % 1 x original pnts
for r = 1:size(RejRegions, 1)
    CleanSampMask(RejRegions(r,1):RejRegions(r,2)) = false;
end
if sum(CleanChanMask) ~= EEGOut.nbchan
    warning('Channel mask (%d kept) does not match cleaned nbchan (%d).', ...
            sum(CleanChanMask), EEGOut.nbchan);
end
if sum(CleanSampMask) ~= EEGOut.pnts
    warning(['Sample mask (%d kept) does not match cleaned pnts (%d); ' ...
             'vis_artifacts may error.'], sum(CleanSampMask), EEGOut.pnts);
end
EEGOut.etc.clean_channel_mask = CleanChanMask;
EEGOut.etc.clean_sample_mask  = CleanSampMask;

% ---- 8) optional visualisation -----------------------------------------
chunkLen = max(1, round(PlotChunkSec * EEGIn.srate));
WinRej   = buildWinRej(BadChans, RejRegions, nchan, EEGIn.pnts, chunkLen);
if doPlot
    switch Plotter
        case 'vis_artifacts'
            if exist('vis_artifacts', 'file')
                % red = original data, dark = cleaned; removed channels and
                % removed segments show up as red-only.
                vis_artifacts(EEGOut, EEGIn);
            else
                warning(['vis_artifacts (clean_rawdata plugin) not found on ' ...
                         'path; falling back to eegplot.']);
                eegplot(EEGIn.data, 'srate', EEGIn.srate, 'winrej', WinRej, ...
                        'title', 'InterpolationCleaning\_claude');
            end
        case 'eegplot'
            eegplot(EEGIn.data, 'srate', EEGIn.srate, 'winrej', WinRej, ...
                    'title', 'InterpolationCleaning\_claude: red traces = bad channels, green = bad sections');
        case 'none'
            % no plot requested
        otherwise
            error('Unknown Plotter option ''%s'' (use vis_artifacts, eegplot or none).', Plotter);
    end
end

% ---- collect diagnostics ----------------------------------------------
Info = struct();
Info.Mask       = Mask;
Info.BadChans   = BadChans;
Info.BadEpochs  = BadEpochs;
Info.R          = R;
Info.SD         = SD;
Info.EpochStart = EpochStart;
Info.EpochEnd   = EpochEnd;
Info.UsedPct    = UsedPct;
Info.RemainFrac = RemainFrac;
Info.RejRegions = RejRegions;
Info.WinRej     = WinRej;
Info.BadFracPerChan = BadFracPerChan;   % per-channel fraction of bad epochs
Info.MaskFrac       = mean(Mask(:));    % overall fraction of bad chan/epoch combos
Info.CleanChannelMask = CleanChanMask;  % 1 x nbchan_orig, true = channel kept
Info.CleanSampleMask  = CleanSampMask;  % 1 x pnts_orig,  true = sample kept

end % main function


% =======================================================================
function merged = mergeRegions(rej, pnts)
% Sort and merge overlapping/adjacent [start end] sample ranges, clipped to
% the valid data range [1 pnts]. Returns a 0x2 matrix when rej is empty.
if isempty(rej)
    merged = zeros(0, 2);
    return;
end
rej(:,1) = max(rej(:,1), 1);
rej(:,2) = min(rej(:,2), pnts);
rej = sortrows(rej, 1);
merged = rej(1,:);
for i = 2:size(rej,1)
    if rej(i,1) <= merged(end,2) + 1        % overlapping or adjacent
        merged(end,2) = max(merged(end,2), rej(i,2));
    else
        merged = [merged; rej(i,:)]; %#ok<AGROW>
    end
end
end


% =======================================================================
function W = buildWinRej(BadChans, RejRegions, nchan, pnts, chunkLen)
% Build the eegplot 'winrej' matrix: [start end R G B e1 e2 ... enchan].
%
% eegplot renders winrej in two independent ways on continuous data:
%   (1) a FULL-HEIGHT background patch per row, coloured by columns 3:5 and
%       ignoring the channel flags;
%   (2) the flagged channels' traces redrawn in RED, but only in the scroll
%       pages where the row's start OR end is visible.
% To show a bad *channel* (not a full-height band) we therefore give the
% channel rows an invisible (white) background and tile them in short
% chunks so the red trace is redrawn across the whole recording. Bad time
% sections keep a full-height green band with no channel flags.
W = zeros(0, 5 + nchan);

% --- bad channels: invisible band + red trace, tiled in short chunks ----
badMask = double(BadChans(:)');               % 1 x nchan channel flags
if any(badMask)
    starts = (1:chunkLen:pnts)';
    ends   = min(starts + chunkLen - 1, pnts);
    for s = 1:numel(starts)
        % white (invisible) band so only the red channel traces show
        W(end+1, :) = [starts(s), ends(s), 1.0, 1.0, 1.0, badMask]; %#ok<AGROW>
    end
end

% --- bad sections: full-height green band, no per-channel flags ---------
for r = 1:size(RejRegions, 1)
    W(end+1, :) = [RejRegions(r,1), RejRegions(r,2), 0.0, 1.0, 0.0, zeros(1, nchan)]; %#ok<AGROW>
end
end
