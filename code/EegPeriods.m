classdef EegPeriods
% EegPeriods  Shared machinery for the four 'periods' cleaning steps:
% Flat periods, Excessive periods, EMG periods and Interpolation clean.
%
% Every step runs the same pipeline; only step 2 differs per method.
%
%   1. epochs     cut the continuous data into (overlapping) epochs
%   2. mask*      the method's criterion -> channel x epoch logical, true = bad
%   3. decide     channels bad in more than 'maxbadtime' of the epochs are
%                 removed; an epoch is removed when any REMAINING channel is
%                 bad in it; bad epochs become merged sample ranges, which
%                 takes care of the overlap between epochs; while the removed
%                 time exceeds 'maxbadtime', more channels are removed
%   4. apply      pop_select the channels, eeg_eegrej the ranges
%
% Settings every step has (per button, in EegWorkflow_parameters.xlsx):
%   epochlen   epoch length (s)
%   overlap    overlap of successive epochs (proportion, 0 to <1)
%   maxbadtime a channel bad in more than this proportion of the epochs is
%              removed; and while more than this proportion of the recording
%              would be removed, the noisiest channel that helps is removed too
%   mergegap   removed periods closer together than this (s) are joined, and
%              any stretch of data shorter than this between removed periods
%              and/or boundary events is removed too
%
% Only EEG channels with a location, excluding *eog*, are tested
% (eegChannels). The others are never removed.

    methods (Static)

        % ------------------------------------------------------------------
        % Channels that are tested: located, not *eog*. All channels when
        % none qualify.
        function chans = eegChannels(EEG)
            eog     = FindSetNdx({EEG.chanlocs.labels}, '*eog*', 'match', 'pattern');
            located = find(~cellfun(@isempty, {EEG.chanlocs.X}));
            chans   = setdiff(located, eog);
            if isempty(chans)
                chans = 1:EEG.nbchan;
            end
            chans = chans(:)';
        end

        % ------------------------------------------------------------------
        % Epoch boundaries in samples. The length is rounded to an even
        % number of samples (pfft needs that) and at least 2.
        function ep = epochs(EEG, lenSec, overlap)
            len = max(2, 2*round(lenSec*EEG.srate/2));
            if len > EEG.pnts
                error('EegPeriods:epochTooLong', ...
                      'Epoch length (%d samples) exceeds the data (%d samples).', len, EEG.pnts);
            end
            step = max(1, round(len*(1 - overlap)));
            ep.start = (1:step:(EEG.pnts - len + 1))';
            ep.stop  = ep.start + len - 1;
            ep.len   = len;
            ep.n     = numel(ep.start);
        end

        % ------------------------------------------------------------------
        % Surrogate for each tested channel, spherically interpolated from
        % all other located channels. Same result as pop_interp per channel,
        % in one matrix operation (eeg_fastinterpolate).
        function imp = interpolated(EEG, chans)
            imp = eeg_fastinterpolate(EEG, chans);
            imp = reshape(imp, numel(chans), []);
        end

        % ==================================================================
        % Criteria. Each returns numel(chans) x ep.n logical, true = bad.
        % ==================================================================

        % Flat: the (linearly detrended) SD of the epoch is below sdThr uV.
        function M = maskFlat(EEG, chans, ep, sdThr)
            M = false(numel(chans), ep.n);
            for e = 1:ep.n
                x = detrend(double(EEG.data(chans, ep.start(e):ep.stop(e)))');
                M(:, e) = std(x)' < sdThr;
            end
        end

        % Excessive: the epoch SD is unusually HIGH for that channel (z-scored
        % across epochs, upper tail only, FDR corrected at alpha), or the
        % mean-centred amplitude exceeds maxAmp uV anywhere in the epoch.
        % Unusually low SD is left to Flat periods.
        function M = maskExcessive(EEG, chans, ep, alpha, maxAmp)
            sds    = nan(numel(chans), ep.n);
            maxval = sds;
            for e = 1:ep.n
                x = double(EEG.data(chans, ep.start(e):ep.stop(e)));
                sds(:, e)    = std(x, [], 2);
                maxval(:, e) = max(abs(x - mean(x, 2)), [], 2);
            end
            p = normcdf(-zscore(sds, 0, 2));   % P(Z >= z)
            [~, ~, p_fdr] = fdr_bh(p(:));
            M = reshape(p_fdr, size(p)) < alpha | maxval > maxAmp;
        end

        % EMG. Muscle activity is local and broadband: interpolation from the
        % neighbours does not reproduce it, it flattens the spectrum in the
        % band [lo hi) Hz, and it raises the power in that band. Three tests,
        % each switched on separately:
        %
        %   dB test     the interpolated signal has less power in the band than
        %               the recorded one, by more than -thrDb dB (thrDb < 0)
        %   slope test  the spectral slope in the band (least-squares line
        %               through 10*log10(power) against Hz, in dB/Hz) is above
        %               slopeAbs, AND the recorded slope minus the interpolated
        %               slope is above slopeDelta (below it when slopeDelta < 0)
        %   z tests     the band power in dB, z-scored per channel across all
        %               epochs (absolute power, no interpolation), turned into a
        %               one-sided p = P(Z >= z). Two ways to judge it, each
        %               switched on separately:
        %                 p-value   p < zP (fixed, uncorrected)
        %                 FDR       significant after Benjamini-Hochberg over all
        %                           channel x epoch cells at q = zQ, which adapts
        %                           to the number of channels and epochs tested
        %
        % An epoch is bad in a channel when any switched-on test flags it.
        % Slopes are noisy in short epochs; 2 s epochs give steadier estimates.
        function M = maskEMG(EEG, chans, ep, lo, hi, thrDb, useDb, useSlope, slopeAbs, slopeDelta, useZp, zP, useZfdr, zQ)
            if nargin < 7,  useDb = true; end
            if nargin < 8,  useSlope = false; end
            if nargin < 11, useZp = false; end
            if nargin < 13, useZfdr = false; end
            useZ = useZp || useZfdr;
            needImp = useDb || useSlope;
            if needImp
                imp = EegPeriods.interpolated(EEG, chans);
            end
            win = ones(1, ep.len);
            M = false(numel(chans), ep.n);
            bandDb = nan(numel(chans), ep.n);              % for the z test
            for e = 1:ep.n
                idx = ep.start(e):ep.stop(e);
                [Po, fs] = pfft(double(EEG.data(chans, idx))', EEG.srate, win, 0);
                band = fs >= lo & fs < hi;
                if nnz(band) < 2
                    error('EegPeriods:emptyBand', ...
                          'Too few frequencies in %g-%g Hz at %g Hz resolution; lengthen the epochs.', ...
                          lo, hi, fs(2)-fs(1));
                end
                if needImp
                    Pi = pfft(double(imp(:, idx))', EEG.srate, win, 0);
                end
                bad = false(numel(chans), 1);
                if useDb
                    delta = 10*log10(mean(Pi(band, :), 1)) - 10*log10(mean(Po(band, :), 1));
                    bad = bad | delta' < thrDb;
                end
                if useSlope
                    so = EegPeriods.spectralSlope(Po(band, :), fs(band));
                    si = EegPeriods.spectralSlope(Pi(band, :), fs(band));
                    d  = so - si;
                    if slopeDelta >= 0
                        deltaHit = d > slopeDelta;
                    else
                        deltaHit = d < slopeDelta;
                    end
                    bad = bad | (so(:) > slopeAbs & deltaHit(:));
                end
                if useZ
                    bandDb(:, e) = 10*log10(mean(Po(band, :), 1))';
                end
                M(:, e) = bad;
            end
            if useZ
                bandDb(~isfinite(bandDb)) = NaN;           % flat channels: no z
                z = (bandDb - mean(bandDb, 2, 'omitnan')) ./ std(bandDb, 0, 2, 'omitnan');
                if useZp
                    p = 0.5 * erfc(z / sqrt(2));               % NaN z -> NaN p -> not flagged
                    M = M | p < zP;
                end
                if useZfdr
                    M = M | EegPeriods.fdrHigh(z, zQ);
                end
            end
        end

        % ------------------------------------------------------------------
        % Cells whose z is significantly HIGH after Benjamini-Hochberg FDR at
        % level q over all finite cells. One-sided p = P(Z >= z), from erfc so
        % no toolbox is needed.
        function sig = fdrHigh(z, q)
            sig = false(size(z));
            ok  = isfinite(z);
            p   = 0.5 * erfc(z(ok) / sqrt(2));
            m   = numel(p);
            if m == 0
                return
            end
            ps  = sort(p(:));
            k   = find(ps <= (1:m)' / m * q, 1, 'last');
            if isempty(k)
                return
            end
            hit = false(size(p));
            hit(p <= ps(k)) = true;
            sig(ok) = hit;
        end

        % Least-squares slope (dB/Hz) of 10*log10(P) against f, per column of P
        % (freq x chan). Channels with zero power give NaN, which never passes
        % a test.
        function b = spectralSlope(P, f)
            D  = 10*log10(P);
            D(~isfinite(D)) = NaN;
            fc = f(:) - mean(f);
            b  = (fc' * (D - mean(D, 1))) / (fc' * fc);
        end

        % Interpolation: the recorded and interpolated signal correlate
        % below rcrit, or the SD of their difference (scaled to the median
        % over all channels and epochs) exceeds sdcrit.
        % method 'Leave-one-out' predicts each channel from all other
        % channels (eeg_fastinterpolate); 'RANSAC' uses the median over
        % random channel subsets (eeg_ransac, PREP). The tests are the same:
        % per epoch, detrended recorded vs predicted signal, Pearson r below
        % rcrit or SD of the difference above sdcrit times its median.
        function M = maskInterpolation(EEG, chans, ep, sdcrit, rcrit, method, draws, fraction)
            if nargin >= 6 && strcmpi(method, 'RANSAC')
                [R, info] = eeg_ransac(EEG, 'Channels', chans, ...
                    'Windows', [ep.start(:), ep.stop(:)], 'Draws', draws, ...
                    'Fraction', fraction, 'Style', 'detrend');
                SD = info.sd;
            else
                imp = EegPeriods.interpolated(EEG, chans);
                R  = zeros(numel(chans), ep.n);
                SD = R;
                for e = 1:ep.n
                    idx  = ep.start(e):ep.stop(e);
                    real = detrend(double(EEG.data(chans, idx))');
                    surr = detrend(double(imp(:, idx))');
                    R(:, e)  = diag(corr(surr, real));
                    SD(:, e) = std(surr - real)';
                end
            end
            SD = SD ./ median(SD(:));
            M = R < rcrit | SD > sdcrit;
        end

        % ==================================================================
        % Shared decision and removal
        % ==================================================================

        % Channels and periods to remove, given the channel x epoch mask M.
        %
        %   1. a channel bad in more than maxbadtime of the epochs is removed
        %   2. an epoch is bad when any remaining channel is bad in it; the bad
        %      epochs become merged sample ranges (overlap, mergegap)
        %   3. while those ranges still cover more than maxbadtime of the
        %      recording, one more channel is removed: the noisiest remaining
        %      channel (highest proportion of bad epochs) whose removal lowers
        %      the removed time. Stops when the removed time is within
        %      maxbadtime, when no channel helps, or when one channel is left.
        %
        % D.badChans     logical over chans, all channels removed
        % D.extraChans   logical over chans, the ones added by step 3
        % D.fracBad      proportion of bad epochs per channel
        % D.badEpochs    logical over epochs, removed
        % D.regions      [start stop] sample ranges to remove, merged
        % D.remain       proportion of the recording that survives
        % D.maxbadtime   the criterion used
        function D = decide(M, ep, pnts, srate, maxbadtime, mergegap)
            gap     = round(mergegap*srate);
            fracBad = mean(M, 2);
            nChan   = size(M, 1);

            badChans = fracBad > maxbadtime;
            extra    = false(nChan, 1);
            [regions, remain, badEpochs] = EegPeriods.periodsFor(M, badChans, ep, pnts, gap);

            while (1 - remain) > maxbadtime && sum(~badChans) > 1
                cand = find(~badChans);
                [~, order] = sort(fracBad(cand), 'descend');
                helped = false;
                for c = cand(order)'
                    trial = badChans;
                    trial(c) = true;
                    [r2, rem2, be2] = EegPeriods.periodsFor(M, trial, ep, pnts, gap);
                    if rem2 > remain
                        badChans = trial;
                        extra(c) = true;
                        regions = r2; remain = rem2; badEpochs = be2;
                        helped = true;
                        break
                    end
                end
                if ~helped
                    break
                end
            end

            D.badChans   = badChans;
            D.extraChans = extra;
            D.fracBad    = fracBad;
            D.badEpochs  = badEpochs;
            D.regions    = regions;
            D.remain     = remain;
            D.maxbadtime = maxbadtime;
        end

        % ------------------------------------------------------------------
        % Periods removed when the channels in removed are left out.
        function [regions, remain, badEpochs] = periodsFor(M, removed, ep, pnts, gap)
            if all(removed)
                badEpochs = false(1, ep.n);
            else
                badEpochs = any(M(~removed, :), 1);
            end
            regions = EegPeriods.mergeRegions([ep.start(badEpochs) ep.stop(badEpochs)], pnts, gap);
            remain  = 1 - sum(diff(regions, 1, 2) + 1)/pnts;
        end

        % ------------------------------------------------------------------
        % Remove short islands of data: a stretch shorter than gap samples with
        % a removed period or an existing boundary event on BOTH sides is
        % removed as well. mergeRegions only joins this step's own periods;
        % this also catches islands next to boundaries left by earlier steps.
        % The start and end of the recording count as a boundary.
        function D = absorbShortIslands(D, EEG, gap)
            D.islands = 0;
            pnts = EEG.pnts;
            if gap < 1 || isempty(D.regions)
                return
            end
            removed = false(1, pnts);
            for r = 1:size(D.regions, 1)
                removed(D.regions(r,1):D.regions(r,2)) = true;
            end
            [removed, D.islands] = EegPeriods.markShortStretches(removed, EEG, gap);
            if D.islands > 0
                d = diff([false removed false]);
                D.regions = [find(d == 1)' find(d == -1)' - 1];
                D.remain  = 1 - nnz(removed)/pnts;
            end
        end

        % ------------------------------------------------------------------
        % Remove every stretch of data shorter than gap samples that has a
        % boundary event on both sides (e.g. left between bursts that
        % clean_rawdata cut out). The start and end of the recording count as
        % a boundary. n = number of stretches, sec = seconds removed.
        function [EEG, n, sec] = removeShortStretches(EEG, gap)
            n = 0; sec = 0;
            if gap < 1 || EEG.trials > 1
                return
            end
            [removed, n] = EegPeriods.markShortStretches(false(1, EEG.pnts), EEG, gap);
            if n == 0
                return
            end
            d = diff([false removed false]);
            regions = [find(d == 1)' find(d == -1)' - 1];
            sec = nnz(removed) / EEG.srate;
            EEG.icaact = [];
            EEG = eeg_eegrej(EEG, regions);
            EEG.history = [EEG.history newline ...
                sprintf('EEG = eeg_eegrej(EEG, %s); %% stretches shorter than %g s between boundaries', ...
                mat2str(regions), gap/EEG.srate)];
        end

        % ------------------------------------------------------------------
        % Mark (in removed) every run of kept samples shorter than gap that is
        % bounded on both sides by a removed sample, a boundary event, or the
        % start/end of the recording. Nothing is marked if that would remove
        % all data (e.g. a recording shorter than gap without boundaries).
        function [removed, n] = markShortStretches(removed, EEG, gap)
            pnts = numel(removed);
            % cut(k): the data are discontinuous between sample k-1 and k
            cut = false(1, pnts);
            if isfield(EEG, 'event') && ~isempty(EEG.event)
                isB = strcmpi({EEG.event.type}, 'boundary');
                k   = floor([EEG.event(isB).latency]) + 1;
                cut(k(k > 1 & k <= pnts)) = true;
            end
            startsRun = ~removed & [true, removed(1:end-1) | cut(2:end)];
            endsRun   = ~removed & [removed(2:end) | cut(2:end), true];
            s0 = find(startsRun);
            s1 = find(endsRun);
            short = (s1 - s0 + 1) < gap;
            if all(short) || nnz(~removed) == 0
                n = 0;                                  % would leave no data
                return
            end
            for i = find(short)
                removed(s0(i):s1(i)) = true;
            end
            n = nnz(short);
        end

        % ------------------------------------------------------------------
        % Sort and join [start stop] ranges that overlap, touch, or lie
        % within gap samples of each other; clipped to 1..pnts.
        function merged = mergeRegions(rej, pnts, gap)
            merged = zeros(0, 2);
            if isempty(rej)
                return
            end
            rej(:,1) = max(rej(:,1), 1);
            rej(:,2) = min(rej(:,2), pnts);
            rej = sortrows(rej, 1);
            merged = rej(1,:);
            for i = 2:size(rej,1)
                if rej(i,1) <= merged(end,2) + 1 + gap
                    merged(end,2) = max(merged(end,2), rej(i,2));
                else
                    merged(end+1,:) = rej(i,:); %#ok<AGROW>
                end
            end
        end

        % ------------------------------------------------------------------
        % Remove the decided channels and periods.
        function EEG = apply(EEG, chans, D)
            if ~isempty(D.regions)
                % a stored icaact no longer matches once samples go, and
                % eeg_checkset (inside eeg_eegrej) then errors; it is
                % recomputed from the weights when needed
                EEG.icaact = [];
                EEG = eeg_eegrej(EEG, D.regions);
                EEG.history = [EEG.history newline ...
                    sprintf('EEG = eeg_eegrej(EEG, %s);', mat2str(D.regions))];
            end
            if any(D.badChans)
                [EEG, cmd] = pop_select(EEG, 'nochannel', chans(D.badChans));
                EEG.history = [EEG.history newline cmd];
            end
            EEG = eeg_checkset(EEG);
        end

    end
end
