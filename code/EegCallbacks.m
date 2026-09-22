classdef EegCallbacks
% EegCallbacks  Static class holding callback/helper logic for
% guiEegAutoflow_App. In App Designer each callback stub becomes:
%     EegCallbacks.pushbuttonFilter_Callback(app, event);
% hObject is always app.eeg_workflow (the main uifigure handle).

    % Colours of the workflow buttons in the left-hand column.
    %   idle / done : light red   - not run yet, or finished
    %   ready       : light green - a file is open, this step can be run
    %   busy        : dark green  - running right now
    % Open starts out ready and every other left button starts out idle;
    % opening a file flips that round. Buttons outside the left column keep
    % the colour they were given in App Designer.
    properties (Constant)
        ColIdle  = [1 .6 .6];
        ColReady = [.6 1 .6];
        ColBusy  = [.3 .6 .3];
    end

    methods (Static)

        % Button pushed function: pushbuttonFilter
        function pushbuttonFilter_Callback(app, event)
            hObject = app.eeg_workflow;

            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            app.pushbuttonFilter.BackgroundColor = EegCallbacks.ColBusy;
            pause(0.005)

            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.AddToListbox(app, app.listboxStdout,'*** error *** No EEG data avilable');
                EegCallbacks.abortStep(app, event, 'No data available');
                app.pushbuttonFilter.BackgroundColor = EegCallbacks.ColIdle;
                return
            end

            tmp = data.EEG;
            hp  = P('filter','low');      % 0 = no highpass
            lp  = P('filter','high');
            fsr = tmp.srate;
            if lp >= fsr/2
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                    ' *** warning *** %.1f Hz is not below Nyquist (%.1f Hz): no lowpass', lp, fsr/2));
                lp = 0;
            end
            if hp <= 0 && lp <= 0 || (hp > 0 && lp > 0 && hp >= lp)
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                    ' *** error *** no valid band from %.1f to %.1f Hz; nothing filtered', hp, lp));
                EegCallbacks.abortStep(app, event, 'Filter: lower bound must be below upper bound.');
                app.pushbuttonFilter.BackgroundColor = EegCallbacks.ColIdle;
                return
            end
            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf('Filter: highpass %s, lowpass %s', ...
                ifthen(hp > 0, sprintf('%.2f Hz', hp), 'none'), ifthen(lp > 0, sprintf('%.2f Hz', lp), 'none')));

            if strcmpi(P('filter','type'), 'Butterworth')
                % Zero-phase IIR (filtfilt). HalfPowerFrequency is the -3 dB
                % point of one pass, so -6 dB after the forward-backward pass.
                % Steepness comes from the order; 'transition' does not apply.
                bOrd = 4;
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                    '- zero-phase Butterworth, order %d (%d effective), edges at -6 dB', bOrd, 2*bOrd));
                tmp = filter_butter(tmp, fsr, hp, lp, bOrd, true, false, false);
                tmp.history = [tmp.history newline sprintf( ...
                    '%% filter_butter %.2f-%.2f Hz, order %d, zero phase', hp, lp, bOrd)];
            else
                % Linear-phase FIR (Hamming windowed sinc, pop_eegfiltnew), one
                % pass with group-delay correction: cutoffs are the -6 dB points,
                % centred in the transition band. The order follows from the
                % transition bandwidth (Widmann et al. 2015). Highpass and
                % lowpass are separate filters: one windowed-sinc kernel has a
                % single length and so a single transition width, and a narrow
                % highpass band would otherwise make the lowpass needlessly long.
                if hp > 0
                    df   = P('filter','transition_hp');
                    dfHp = min(df, hp);             % stopband edge stays >= 0 Hz
                    n = EegCallbacks.firOrderHamming(fsr, dfHp);
                    [tmp, cmd] = pop_eegfiltnew(tmp, hp, [], n, 0);
                    tmp.history = [tmp.history newline cmd];
                    EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                        '- FIR highpass %.2f Hz, transition %.2f Hz (%.2f-%.2f), order %d%s', ...
                        hp, dfHp, hp-dfHp/2, hp+dfHp/2, n, ifthen(dfHp < df, ' (transition capped at cutoff)', '')));
                end
                if lp > 0
                    df   = P('filter','transition_lp');
                    dfLp = min(df, fsr/2 - lp);     % stopband edge stays <= Nyquist
                    n = EegCallbacks.firOrderHamming(fsr, dfLp);
                    [tmp, cmd] = pop_eegfiltnew(tmp, [], lp, n, 0);
                    tmp.history = [tmp.history newline cmd];
                    EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                        '- FIR lowpass %.2f Hz, transition %.2f Hz (%.2f-%.2f), order %d%s', ...
                        lp, dfLp, lp-dfLp/2, lp+dfLp/2, n, ifthen(dfLp < df, ' (transition capped at Nyquist)', '')));
                end
            end

            % push existing data onto stack. Update <data.EEG> to tmp.
            % remember the passband; the spectrum tracker limits its x-axis to it
            tmp.etc.filterBand = [hp, ifthen(lp > 0, lp, Inf)];
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Filter');
            data.EEG = tmp;

            guidata(hObject, data);

            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'Filter');
            guidata(hObject, data);

            app.pushbuttonFilter.BackgroundColor = EegCallbacks.ColIdle;
        end

        % Code that executes after component creation
        function guiEegAutoflow_OpeningFcn(app, varargin)
            hObject = app.eeg_workflow;
            movegui(hObject, 'onscreen');

            data = guidata(hObject);

            EegCallbacks.AddToListbox(app, app.listboxStdout, '  *** Warning *** EEGLAB with specific plugins is required')
            EegCallbacks.AddToListbox(app, app.listboxStdout, '   - EEGLAB V2020 has been tested, requires signal processing toolbox')
            EegCallbacks.AddToListbox(app, app.listboxStdout, '   - AAR')
            EegCallbacks.AddToListbox(app, app.listboxStdout, '   - CleanRawdata')
            EegCallbacks.AddToListbox(app, app.listboxStdout, '   - ICLabel')
            EegCallbacks.AddToListbox(app, app.listboxStdout, '   - file import plugins (ANT, Biosemi, EDF)')
            EegCallbacks.AddToListbox(app, app.listboxStdout, '   - and several support functions')

            if isempty(which("eeglab.m"))
                EegCallbacks.AddToListbox(app, app.listboxStdout, '*** warning *** cannot find EEGLAB. Please activate locate.')
                filepath = uigetdir(pwd, "Locate EEGLAB");
                EegCallbacks.bringToFront(app.eeg_workflow);
                addpath(filepath)
                eeglab;
            end

            % SETTINGSDIR: the system's application-data folder, where every
            % settings file lives. DEFAULTDIR: the folder file dialogs (Open,
            % batch file selection) start in; updated whenever files are chosen.
            data.SETTINGSDIR = EegCallbacks.settingsDir();
            [nMoved, oldDir] = EegCallbacks.migrateSettings(data.SETTINGSDIR);
            if nMoved > 0
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                    '- copied %d settings files from %s to %s', nMoved, oldDir, data.SETTINGSDIR));
            end
            data.DEFAULTDIR = EegCallbacks.loadDefaultDir(data.SETTINGSDIR);

            % Read the parameter definitions from EegWorkflow_parameters.xlsx,
            % merge in the values saved last session and show them next to the
            % buttons they belong to. On a first run there is no settings file
            % yet and the initialpardefault column of the spreadsheet is used.
            data.params = EegParams.load(data.SETTINGSDIR);
            EegParams.refresh(app, data.params);

            % nothing loaded yet: Open is ready, every other step is idle
            EegCallbacks.resetButtonColors(app);

            % restore the value controls of the main window (ICA type, ICA
            % number, code text area, ...) saved on the last close
            strlist = table();
            try
                FN = sprintf('%s/%s.ini', data.SETTINGSDIR, get(hObject,'name'));
                % tab only: writetable writes tabs, and letting MATLAB guess
                % would split texts containing commas or spaces
                opts = detectImportOptions(FN, 'TextType', 'string', 'filetype', 'text', 'Delimiter', '\t');
                opts.DataLines = [2 Inf];
                opts.VariableTypes(:) = {'string'};
                strlist = readtable(FN, opts);
                EegCallbacks.SetUIControlData(app, hObject, strlist);
            catch
                warning('Initialization file not found. Will be created on close.')
            end

            if ispc
                data.fontsize = 9;
            else
                data.fontsize = 10;
            end
            % font size of the last session: an extra 'fontsize' line in the
            % same ini file (SetUIControlData skips keys that are no control)
            if ismember('key', strlist.Properties.VariableNames)
                ndx = find(strcmpi(strlist.key, 'fontsize'), 1);
                fs = str2double(strlist.val(ndx));
                if ~isempty(ndx) && ~isnan(fs)
                    data.fontsize = min(max(round(fs), 3), 16);   % range of the up/down buttons
                end
            end
            % hObject, not hObject.Parent: the parent of the figure is the
            % graphics root, and walking that restyles every open figure.
            EegCallbacks.setFontSize(app, hObject, data.fontsize)

            guidata(hObject, data);
        end

        % Close request function: fig_eeg_workflow
        function eeg_workflow_CloseRequestFcn(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            try
                data = guidata(hObject);
                if ~exist(data.SETTINGSDIR,'dir')
                    mkdir(data.SETTINGSDIR)
                end
                if isfield(data,'params')
                    EegParams.save(data.params, data.SETTINGSDIR);
                end
                strlist = EegCallbacks.GetUIControlData(app, hObject);
                if isfield(data, 'fontsize')
                    strlist = [strlist; table({'fontsize'}, {sprintf('%d', data.fontsize)}, ...
                                              'VariableNames', {'key','val'})];
                end
                writetable(strlist,sprintf('%s/%s.ini', data.SETTINGSDIR, get(hObject,'name')), 'delimiter','\t','filetype','text')
                files = dir('runica*');
                for f = 1:length(files)
                    delete(fullfile(files(f).folder, files(f).name));
                end
            catch
                warndlg('Error saving settings.')
            end
            delete(findall(groot, 'Type', 'figure', 'Tag', 'EegTrackViewer'));   % hidden viewer
            delete(hObject);
        end

        % ---- Settings listboxes ------------------------------------------
        % Every listbox next to a workflow button shows that button's settings
        % as read-only '<key>=<value>' lines and opens the editor when clicked.
        % In App Designer both Clicked callbacks (ListBoxOpenClicked and
        % ListBoxFlatlineClicked) only need the single line:
        %     EegCallbacks.ListBoxClicked(app, event);
        function ListBoxClicked(app, event)
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            if ~isfield(data, 'params')
                data.params = EegParams.load(data.SETTINGSDIR);
            end

            % Which listbox was clicked, and so which button it belongs to.
            % event.Source is what a Clicked callback provides; CurrentObject
            % is the fallback for callers that pass no usable event.
            m = EegParams.listboxMap();
            src = [];
            try
                src = event.Source;
            catch
            end
            [lbName, buttonKey, dlgTitle] = EegCallbacks.matchListbox(app, m, src);
            if isempty(lbName)
                try
                    src = hObject.CurrentObject;
                catch
                    src = [];
                end
                [lbName, buttonKey, dlgTitle] = EegCallbacks.matchListbox(app, m, src);
            end
            if isempty(lbName)
                % Never fail silently: a settings listbox that does nothing
                % when clicked is otherwise indistinguishable from a missing
                % callback body in the .mlapp.
                EegCallbacks.AddToListbox(app, app.listboxStdout, ...
                    '*** warning *** clicked listbox not recognised as a settings listbox');
                return
            end

            % clicking a line selects it; the listbox is display-only
            EegParams.refreshOne(app, data.params, lbName, buttonKey);

            % The dialog returns immediately; everything below happens in
            % applyParamEdit when the user presses OK.
            EegParams.editDialog(hObject, data.params, buttonKey, dlgTitle, ...
                @(Tnew) EegCallbacks.applyParamEdit(app, lbName, buttonKey, Tnew));
        end

        % ------------------------------------------------------------------
        % Called by the settings dialog when OK is pressed. guidata is read
        % fresh here rather than captured, so nothing goes stale while the
        % dialog is open.
        function applyParamEdit(app, lbName, buttonKey, Tedited)
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            data.params = EegParams.mergeButton(data.params, Tedited, buttonKey);
            guidata(hObject, data);
            EegParams.refreshOne(app, data.params, lbName, buttonKey);
            EegParams.save(data.params, data.SETTINGSDIR);
        end

        % ------------------------------------------------------------------
        function [lbName, buttonKey, dlgTitle] = matchListbox(app, m, src)
            lbName = ''; buttonKey = ''; dlgTitle = '';
            if isempty(src) || ~isgraphics(src)
                return
            end
            for k = 1:size(m,1)
                if isprop(app, m{k,1}) && isvalid(app.(m{k,1})) && app.(m{k,1}) == src
                    lbName    = m{k,1};
                    buttonKey = m{k,2};
                    dlgTitle  = m{k,3};
                    return
                end
            end
        end

        % ---- Font size ----------------------------------------------------
        function pushbuttonUp_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            if data.fontsize>15, return; end
            data.fontsize = data.fontsize + 1;
            EegCallbacks.setFontSize(app, hObject, data.fontsize)
            guidata(hObject, data);
        end

        function pushbuttonDOWN_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            if data.fontsize<4, return; end
            data.fontsize = data.fontsize - 1;
            EegCallbacks.setFontSize(app, hObject, data.fontsize)
            guidata(hObject, data);
        end

        % ---- Simple data ops ---------------------------------------------
        function pushbuttonRemoveFirstSec_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            EegCallbacks.AddToListbox(app, app.listboxStdout, '- removing first 1 second of data');
            [tmp, cmd] = pop_select(data.EEG, 'rmtime', [0 1]);
            tmp.history = [tmp.history char(uint8(10)) cmd];
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Rm first 1s');
            data.EEG = tmp;
            guidata(hObject, data);
            EegCallbacks.trackSpectrum(app, data.EEG, 'Remove 1st s');
        end

        function pushbuttonLastSecButtonPushed(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            EegCallbacks.AddToListbox(app, app.listboxStdout, '- removing last 1 second of data');
            [tmp, cmd] = pop_select(data.EEG, 'rmtime', [data.EEG.xmax-1 data.EEG.xmax]);
            tmp.history = [tmp.history char(uint8(10)) cmd];
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Rm last 1s');
            data.EEG = tmp;
            guidata(hObject, data);
            EegCallbacks.trackSpectrum(app, data.EEG, 'Remove last s');
        end

        function pushbuttonLast4minButtonPushed(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Extracting last 4 minutes of data');
            if data.EEG.xmax>60*4
                [tmp, cmd] = pop_select(data.EEG, 'time', [data.EEG.xmax-60*4 data.EEG.xmax]);
                tmp.history = [tmp.history char(uint8(10)) cmd];
                data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Keep last 4 min');
                data.EEG = tmp;
                guidata(hObject, data);
            else
                EegCallbacks.AddToListbox(app, app.listboxStdout, ' - data too short. No change.');
            end
            EegCallbacks.trackSpectrum(app, data.EEG, 'Keep last 4 min');
        end

        function pushbuttonSaveMemory_Callback(app, event) %#ok<INUSD>
            data = guidata(app.eeg_workflow);
            global GlobEEG %#ok<GVMIS>
            GlobEEG = data.EEG;
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Saving data into global variable GlobEEG');
        end

        function pushbuttonMemoryBack_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            global GlobEEG %#ok<GVMIS>
            data.EEG = GlobEEG;
            guidata(hObject, data);
        end

        function pushbuttonPlot2D_Callback(app, event) %#ok<INUSD>
            data = guidata(app.eeg_workflow);
            figure;
            topoplot([],data.EEG.chanlocs, 'style', 'blank', ...
                'electrodes', 'labelpoint', ...
                'chaninfo', data.EEG.chaninfo);
        end

        function pbUndo_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            if isempty(data.Stack)
                EegCallbacks.abortStep(app, event, 'No more saved datasets.');
                EegCallbacks.AddToListbox(app, app.listboxStdout, 'No more saved datasets.');
            else
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf('Undo %s',data.StackLabel{end}));
                EegCallbacks.RemoveLastHistory(app);
                data.EEG = data.Stack{end};
                if length(data.Stack)==1
                    data.Stack = {};
                    data.StackLabel = {};
                else
                    data.Stack  = data.Stack(1:end-1);
                    data.StackLabel  = data.StackLabel(1:end-1);
                end
                guidata(hObject, data);
                EegCallbacks.untrackSpectrum(app, data.EEG, numel(data.Stack));
            end
            guidata(hObject, data);
        end

        % Alias for the typo in the .mlapp stub of pbEpoch, which calls
        % EegCallbacks.bEpoch_Callback. Delete once the stub is corrected.
        function bEpoch_Callback(app, event)
            EegCallbacks.pbEpoch_Callback(app, event);
        end

        function pbEpoch_Callback(app, event) %#ok<INUSD>
            data = guidata(app.eeg_workflow);
            h = guiEpoch(data.EEG, gcf);
            uiwait(h);
        end

        % Button pushed function: pushbuttonFlatline
        function pushbuttonFlatline_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            set(app.pushbuttonFlatline, 'BackgroundColor', EegCallbacks.ColBusy);
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                set(app.pushbuttonFlatline, 'BackgroundColor', EegCallbacks.ColIdle);
                return
            end
            tmp = data.EEG;
            crit = P('flatline','sd');
            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf('Removing channels with <%.1f stdev.', crit));
            SD = std(data.EEG.data(:,:)');
            ndx = find(SD < crit);
            if ~isempty(ndx)
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf('- Removing %d channels ', length(ndx)));
                [tmp, cmd] = pop_select(tmp,'nochannel',ndx);
                tmp.history = [tmp.history char(uint8(10)) cmd];
                data.tabLine.flatline = length(ndx);
            else
                EegCallbacks.AddToListbox(app, app.listboxStdout, '- NO channels removed');
                data.tabLine.flatline = 0;
            end
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Flatline');
            data.EEG = tmp;
            guidata(hObject, data);
            set(app.pushbuttonFlatline, 'BackgroundColor', EegCallbacks.ColIdle);
            EegCallbacks.trackSpectrum(app, data.EEG, 'Flatline');
        end

        % Button pushed function: pushbuttonExcessive
        function pushbuttonExcessive_Callback(app, event)
            % Bad-channel removal. Every test has an on/off flag and its own
            % threshold, so the step can run the original SD z-score test, a
            % fixed SD limit, the clean_rawdata channel tests DISCOVER-EEG uses
            % for its bad-channel pass (flatline, correlation, line noise, drift
            % highpass; thresholds default to DISCOVER-EEG's), or any mix.
            % The clean_rawdata tests run first; the SD tests then run on what
            % is left (after the drift highpass, when that is on). RANSAC (PREP's
            % bad-channel test, eeg_ransac) runs last, as in PREP: channels the
            % other tests removed are no longer used to predict the others.
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(p) EegParams.get(data.params, 'excessive signal', p);
            say = @(varargin) EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(varargin{:}));
            set(app.pushbuttonExcessive, 'BackgroundColor', EegCallbacks.ColBusy);
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                set(app.pushbuttonExcessive, 'BackgroundColor', EegCallbacks.ColIdle);
                return
            end
            tmp = data.EEG;
            before = {tmp.chanlocs.labels};
            say('Removing bad channels');

            % ---- clean_rawdata channel tests ---------------------------------
            useFlat  = P('use_flatline');
            useMinr  = P('use_minr');
            useNoise = P('use_noise');
            useDrift = P('use_drift');
            if useFlat || useMinr || useNoise || useDrift
                band = sscanf(strrep(P('drift'), '-', ' '), '%f')';
                chanArg = 'off';
                if useMinr
                    chanArg = P('minr');
                elseif useNoise
                    chanArg = -1;          % the line-noise test needs a numeric criterion; -1 never flags
                end
                say('- clean_rawdata: flatline %s | min r %s | line noise %s | drift highpass %s | max bad time %.2f', ...
                    ifthen(useFlat, sprintf('%g s', P('flatline')), 'off'), ...
                    ifthen(useMinr, sprintf('%.2f', P('minr')), 'off'), ...
                    ifthen(useNoise, sprintf('%g SD', P('noise')), 'off'), ...
                    ifthen(useDrift, sprintf('%g-%g Hz', band), 'off'), P('maxbadtime'));
                [tmp, cmd] = pop_clean_rawdata(tmp, ...
                    'FlatlineCriterion', ifthen(useFlat, P('flatline'), 'off'), ...
                    'ChannelCriterion', chanArg, ...
                    'LineNoiseCriterion', ifthen(useNoise, P('noise'), 'off'), ...
                    'Highpass', ifthen(useDrift, band, 'off'), ...
                    'BurstCriterion', 'off', 'WindowCriterion', 'off', 'BurstRejection', 'off', ...
                    'Distance', 'Euclidian', ...
                    'ChannelCriterionMaxBadTime', P('maxbadtime'));
                tmp.history = [tmp.history newline cmd];
                gone = setdiff(before, {tmp.chanlocs.labels}, 'stable');
                say('  removed %d: %s', numel(gone), strjoin(gone, ' '));
            end

            % ---- SD tests -----------------------------------------------------
            SD = std(double(tmp.data(:,:)), [], 2)';
            bad = false(1, tmp.nbchan);
            if P('use_z')
                Z = (SD - mean(SD)) ./ std(SD);
                hit = Z > P('excessive_z');
                say('- SD z-score > %.1f: %d channels %s', P('excessive_z'), sum(hit), strjoin({tmp.chanlocs(hit).labels}, ' '));
                bad = bad | hit;
            end
            if P('use_sd')
                hit = SD > P('excessive_sd');
                say('- SD > %g uV: %d channels %s', P('excessive_sd'), sum(hit), strjoin({tmp.chanlocs(hit).labels}, ' '));
                bad = bad | hit;
            end
            if all(bad)
                say(' *** error *** every channel failed the SD tests; those are not removed');
                bad(:) = false;
            end
            if any(bad)
                [tmp, cmd] = pop_select(tmp, 'nochannel', find(bad));
                tmp.history = [tmp.history newline cmd];
            end

            % ---- RANSAC (PREP) ------------------------------------------------
            if P('use_ransac')
                hit = EegCallbacks.ransacBadChannels(app, tmp, P);
                if ~isempty(hit)
                    [tmp, cmd] = pop_select(tmp, 'nochannel', hit);
                    tmp.history = [tmp.history newline cmd];
                end
            end

            gone = setdiff(before, {tmp.chanlocs.labels}, 'stable');
            data.tabLine.excessive = numel(gone);
            if isempty(gone) && isequal(tmp.data, data.EEG.data)
                say('- NO channels removed');
            else
                say('- %d channels removed in total', numel(gone));
                data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Bad chans');
                data.EEG = tmp;
            end
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'Excessive');
            guidata(hObject, data);
            set(app.pushbuttonExcessive, 'BackgroundColor', EegCallbacks.ColIdle);
        end

        % Button pushed function: pushbuttonResample
        function pushbuttonResample_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            set(app.pushbuttonResample,'backgroundcolor', EegCallbacks.ColBusy)
            pause(0.005);
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            EegCallbacks.AddToListbox(app, app.listboxStdout,'Resampling data');
            tmp = data.EEG;
            fs_new = str2double(P('resample','srate'));
            fs_old = tmp.srate;
            if isnan(fs_new) || fs_new==fs_old
                EegCallbacks.AddToListbox(app, app.listboxStdout,'*** warning *** no change in sampling rate.');
                set(app.pushbuttonResample,'backgroundcolor', EegCallbacks.ColIdle)
                return
            end
            isint = (fs_new/fs_old==round(fs_new/fs_old)) | (fs_old/fs_new==round(fs_old/fs_new));
            % pop_resample low-pass filters (anti-aliasing) for any ratio. The
            % 'spline' method is the older route for non-integer ratios: it
            % interpolates the UNFILTERED data, so content above the new
            % Nyquist frequency aliases.
            useSpline = ~isint && strcmpi(P('resample','method'), 'spline');
            if ~isint
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                    '- non-integer ratio %g -> %g Hz, using %s', fs_old, fs_new, ifthen(useSpline, 'spline interpolation (no anti-aliasing)', 'pop_resample')));
            end
            if ~useSpline
                [tmp, cmd] = pop_resample(tmp, fs_new);
            else
                warning('Resampling by spline interpolation!')
                dummy = pop_resample(tmp, fs_new);
                t_old = (0:tmp.pnts-1)/fs_old;
                t_new = (0:dummy.pnts-1)/fs_new;
                tmp = dummy;
                for ch=1:tmp.nbchan
                    tmp.data(ch,:) = interp1(t_old, data.EEG.data(ch,:), t_new, 'spline');
                end
                cmd = '% custom spline resampling';
            end
            tmp.history = [tmp.history char(uint8(10)) cmd];
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Resample');
            data.EEG = tmp;
            guidata(hObject, data);
            set(app.pushbuttonResample,'backgroundcolor', EegCallbacks.ColIdle)
            EegCallbacks.trackSpectrum(app, data.EEG, 'Resample');
        end

        % Button pushed function: pushbuttonRereference
        function pushbuttonRereference_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            set(app.pushbuttonRereference, 'BackgroundColor', EegCallbacks.ColBusy);
            pause(.005);
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                app.pushbuttonRereference.BackgroundColor = EegCallbacks.ColIdle;
                return
            end
            tmp = data.EEG;
            if size(tmp.data,3) ~= tmp.trials
                tmp.trials = size(tmp.data,3);
            end
            if tmp.trials==1 && tmp.pnts~=size(tmp.data,2)
                tmp.pnts = size(tmp.data,2);
            end
            refchans = P('rereference','refchans');
            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf('Rereferencing data (%s)', refchans));
            cmd = '';
            switch upper(refchans)
                case 'CPZ', [tmp, cmd] = pop_reref(tmp, find(ismember(upper({tmp.chanlocs.labels}),'CPZ')));
                case 'M1/M2', [tmp, cmd] = pop_reref(tmp, find(ismember(upper({tmp.chanlocs.labels}),{'M1','M2'})));
                case 'AVERAGE'
                    [tmp, cmd] = EegCallbacks.averageReference(app, tmp, P('rereference','interpremoved'));
                case 'REST'
                    tmp = eeg_REST_reref(tmp);
                    cmd = 'eeg_REST_reref(tmp);';
                case 'A1/A2', [tmp, cmd] = pop_reref(tmp, find(ismember(upper({tmp.chanlocs.labels}),{'A1','A2'})));
                case 'CSD'
                    trodes = {};
                    for ch=1:tmp.nbchan
                        if ~isempty(tmp.chanlocs(ch).X)
                            trodes = cat(1, trodes, {tmp.chanlocs(ch).labels});
                        end
                    end
                    if length(trodes)<tmp.nbchan
                        warning('DOWNSIZING DATA. Only channels with  location info can be used for CSD')
                        tmp = pop_select(tmp,'channel',trodes);
                    end
                    Montage_64 = ExtractMontage('/data/damitsea/Matlab/CSDtoolbox/resource/10-5-System_Mastoids_EGI129.csd',trodes);
                    [G,H] = GetGH(Montage_64);
                    [s1,s2,s3] = size(tmp.data);
                    tmp.data = CSD(tmp.data(:,:),G,H);
                    if s3>1
                        tmp.data = reshape(tmp.data, [s1,s2,s3]);
                    end
                    cmd = '%% CSD';
                case 'TP9/10'
                    ndx = find(ismember(upper({tmp.chanlocs.labels}),{'TP9','TP10'}));
                    if ~isempty(ndx)
                        [tmp, cmd] = pop_reref(tmp, ndx);
                    else
                        EegCallbacks.AddToListbox(app, app.listboxStdout, 'Channels TP9/10 not found');
                        cmd = '';
                    end
                case 'P9/10'
                    ndx = find(ismember(upper({tmp.chanlocs.labels}),{'P9','P10'}));
                    if ~isempty(ndx)
                        [tmp, cmd] = pop_reref(tmp, ndx);
                    else
                        EegCallbacks.AddToListbox(app, app.listboxStdout, 'Channels P9/10 not found');
                        cmd = '';
                    end
                case 'E56/E57/E100/E107'
                    ndx = find(ismember(upper({tmp.chanlocs.labels}),{'E56','E57','E100','E107'}));
                    if length(ndx)==4
                        [tmp, cmd] = pop_reref(tmp, ndx);
                    else
                        ndx = find(ismember(upper({tmp.chanlocs.labels}),{'E56','M1','M2','E107'}));
                        if length(ndx)==4
                            [tmp, cmd] = pop_reref(tmp, ndx);
                        else
                            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Channels E56/E57/E100/E107 not found');
                            cmd = '';
                        end
                    end
                otherwise
                    EegCallbacks.AddToListbox(app, app.listboxStdout, ...
                        sprintf('*** warning *** unknown reference ''%s''. Nothing done.', refchans));
                    set(app.pushbuttonRereference, 'BackgroundColor', EegCallbacks.ColIdle);
                    return
            end
            if ~isempty(cmd)
                tmp.history = [tmp.history char(uint8(10)) cmd];
            end
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Rereference');
            data.EEG = tmp;
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'Reref');
            guidata(hObject, data);
            set(app.pushbuttonRereference, 'BackgroundColor', EegCallbacks.ColIdle);
            pause(0.005);
        end

        % Button pushed function: pushbuttonChanlocs
        function pushbuttonChanlocs_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            set(app.pushbuttonChanlocs, 'BackgroundColor', EegCallbacks.ColBusy);
            pause(0.005)
            if ~isfield(data,'EEG')
                EegCallbacks.abortStep(app, event, 'No data available');
                app.pushbuttonChanlocs.BackgroundColor = EegCallbacks.ColIdle;
                return
            end
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Reading channel locations.');
            tmp = data.EEG;

            switch lower(P('lookup','locs'))
                case 'standard'
                    EegCallbacks.AddToListbox(app, app.listboxStdout, ' - Looking up channels in standard-10-5-cap385.elp.');
                    tmp = pop_chanedit(tmp, 'lookup', EegCallbacks.resourceFile('standard-10-5-cap385.elp'));
                    tmp.history = [tmp.history char(uint8(10)) 'pop_chanedit(tmp, "lookup","standard-10-5-cap385.elp")' ];
                case 'biosemi 64'
                    EegCallbacks.AddToListbox(app, app.listboxStdout, ' - Renaming Biosemi channels to 10/10 and looking up channels in standard-10-5-cap385.elp.');
                    labs = readtable(EegCallbacks.resourceFile('BioSemi68_labels.txt'));
                    for ch=1:length(labs.Label)
                        ndx = find(strcmp(labs.Label{ch}, {tmp.chanlocs.labels}));
                        if length(ndx)==1
                            tmp.chanlocs(ndx).labels = labs.ten10{ch};
                        elseif length(ndx)>2
                            pause;
                        else
                            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( '   channel %s not found', labs.Label{ch}));
                        end
                    end
                    tmp = pop_chanedit(tmp, 'lookup', EegCallbacks.resourceFile('standard-10-5-cap385.elp'));
                    tmp.history = [tmp.history char(uint8(10)) 'pop_chanedit(tmp, "lookup","standard-10-5-cap385.elp")' ];
                case 'biosemi 128'
                    EegCallbacks.AddToListbox(app, app.listboxStdout, ' - Copying channel locations from a 128 channel EEGLAB dataset.');
                    lookup = pop_loadset(EegCallbacks.resourceFile('BioSemi128.set'));
                    for ch=1:tmp.nbchan
                        ndx = find(strcmp(tmp.chanlocs(ch).labels, {lookup.chanlocs.labels}));
                        if length(ndx)==1
                            tmp.chanlocs(ch) = lookup.chanlocs(ndx);
                        else
                            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( '   channel %s not found', tmp.chanlocs(ch).labels));
                        end
                    end
                    tmp.history = [tmp.history char(uint8(10)) '%% added biosemi 128 channel locations' ];
                case 'biosemi 128 (10/20 names)'
                    EegCallbacks.AddToListbox(app, app.listboxStdout, ' - Lookup BioSemi 128 channel (rename known 10/10)');
                    lookup = pop_loadset(EegCallbacks.resourceFile('BioSemi128.set'));
                    for ch=1:tmp.nbchan
                        ndx = find(strcmp(tmp.chanlocs(ch).labels, {lookup.chanlocs.labels}));
                        if length(ndx)==1
                            tmp.chanlocs(ch) = lookup.chanlocs(ndx);
                        else
                            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( '   channel %s not found', tmp.chanlocs(ch).labels));
                        end
                    end
                    tmp.history = [tmp.history char(uint8(10)) '%% added biosemi 128 channel locations' ];
                    lookup1010 = readtable(EegCallbacks.resourceFile('Biosemi128_labels1010.txt'));
                    rencnt = 0;
                    for ch=1:length(tmp.chanlocs)
                        ndx = find(strcmpi(tmp.chanlocs(ch).labels, lookup1010.orig128));
                        if ~isempty(ndx) && length(ndx)==1
                            if ~isempty(lookup1010.new1010{ndx})
                                rencnt = rencnt + 1;
                            end
                            tmp.chanlocs(ch).labels = lookup1010.mixed128{ndx};
                        end
                    end
                    EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(' - Renamed %d channels to 10/10', rencnt));
                case 'tdbrain'
                    EegCallbacks.AddToListbox(app, app.listboxStdout, ' - TD Brain channel locations');
                    lookup = pop_loadset(EegCallbacks.resourceFile('TDBrainLocs.set'));
                    for ch=1:tmp.nbchan
                        ndx = find(strcmp(tmp.chanlocs(ch).labels, {lookup.chanlocs.labels}));
                        if length(ndx)==1
                            tmp.chanlocs(ch) = lookup.chanlocs(ndx);
                        else
                            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( '   channel %s not found', tmp.chanlocs(ch).labels));
                        end
                    end
                    tmp.history = [tmp.history char(uint8(10)) '%% added TDBrain channel locations' ];
                case 'brain products 128'
                    EegCallbacks.AddToListbox(app, app.listboxStdout, '- Brain Products GmbH 128 channel locations');
                    labs = readtable(EegCallbacks.resourceFile('output_montage.txt'), 'Delimiter', '\t');
                    for ch=1:length(labs.Label)
                        ndx = find(strcmp(labs.Label{ch}, {tmp.chanlocs.labels}));
                        if length(ndx)==1
                            tmp.chanlocs(ndx).labels = labs.ten10{ch};
                        elseif length(ndx)>2
                            pause;
                        else
                            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( '   channel %s not found', labs.Label{ch}));
                        end
                    end
                    tmp = pop_chanedit(tmp, 'lookup', EegCallbacks.resourceFile('standard-10-5-cap385.elp'));
                    tmp.history = [tmp.history char(uint8(10)) '%% added Brain Products 128 channel locations' ];
                case 'temple uni'
                    EegCallbacks.AddToListbox(app, app.listboxStdout, '- Temple university channel locations');
                    for ch=1:tmp.nbchan
                        s = strsplit(tmp.chanlocs(ch).labels, '-');
                        tmp.chanlocs(ch).labels = s{1};
                    end
                    tmp = pop_chanedit(tmp, 'lookup', EegCallbacks.resourceFile('standard-10-5-cap385.elp'));
                    tmp.history = [tmp.history char(uint8(10)) 'pop_chanedit(tmp, "lookup","standard-10-5-cap385.elp")' ];
                case 'egi 128'
                    EegCallbacks.AddToListbox(app, app.listboxStdout, ' - EGI 128 EEG channels lookup');
                    tmp = pop_chanedit(tmp, {'lookup', EegCallbacks.resourceFile('GSN-HydroCel-128.sfp'), 'filetype','autodetect'});
                    tmp.history = [tmp.history char(uint8(10)) '%% added EGI 128 channel locations' ];
                case 'egi 128 (10/20 names)'
                    EegCallbacks.AddToListbox(app, app.listboxStdout, ' - EGI 128 EEG channels lookup with 10/10 renaming');
                    tmp = pop_chanedit(tmp, {'lookup', EegCallbacks.resourceFile('GSN-HydroCel-129.sfp'), 'filetype','autodetect'});
                    recode = readtable(EegCallbacks.resourceFile('GSN-HydroCel-128_recode1020.txt'));
                    for ch=1:size(recode,1)
                        ndx = find(strcmpi({tmp.chanlocs.labels},recode.OldLabel{ch}));
                        if length(ndx)==1
                            tmp.chanlocs(ndx).labels = recode.NewLabel{ch};
                        elseif length(ndx)>1
                            warning('Multiple mathcing labels found?')
                        end
                    end
                    tmp.history = [tmp.history char(uint8(10)) '%% added EGI 128 channel locations (10/20)' ];
            end

            % Add a zeroed reference channel AFTER lookup/renaming, so the
            % duplicate check tests the FINAL channel names (e.g. a Cz that only
            % appears once the 10/10 renaming has run — prevents a second Cz).
            adding = P('lookup','addref');
            if ~strcmpi(adding,'none')
                if sum(strcmpi({tmp.chanlocs.labels}, adding))==0
                    EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(' - Adding %s as a flatline', adding));
                    tmp.data(end+1,:) = 0;
                    tmp.nbchan = tmp.nbchan + 1;
                    tmp.chanlocs(end+1).labels = adding;
                    % give the new channel a location from the standard montage
                    try
                        loc = readlocs(EegCallbacks.resourceFile('standard-10-5-cap385.elp'));
                        m = find(strcmpi({loc.labels}, adding), 1);
                        if ~isempty(m)
                            fn = fieldnames(loc);
                            for k = 1:numel(fn)
                                if isfield(tmp.chanlocs, fn{k})
                                    tmp.chanlocs(end).(fn{k}) = loc(m).(fn{k});
                                end
                            end
                            tmp.chanlocs(end).labels = adding;
                        else
                            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf('   (no standard location found for %s)', adding));
                        end
                    catch
                        EegCallbacks.AddToListbox(app, app.listboxStdout, '   (location lookup for added channel failed)');
                    end
                    EegCallbacks.AddToListbox(app, app.listboxStdout, ' - Removing ICA decomposition.');
                    tmp.icaweights = [];
                    tmp.icawinv = [];
                    tmp.icasphere = [];
                    tmp = eeg_checkset(tmp);
                else
                    EegCallbacks.AddToListbox(app, app.listboxStdout, 'Warning: ref channel already exists. Skipping...');
                end
            end

            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Lookup');
            data.EEG = tmp;
            guidata(hObject,data)

            EegCallbacks.AddToListbox(app, app.listboxStdout, 'First record of power values.');
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'Chanlocs');
            guidata(hObject, data);

            EegCallbacks.listboxEegProperties_Update(app, hObject)
            set(app.pushbuttonChanlocs, 'BackgroundColor', EegCallbacks.ColIdle);
        end

        % Button pushed function: pushbuttonInitialICA
        function pushbuttonInitialICA_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            set(app.pushbuttonInitialICA,'backgroundcolor',EegCallbacks.ColBusy)
            pause(0.005);
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Running intial ICA of 16 PCs');
            tmp = data.EEG;
            if (tmp.trials==1)
                tmp.data = detrend(tmp.data','constant')';
            end
            % ICLabel topoplots each IC over EEG.chanlocs(icachansind), so ICA
            % must run ONLY on channels that have scalp locations.
            goodchans = find(~cellfun(@isempty, {tmp.chanlocs.X}));
            if isempty(goodchans), goodchans = 1:tmp.nbchan; end
            npca = min(16, numel(goodchans));
            icatype = EegCallbacks.icaType(app);
            try
                if strcmpi(icatype,'jader')
                    warning('using the JADER ICA algorithm. Ncomps/ PCA not used')
                    tmp = pop_runica(tmp, 'icatype', 'jader', 'chanind', goodchans);
                else
                    tmp = pop_runica(tmp,'icatype',icatype,'pca',npca,'chanind',goodchans);
                end
            catch
                tmp = pop_runica(tmp,'icatype','runica','extended',1,'pca',npca,'chanind',goodchans);
            end
            if isempty(tmp.icaact)
                EegCallbacks.AddToListbox(app, app.listboxStdout, ' Recalculate ICA activations.');
                tmp.icaact=icaact(tmp.data(tmp.icachansind,:), tmp.icaweights*tmp.icasphere);
            end
            EegCallbacks.AddToListbox(app, app.listboxStdout, ' Get eye PCs using ICLabel.');
            tmp = pop_iclabel(tmp, 'default');
            eyelabel = FindSetNdx(tmp.etc.ic_classification.ICLabel.classes,'Eye');
            icdeselect = (tmp.etc.ic_classification.ICLabel.classifications(:,eyelabel)')' > P('eog','criterion');
            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(' Removing %d ICs',sum(icdeselect)));
            if sum(icdeselect)>0
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(' - Remove %d eye ICs',sum(icdeselect)));
                fprintf(' - Remove %d eye ICs using subtraction method\n',sum(icdeselect));
                tmp = EegCallbacks.subtractIcs(tmp, icdeselect);
                data.tabLine.removeEOG = sum(icdeselect);
                tmp.history = [tmp.history char(uint8(10)) ...
                    sprintf('%% removed %d eye ICs (',sum(icdeselect)) ...
                    sprintf('%d ', find(icdeselect)) ...
                    ')' ];
            end
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'EOG');
            data.EEG = tmp;
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'InitialEOG');
            guidata(hObject, data);
            set(app.pushbuttonInitialICA,'backgroundcolor',EegCallbacks.ColIdle)
        end

        % Button pushed function: pushbuttonAltEOG
        function pushbuttonAltEOG_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            set(app.pushbuttonAltEOG,'backgroundcolor',EegCallbacks.ColBusy)
            pause(0.005);
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Removing EOG with *EOG* channels');
            tmp = data.EEG;
            tmp.data = detrend(tmp.data','constant')';
            fprintf('\n== EOG ICA Debug Info ==\n');
            fprintf('EEG data size: [%d channels x %d timepoints]\n', size(tmp.data, 1), size(tmp.data, 2));
            fprintf('EEG data rank: %d\n', rank(double(tmp.data)));
            fprintf('Requested PCA components: 12\n');
            eogndx = FindSetNdx({tmp.chanlocs.labels}, '*eog*', 'match','pattern');
            if isempty(eogndx)
                EegCallbacks.AddToListbox(app, app.listboxStdout, 'WARNING! could not match EOG channels. Nothing performed');
                set(app.pushbuttonAltEOG,'backgroundcolor',EegCallbacks.ColIdle)
                return
            end
            icatype = EegCallbacks.icaType(app);
            try
                if strcmpi(icatype,'jader')
                    warning('using the JADER ICA algorithm. Ncomps/ PCA not used')
                    [tmp, cmd] = pop_runica(tmp, 'icatype', 'jader'); %#ok<ASGLU>
                else
                    [tmp, cmd] = pop_runica(tmp,'icatype',icatype,'pca',12); %#ok<ASGLU>
                end
            catch
                [tmp, cmd] = pop_runica(tmp,'icatype','runica','extended',1,'pca',12); %#ok<ASGLU>
            end
            if isempty(tmp.icaact)
                EegCallbacks.AddToListbox(app, app.listboxStdout, '- Recalculate ICA activations.');
                tmp.icaact = icaact(tmp.data(tmp.icachansind,:), tmp.icaweights*tmp.icasphere);
            end
            EegCallbacks.AddToListbox(app, app.listboxStdout, '- Match eye ICs to EOG channels');
            EYE = filter_fir(tmp.data(eogndx,:), tmp.srate, 0, 20, 3.0, true);
            ICs = filter_fir(tmp.icaact(:,:), tmp.srate, 0, 20, 3.0, true);
            R = abs(corr(EYE',ICs'));
            icdeselect = any(R >= P('alt eog','criterion'), 1);
            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf('- Removing %d ICs with subtraction method',sum(icdeselect)));
            if sum(icdeselect)>0
                tmp = EegCallbacks.subtractIcs(tmp, icdeselect);
                data.tabLine.removeEOG = sum(icdeselect);
                tmp.history = [tmp.history char(uint8(10)) sprintf('%% removed %d eye ICs', sum(icdeselect))];
            end
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Alt EOG');
            data.EEG = tmp;
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'AltInitialEOG');
            guidata(hObject, data);
            set(app.pushbuttonAltEOG,'backgroundcolor',EegCallbacks.ColIdle)
        end

        % Button pushed function: pushbuttonEMG
        function pushbuttonEMG_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            set(app.pushbuttonEMG,'backgroundcolor',EegCallbacks.ColBusy)
            pause(0.005);
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Clean data of muscle actvit using AAR in 40s windows.');
            tmp = data.EEG;
            ws = P('emg','winlen');
            ss = P('emg','winshift');
            tmp = pop_autobssemg(tmp, ws, ss, 'bsscca', {'eigratio', [1000000]}, ...
                'emg_psd', {'ratio', [10],'fs', [256],'femg', [15],...
                'estimator', spectrum.welch, 'range', [0  34]});
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'EMG');
            data.EEG = tmp;
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'AAR');
            guidata(hObject, data);
            set(app.pushbuttonEMG,'backgroundcolor',EegCallbacks.ColIdle)
            pause(0.005);
        end

        % Button pushed function: pushbuttonAltEMG
        function pushbuttonAltEMG_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            set(app.pushbuttonAltEMG,'backgroundcolor',EegCallbacks.ColBusy)
            pause(0.005);
            tmp = data.EEG;
            [~,~,~,~,~,int_data] = InterpolationCleaning(tmp);
            len = tmp.srate*5;
            startpnts = 1:len:(tmp.pnts-len+1);
            Forig = nan(len,tmp.nbchan,length(startpnts));
            Fint  = Forig;
            Fnew  = Forig;
            fft_fs = linspace(0,tmp.srate,len);
            fft_fs = fft_fs(1:end-1);
            Porig = nan(tmp.srate/2+1,tmp.nbchan,length(startpnts));
            Pint  = Porig;
            cnt = 0;
            for start=startpnts
                cnt=cnt+1;
                Forig(:,:,cnt) = fft(tmp.data(:, start:start+len-1)');
                Fint(:,:,cnt)  = fft(int_data(:, start:start+len-1)');
                [Porig(:,:,cnt), ~] = pfft(tmp.data(:, start:start+len-1)', tmp.srate, ones(1,tmp.srate), 0);
                [Pint(:,:,cnt), fs] = pfft(int_data(:, start:start+len-1)', tmp.srate, ones(1,tmp.srate), 0);
                ndx = fs>13&fs<35;
                [~,~,~,reg] = ttest(db(Pint(ndx,:,1)),db(Porig(ndx,:,1)));
                select = reg.tstat < -P('alt emg','type');
                data_new(:,:,cnt) = tmp.data(:, start:start+len-1); %#ok<AGROW>
                Fnew(:,:,cnt) = Forig(:,:,cnt);
                if sum(select)
                    weight = zeros(1, len);
                    weight(fft_fs>=35) = 1;
                    weight(fft_fs>=13 & fft_fs<=35) = linspace(0,1,sum(fft_fs>=13 & fft_fs<=35));
                    weight = [weight(1:end/2) 1 weight(end/2:-1:2)];
                    Fnew(:,select,cnt) = Forig(:,select,cnt) .* repmat((1-weight)',1, sum(select)) + Fint(:,select,cnt) .* repmat(weight', 1, sum(select));
                    data_new(select,:,cnt) = ifft(Fnew(:,select,cnt))'; %#ok<AGROW>
                end
            end
            tmp.data = data_new(:,:);
            tmp.pnts = size(tmp.data,2);
            tmp.xmax = (tmp.pnts-1)/tmp.srate;
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Alt EMG');
            data.EEG = tmp;
            set(app.pushbuttonAltEMG,'backgroundcolor',EegCallbacks.ColIdle)
            guidata(hObject, data)
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'AltEMG');
            guidata(hObject, data);
        end

        % Button pushed function: pushbuttonLineNoise ('Line noise')
        % Line noise removal. Method (setting cleanline/method):
        %   eeg_linenoise  sinusoids of estimated frequency regressed out of all
        %                  channels at once, in pieces of segmentlength seconds
        %   CleanLine      Tim Mullen's plugin, called as in DISCOVER-EEG:
        %                  pop_cleanline(EEG, 'linefreqs', f, 'newversion', 1)
        % Both remove the line frequency and its 2nd-4th harmonics below
        % Nyquist, or every harmonic below Nyquist with 'harmonics' on.
        function pushbuttonLineNoise_Callback(app, event)
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(p) EegParams.get(data.params, 'cleanline', p);
            say = @(varargin) EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(varargin{:}));
            set(app.pushbuttonLineNoise, 'BackgroundColor', EegCallbacks.ColBusy); pause(0.005)
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                set(app.pushbuttonLineNoise, 'BackgroundColor', EegCallbacks.ColIdle)
                return
            end
            tmp    = data.EEG;
            f0     = str2double(P('linefreq'));
            method = P('method');
            if P('harmonics')
                mult = 1:floor((tmp.srate/2 - 1) / f0);
            else
                mult = 1:4;
            end
            mult  = mult(f0*mult < tmp.srate/2);
            fstr  = strjoin(arrayfun(@(x) sprintf('%g', x), f0*mult, 'uni', 0), ', ');

            switch lower(method)
                case 'cleanline'
                    if isempty(which('pop_cleanline'))
                        say(' *** error *** CleanLine plugin not found (EEGLAB plugins: Cleanline)');
                        EegCallbacks.abortStep(app, event, 'CleanLine plugin not found.');
                        set(app.pushbuttonLineNoise, 'BackgroundColor', EegCallbacks.ColIdle)
                        return
                    end
                    % the plugin adds its helper folders only when EEGLAB builds
                    % its menus; without that (eeglab nogui, batch) it cannot run
                    if isempty(which('hlp_varargin2struct')) || isempty(which('arg_define'))
                        addpath(genpath(fileparts(which('pop_cleanline'))));
                    end
                    % a single frequency makes the new CleanLine take f..4f itself
                    freqs = f0;
                    if P('harmonics'), freqs = f0*mult; end
                    say('Removing line noise with CleanLine at %s Hz', fstr);
                    % 'sigtype' Channels: with an ICA decomposition present
                    % CleanLine defaults to components, which its new version refuses
                    [tmp, cmd] = pop_cleanline(tmp, 'linefreqs', freqs, 'newversion', 1, 'computepower', 0, ...
                        'sigtype', 'Channels', 'chanlist', 1:tmp.nbchan);
                    if isempty(cmd)
                        cmd = sprintf('EEG = pop_cleanline(EEG, ''linefreqs'', %s, ''newversion'', 1);', mat2str(freqs));
                    end
                    tmp.history = [tmp.history newline cmd];

                otherwise   % eeg_linenoise
                    L = P('segmentlength');
                    say('Removing line noise with eeg_linenoise at %s Hz, pieces of about %g s%s', ...
                        fstr, L, ifthen(P('drift'), ', frequency drift fitted', ''));
                    [tmp, info] = eeg_linenoise(tmp, 'LineFreq', f0, 'Harmonics', mult, ...
                        'SegmentLength', L, 'MinPeakDb', P('minpeak'), 'Drift', logical(P('drift')), ...
                        'Verbose', false);
                    for f = unique([info.nominal])
                        k    = [info.nominal] == f;
                        done = k & ~[info.skipped];
                        if ~any(done)
                            say('- %g Hz: no clear peak (max %.1f dB < %g dB) in any of %d pieces; nothing removed', ...
                                f, max([info(k).peakDb]), P('minpeak'), sum(k));
                            continue
                        end
                        est = [info(done).freq];
                        say('- %g Hz: removed in %d of %d pieces (peak %.1f-%.1f dB), estimated %.3f-%.3f Hz, median %.2f%% of channel variance', ...
                            f, sum(done), sum(k), min([info(done).peakDb]), max([info(done).peakDb]), ...
                            min(est), max(est), 100*median([info(done).removed]));
                    end
            end

            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Line noise');
            data.EEG = tmp;
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'Line noise');
            guidata(hObject, data);
            set(app.pushbuttonLineNoise, 'BackgroundColor', EegCallbacks.ColIdle)
        end

        % Button pushed function: pushbuttonASR
        function pushbuttonASR_Callback(app, event)
            EegCallbacks.runCleanRawdata(app, event, 'asr', 'pushbuttonASR', 'Clean', 'CleanRawdata');
        end

        % Button pushed function: pushbuttonICA
        % ICA + ICLabel, subtracting the artefact components. With 'ICA number'
        % above 1 the decomposition is repeated that many times and the most
        % representative run is kept, following DISCOVER-EEG (see
        % selectIcaRun). With 1, a single run. With the discovereeg setting on,
        % the bad segments of the kept run (clean_rawdata burst + window
        % criteria, DISCOVER-EEG step 6) are removed as well.
        function pushbuttonICA_Callback(app, event)
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            say = @(varargin) EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(varargin{:}));
            set(app.pushbuttonICA, 'BackgroundColor', EegCallbacks.ColBusy)
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                set(app.pushbuttonICA, 'BackgroundColor', EegCallbacks.ColIdle)
                return
            end
            tmp = data.EEG;
            ncomps = P('ICA','pca');
            if ncomps<10
                ncomps = 10;
                say('Too few components selected for ICA. Taking the minimum of 10.')
            end
            goodchans = find(~cellfun(@isempty, {tmp.chanlocs.X}));
            if isempty(goodchans), goodchans = 1:tmp.nbchan; end
            if ncomps>numel(goodchans)
                ncomps=numel(goodchans);
                say('Too many components selected for ICA. Taking the maximum.')
            end
            icatype = EegCallbacks.icaType(app);
            nRuns   = EegCallbacks.icaRuns(app);
            if nRuns > 1 && strcmpi(icatype, 'jader')
                say(' *** warning *** jader is deterministic: repeating it gives the same result. Running once.');
                nRuns = 1;
            end
            if nRuns > 1 && tmp.trials > 1
                say(' *** warning *** choosing between ICA runs needs continuous data. Running once.');
                nRuns = 1;
            end
            removeSegs = logical(P('ICA','discovereeg'));
            if removeSegs && tmp.trials > 1
                say(' *** warning *** bad segments can only be removed from continuous data. Skipped.');
                removeSegs = false;
            end
            scoreRuns = nRuns > 1 || removeSegs;
            say('ICA (%s, %d components, %d run%s)', icatype, ncomps, nRuns, ifthen(nRuns>1, 's', ''));
            rules = EegCallbacks.icRules(data.params);
            if isempty(rules)
                say(' *** warning *** no ICLabel class switched on: no components will be removed');
            else
                say('- remove a component when ICLabel gives %s', strjoin(cellfun(@(c, t) sprintf('%s >= %.1f', c, t), ...
                    rules(:,1)', rules(:,2)', 'uni', 0), ' or '));
            end

            % every run: decomposition + ICLabel flags; keep only the ICA fields
            runs  = struct('icaweights',{},'icasphere',{},'icawinv',{},'icachansind',{}, ...
                           'classification',{},'reject',{},'cmd',{});
            masks = [];
            if nRuns > 1
                say('- each run is scored with a clean_rawdata bad-segment pass on a COPY of the cleaned data;');
                if removeSegs
                    say('  the bad segments of the kept run are removed afterwards (DISCOVER-EEG)');
                else
                    say('  this only picks the run and removes nothing (the ASR button does the actual cleaning)');
                end
            end
            for r = 1:nRuns
                if nRuns > 1
                    fprintf('\n===== ICA run %d of %d =====\n', r, nRuns);
                end
                [one, cmd] = EegCallbacks.decomposeIca(app, tmp, icatype, ncomps, goodchans, r);
                [one, reject] = EegCallbacks.flagIcs(one, rules);
                runs(r).icaweights     = one.icaweights;
                runs(r).icasphere      = one.icasphere;
                runs(r).icawinv        = one.icawinv;
                runs(r).icachansind    = one.icachansind;
                runs(r).classification = one.etc.ic_classification;
                runs(r).reject         = reject;
                runs(r).cmd            = cmd;
                if scoreRuns
                    fprintf(['\n----- ICA run %d of %d: clean_rawdata bad-segment mask to score this run.\n' ...
                             '      Runs on a copy; the data are NOT cleaned here. -----\n'], r, nRuns);
                    masks(r,:) = EegCallbacks.badSegmentMask(EegCallbacks.subtractIcs(one, reject), ...
                        P('ICA','repburst'), P('ICA','repwindow')); %#ok<AGROW>
                    say('- run %d: %d components flagged, %.1f s would be marked bad', ...
                        r, sum(reject), sum(~masks(r,:))/tmp.srate);
                end
            end
            sel = 1;
            if nRuns > 1
                [sel, dist] = EegCallbacks.selectIcaRun(masks);
                say('- keeping run %d, closest to the average bad-segment mask (%s samples from it)', ...
                    sel, mat2str(round(dist(:)')));
            end

            % apply the chosen decomposition and subtract its flagged components
            R = runs(sel);
            tmp.icaweights  = R.icaweights;
            tmp.icasphere   = R.icasphere;
            tmp.icawinv     = R.icawinv;
            tmp.icachansind = R.icachansind;
            tmp.icaact      = [];
            tmp.etc.ic_classification = R.classification;
            tmp.history = [tmp.history newline R.cmd];
            if nRuns > 1
                tmp.history = [tmp.history newline sprintf( ...
                    '%% ICA repeated %d times, run %d kept (DISCOVER-EEG selection)', nRuns, sel)];
            end
            removendx = find(R.reject);
            fprintf('Removing components: %s\n', mat2str(removendx(:)'));
            if isempty(removendx)
                say('- no components met a removal rule');
            end
            if ~isempty(removendx)
                say('- removing %d components: %s', numel(removendx), mat2str(removendx(:)'));
                tmp = EegCallbacks.subtractIcs(tmp, R.reject);
                data.tabLine.iclabelDeselect = numel(removendx);
                tmp.history = [tmp.history newline ...
                    sprintf('%% subtracting %d components ( %s)', numel(removendx), sprintf('%d ', removendx))];
            end
            if removeSegs
                bad = ~masks(sel, :);
                if all(bad)
                    say(' *** warning *** clean_rawdata marks all data as bad; no segments removed');
                elseif ~any(bad)
                    say('- DISCOVER-EEG bad segments (burst %g SD, window %g): none found', ...
                        P('ICA','repburst'), P('ICA','repwindow'));
                else
                    d = diff([0 bad 0]);
                    regions = [find(d == 1)' find(d == -1)' - 1];
                    say('- DISCOVER-EEG bad segments (burst %g SD, window %g): removing %d segments, %.1f s (%.1f%%)', ...
                        P('ICA','repburst'), P('ICA','repwindow'), size(regions, 1), ...
                        sum(bad)/tmp.srate, 100*mean(bad));
                    tmp.icaact = [];
                    tmp = eeg_eegrej(tmp, regions);
                    tmp.history = [tmp.history newline ...
                        sprintf('EEG = eeg_eegrej(EEG, %s); %% DISCOVER-EEG bad segments of the kept ICA run', mat2str(regions))];
                end
            end
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'ICA');
            data.EEG = tmp;
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'IClabel');
            guidata(hObject, data);
            set(app.pushbuttonICA, 'BackgroundColor', EegCallbacks.ColIdle)
        end

        % Button pushed function: pushbuttonRemoveEOG
        function pushbuttonRemoveEOG_Callback(app, event)
            % 'Remove chans': removes every channel whose label matches an entry
            % in TextAreaRemoveChans (one per line, or separated by commas or
            % spaces). Entries with * or ? are wildcards (* any text, ? one
            % character) matched against the whole label; other entries must
            % equal the label. Case-insensitive, as before. The text is saved
            % and restored with the other window controls (GetUIControlData).
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            say = @(varargin) EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(varargin{:}));
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                return
            end
            patterns = {};
            if isprop(app, 'TextAreaRemoveChans') && isvalid(app.TextAreaRemoveChans)
                patterns = cellstr(app.TextAreaRemoveChans.Value);
            elseif isprop(app, 'ListBoxRemoveChans') && isvalid(app.ListBoxRemoveChans)
                patterns = cellstr(app.ListBoxRemoveChans.Items);
            end
            patterns = strsplit(strjoin(patterns, ' '), {' ', ',', ';', char(9)});
            patterns = patterns(~cellfun(@isempty, patterns));
            if isempty(patterns)
                say('Remove chans: the channel list is empty; nothing removed');
                return
            end
            tmp = data.EEG;
            ndx = EegCallbacks.matchLabels({tmp.chanlocs.labels}, patterns);
            say('Removing channels matching: %s', strjoin(patterns, ' '));
            if isempty(ndx)
                say('- no matching channels');
                return
            end
            if numel(ndx) == tmp.nbchan
                say(' *** error *** every channel matches; nothing removed');
                EegCallbacks.abortStep(app, event, 'Remove chans would remove every channel.');
                return
            end
            say('- removing %d channels: %s', numel(ndx), strjoin({tmp.chanlocs(ndx).labels}, ' '));
            [tmp, cmd] = pop_select(tmp, 'nochannel', ndx);
            tmp.history = [tmp.history newline cmd];
            if strcmpi(tmp.ref, 'average')
                say('- redo avg reference');
                [tmp, cmd] = pop_reref(tmp, []);
                tmp.history = [tmp.history newline cmd];
            end
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Remove chans');
            data.EEG = tmp;
            guidata(hObject, data);
            pause(0.01)
            EegCallbacks.trackSpectrum(app, data.EEG, 'Remove chans');
        end

        % ------------------------------------------------------------------
        % Indices of the labels matching any of the patterns. A pattern with
        % * or ? is a wildcard over the whole label, otherwise an exact match;
        % both ignore case. Blank lines are skipped.
        function ndx = matchLabels(labels, patterns)
            hit = false(1, numel(labels));
            for k = 1:numel(patterns)
                pat = strtrim(char(patterns{k}));
                if isempty(pat)
                    continue
                end
                if any(pat == '*' | pat == '?')
                    rx = ['^' regexptranslate('wildcard', pat) '$'];
                    hit = hit | ~cellfun(@isempty, regexpi(labels, rx, 'once'));
                else
                    hit = hit | strcmpi(labels, pat);
                end
            end
            ndx = find(hit);
        end

        % Button pushed function: pushbuttonRemoveNoEEG
        function pushbuttonRemoveNoEEG_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            tmp = data.EEG;
            notEEG = find(cellfun(@isempty, {tmp.chanlocs.X}));
            eog = FindSetNdx({tmp.chanlocs.labels},'*eog*','match','pattern');
            notEEG = setdiff(notEEG, eog);
            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf('Removing %d non-EEG channels', length(notEEG)))
            [tmp, cmd] = pop_select(tmp, 'nochannel', notEEG);
            tmp.history = [tmp.history char(uint8(10)) cmd];
            if strcmpi(tmp.ref, 'average')
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf('- redo avg reference'))
                [tmp, cmd] = pop_reref(tmp, []);
                tmp.history = [tmp.history char(uint8(10)) cmd];
            end
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Rm ~EEG ch');
            data.EEG = tmp;
            guidata(hObject, data);
            pause(0.01)
            EegCallbacks.trackSpectrum(app, data.EEG, 'Remove ~EEG');
        end

        % Button pushed function: pushbuttonMaskEvt
        % Mask every event type down to its low 8 bits, dropping the events
        % that mask to zero. Was a button of its own; it is now the 'mask'
        % setting of Open and runs straight after a file is read.
        function EEG = maskEvents(app, listboxStdout, EEG)
            EegCallbacks.AddToListbox(app, listboxStdout, '- masking all events types to 0-255 (numeric or not)');
            tmp = EEG;
            keep = true(1,length(tmp.event));
            for e=1:length(tmp.event)
                ev = tmp.event(e).type;
                if isnumeric(ev)
                    if ev==round(ev)
                        val = bitand(ev, 255);
                        if val
                            ev = val;
                        else
                            keep(e) = false;
                        end
                    end
                elseif ~isnan(str2double(ev))
                    val = bitand(str2double(ev), 255);
                    if val
                        ev = num2str(bitand(str2double(ev), 255));
                    else
                        keep(e) = false;
                    end
                end
                tmp.event(e).type = ev;
            end
            tmp.event = tmp.event(keep);
            EEG = tmp;
        end

        % Button pushed function: pushbuttonBadChans
        function pushbuttonBadChans_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            data.pushbuttonBadChans.BackgroundColor = [.3 .6 .3];
            pause(0.005)
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                data.pushbuttonBadChans.BackgroundColor = [.6 1 .6];
                return
            end
            savedEEG = data.EEG;
            h = figBadChansModal(app.eeg_workflow);
            uiwait(h);
            % NOTE: original does not re-fetch guidata here (kept faithful).
            if data.EEG.nbchan~=savedEEG.nbchan
                removedndx = ~ismember({savedEEG.chanlocs.labels}, {data.EEG.chanlocs.labels});
                data = EegCallbacks.pushUndo(app, event, data, savedEEG, sprintf('Bad channels removed (%d)', sum(removedndx)));
                data.EEG.history = [data.EEG.history '\npop_select(EEG,{' sprintf('%s',savedEEG.chanlocs(removedndx).labels)];
                guidata(hObject, data);
            end
            data.pushbuttonBadChans.BackgroundColor = [.9 .8 .6];
        end

        % Button pushed function: pushbuttonFlatPeriods
        function pushbuttonFlatPeriods_Callback(app, event)
            EegCallbacks.runPeriods(app, event, 'flat periods', ...
                'Removing flatline periods (epoch SD below a threshold).', ...
                @(EEG, chans, ep, P) EegPeriods.maskFlat(EEG, chans, ep, P('sd')), ...
                'FlatPeriods');
        end

        % Button pushed function: pushbuttonSleepTheta
        function pushbuttonSleepTheta_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Removing excessive theta periods indicating sleep')
            tmp = data.EEG;
            cWinSize = tmp.srate*2;
            cWinSec = cWinSize / tmp.srate;
            cWinSecMove = cWinSec / 2;
            EOG = FindSetNdx({tmp.chanlocs.labels},'*eog*','match','pattern');
            EEG = find(~cellfun(@isempty, {tmp.chanlocs.X}));
            EEG = setdiff(EEG,EOG);
            [~,fs, allP] = pfft(tmp.data(EEG,:)', tmp.srate, ones(1, cWinSize), .5);
            Th = mean(log(squeeze(mean(allP(fs>3&fs<=6, :, :)))'), 2);
            Al = mean(log(squeeze(mean(allP(fs>7&fs<=13, :, :)))'), 2);
            H = nan(size(Th));
            H(:) = (Th(:)./Al(:))>4.0;
            starts = 0:cWinSecMove:tmp.xmax;
            row=0; times=[]; state=0;
            for s=1:length(H)
                if H(s) && state==0
                    state = 1; row = row+1; times(row,1) = starts(s); %#ok<AGROW>
                elseif ~H(s) && state==1
                    state = 0; times(row,2) = starts(s); %#ok<AGROW>
                end
            end
            if state==1 && times(row,2) == 0
                times(row,2) = tmp.xmax;
            end
            for row=1:size(times,1)-1
                if times(row+1,1)<=times(row,2)+1
                    times(row+1,1) = times(row+1,1)-1;
                end
            end
            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf('- Removing %d periods', size(times,1)))
            [tmp, cmd] = pop_select(tmp, 'notime', times);
            tmp.history = [tmp.history char(uint8(10)) cmd];
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Clear flatline periods');
            data.EEG = tmp;
            guidata(hObject, data)
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'SleepTheta');
            guidata(hObject, data);
        end

        % Button pushed function: pushbuttonRemoveTask
        % Removes the task periods: every non-boundary event with the seconds
        % before and after it.
        function pushbuttonRemoveTask_Callback(app, event)
            EegCallbacks.removeEventPeriods(app, event, true);
        end

        % Button pushed function: pushbuttonRemoveResting
        % Removes everything that is not a task period (rest, instruction
        % reading), keeping each non-boundary event with the seconds around it.
        function pushbuttonRemoveResting_Callback(app, event)
            EegCallbacks.removeEventPeriods(app, event, false);
        end

        % ------------------------------------------------------------------
        % Shared logic of Remove task / Remove no-task. The task periods are
        % [latency - before, latency + after] around every non-boundary event.
        % removeTask true deletes those periods, false deletes everything else.
        % Interactively the seconds are asked (pre-filled with the last values)
        % and stored as settings 'task periods' before/after; batch uses the
        % stored values without asking.
        function removeEventPeriods(app, event, removeTask)
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            say = @(varargin) EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(varargin{:}));
            what = ifthen(removeTask, 'task', 'no-task');
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                return
            end
            before = EegParams.get(data.params, 'task periods', 'before');
            after  = EegParams.get(data.params, 'task periods', 'after');
            if ~EegCallbacks.isBatchEvent(event)
                answer = inputdlg({'Seconds before each event', 'Seconds after each event'}, ...
                    sprintf('Remove %s periods', what), [1 35], ...
                    {sprintf('%g', before), sprintf('%g', after)});
                EegCallbacks.bringToFront(app.eeg_workflow);
                if isempty(answer)
                    say('User cancelled');
                    return
                end
                before = str2double(answer{1});
                after  = str2double(answer{2});
                if isnan(before) || isnan(after) || before < 0 || after < 0
                    say(' *** error *** give non-negative numbers of seconds');
                    return
                end
                data.params = EegParams.set(data.params, 'task periods', 'before', before);
                data.params = EegParams.set(data.params, 'task periods', 'after', after);
                EegParams.save(data.params, data.SETTINGSDIR);
            end

            tmp = data.EEG;
            if tmp.trials > 1
                EegCallbacks.abortStep(app, event, 'Task periods can only be removed from continuous data.');
                return
            end
            say('Removing %s periods (events -%g s to +%g s)', what, before, after);
            task = false(1, tmp.pnts);
            nEvents = 0;
            for e = 1:length(tmp.event)
                if ~strcmpi(tmp.event(e).type, 'boundary')
                    s1 = max(1, floor(tmp.event(e).latency - tmp.srate*before));
                    s2 = min(tmp.pnts, ceil(tmp.event(e).latency + tmp.srate*after));
                    if s2 >= s1
                        task(s1:s2) = true;
                        nEvents = nEvents + 1;
                    end
                end
            end
            if nEvents == 0
                say('- no (non-boundary) events found; nothing removed');
                return
            end
            removeMask = ifthen(removeTask, task, ~task);
            if ~any(removeMask)
                say('- nothing to remove');
                return
            end
            if all(removeMask)
                say(' *** error *** this would remove all data; nothing removed');
                EegCallbacks.abortStep(app, event, sprintf('Remove %s would remove all data.', what));
                return
            end
            d = diff([false removeMask false]);
            regions = [find(d == 1)' find(d == -1)' - 1];
            say('- %d events; removing %d periods, %.1f s', nEvents, size(regions,1), sum(removeMask)/tmp.srate);
            [tmp, cmd] = pop_select(tmp, 'nopoint', regions);
            tmp.history = [tmp.history newline cmd];
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, EegCallbacks.callerStepName());
            data.EEG = tmp;
            guidata(hObject, data);
            EegCallbacks.trackSpectrum(app, data.EEG, sprintf('Remove %s', what));
        end

        % Button pushed function: pushbuttonReview
        function pushbuttonReview_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                app.pbView.BackgroundColor = [1 .6 .6];
                return
            end
            % 
            [tmp, cmd] = reviewwindow(data.EEG);

            %global TMPREJ %#ok<GVMIS>
            %TMPREJ = [];
            %pop_eegplot(data.EEG, 1, 1, 0, [], 'title', 'Scroll EEG', ...
            %    'command', 'global TMPREJ', 'winlength', 12,...
            %    'spacing', 50, 'submean', 'on', 'dispchans', min(64,data.EEG.nbchan));
            %pause(.5);
            %uiwait(gcf);
            %pause(.5)
            %[data.EEG, cmd] = eeg_eegrej(data.EEG, eegplot2event(TMPREJ));
            if ~isempty(tmp)
                data.EEG = tmp;
                data.EEG.history = [data.EEG.history char(uint8(10)) cmd];
                guidata(hObject, data);
                data.EEG = EegCallbacks.recordPower(app, data.EEG, 'Review');
                guidata(hObject, data);
            else
                EegCallbacks.AddToListbox(app, app.listboxStdout, '*** review cancelled ***');
            end
        end

        % Button pushed function: pushbuttonImputeAll
        function pushbuttonImputeAll_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Replacing all channels with imputed version');
            tmp = data.EEG;
            % located channels: spherical spline from all other channels, all
            % at once (same result as InterpolationReplace's pop_interp loop)
            located = find(~cellfun(@isempty, {tmp.chanlocs.X}));
            tmp.data(located, :, :) = eeg_fastinterpolate(tmp, located);
            % channels without a location: statistical estimate from the
            % others, as InterpolationReplace does
            for ch = setdiff(1:tmp.nbchan, located)
                others = setdiff(1:tmp.nbchan, ch);
                [~, b] = EigenRotatePlot(double(data.EEG.data(others, :))', 10, 'varimax', false);
                mdl = fitlm(b, double(data.EEG.data(ch, :))');
                tmp.data(ch, :) = predict(mdl);
            end
            tmp.icaact = [];
            tmp.history = [tmp.history newline '% eeg_fastinterpolate: all channels replaced by their interpolation to remove local artefacts'];
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Impute All');
            data.EEG = tmp;
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'ImputeAll');
            guidata(hObject, data);
        end

        % Button pushed function: pushbuttonEMGPeriods
        function pushbuttonEMGPeriods_Callback(app, event)
            data = guidata(app.eeg_workflow);
            if isfield(data, 'params')
                P = @(p) EegParams.get(data.params, 'EMG periods', p);
                if ~P('dBcriterion') && ~P('doPsdSlope') && ~P('zP') && ~P('zFDR')
                    EegCallbacks.AddToListbox(app, app.listboxStdout, ...
                        'EMG periods: the dB, slope and z-score tests are all off; nothing to do');
                    return
                end
                tests = {};
                if P('dBcriterion'), tests{end+1} = sprintf('dB < %g', P('criterion')); end
                if P('doPsdSlope')
                    tests{end+1} = sprintf('slope > %g dB/Hz and slope-interp %s %g dB/Hz', ...
                        P('slopeCriterionAbsolute'), ifthen(P('slopeCriterionInterpolationDelta') >= 0, '>', '<'), ...
                        P('slopeCriterionInterpolationDelta'));
                end
                if P('zP'),   tests{end+1} = sprintf('band power z-score p < %g', P('zPvalue')); end
                if P('zFDR'), tests{end+1} = sprintf('band power z-score significant (FDR q < %g)', P('zFDRalpha')); end
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                    'EMG periods, %g-%g Hz: %s', P('highpass'), P('lowpass'), strjoin(tests, ' OR ')));
            end
            EegCallbacks.runPeriods(app, event, 'EMG periods', ...
                'Removing muscle artefact periods (band power and spectral slope vs the interpolated signal).', ...
                @(EEG, chans, ep, P) EegPeriods.maskEMG(EEG, chans, ep, ...
                    P('highpass'), P('lowpass'), P('criterion'), P('dBcriterion'), ...
                    P('doPsdSlope'), P('slopeCriterionAbsolute'), P('slopeCriterionInterpolationDelta'), ...
                    P('zP'), P('zPvalue'), P('zFDR'), P('zFDRalpha')), ...
                'EMGPeriods');
        end

        % Button pushed function: pushbuttonExcessivePeriods
        function pushbuttonExcessivePeriodsButtonPushed(app, event)
            EegCallbacks.runPeriods(app, event, 'excessive periods', ...
                'Removing excessive periods (outlying epoch SD or amplitude).', ...
                @(EEG, chans, ep, P) EegPeriods.maskExcessive(EEG, chans, ep, ...
                    P('alpha'), P('maxamp')), ...
                'ExcessivePeriods');
        end

        % Button pushed function: pushbuttonRmEyesOpen
        function pushbuttonRmEyesOpenButtonPushed(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            tmp = data.EEG;
            labels = {tmp.chanlocs.labels};
            eoglabels = {'VEOG', 'EOGV', 'Fp1', 'Fp2', 'Fpz', 'F1','Fz','F2'};
            eog = min(find(ismember(eoglabels, labels))); %#ok<MXFND>
            eoglabel = eoglabels{eog};
            blinkch = find(ismember(labels, eoglabel));
            tmp_filt = pop_eegfiltnew(pop_select(tmp, 'channel', blinkch), 'locutoff', 1, 'hicutoff', 10, 'plotfreqz', 0);
            x = double(tmp_filt.data(1,:));
            xa = abs(x);
            med = median(xa);
            madv = median(abs(xa - med));
            sigma_robust = 1.4826 * madv;
            thr = med + 28*sigma_robust;
            minDist = round(0.25 * tmp.srate);
            [pks, locs] = findpeaks(xa, 'MinPeakHeight', thr, 'MinPeakDistance', minDist); %#ok<ASGLU>
            n0 = numel(tmp.event);
            for k = 1:numel(locs)
                tmp.event(n0+k).type    = 'blink';
                tmp.event(n0+k).latency = locs(k);
                tmp.event(n0+k).duration = 0;
            end
            if numel(tmp.event)>0
                [~, idx] = sort([tmp.event.latency]);
                tmp.event = tmp.event(idx);
                tmp = eeg_checkset(tmp, 'eventconsistency');
                EegCallbacks.AddToListbox(app, app.listboxStdout, ...
                    sprintf('Added %d blink events from channel %s (thr=%.2f).\n', numel(locs), eoglabel, thr));
                winSec = 4;
                winSamp = round(winSec * tmp.srate);
                blinkIdx = find(strcmpi({tmp.event.type}, 'blink'));
                if isempty(blinkIdx)
                    warning('No blink events found.');
                    return;
                end
                lat = round([tmp.event(blinkIdx).latency]);
                rej = [lat - winSamp; lat + winSamp]';
                rej(:,1) = max(rej(:,1), 1);
                rej(:,2) = min(rej(:,2), tmp.pnts);
                rej = sortrows(rej, 1);
                merged = rej(1,:);
                for i = 2:size(rej,1)
                    if rej(i,1) <= merged(end,2) + 1
                        merged(end,2) = max(merged(end,2), rej(i,2));
                    else
                        merged = [merged; rej(i,:)]; %#ok<AGROW>
                    end
                end
                tmp = eeg_eegrej(tmp, merged);
                tmp = eeg_checkset(tmp, 'eventconsistency');
                fprintf('Removed %d merged rejection windows (±%gs around each blink).\n', size(merged,1), winSec);
            else
                fprintf('Removed NO blink windows.\n');
            end
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Remove blink based eyes-open periods');
            data.EEG = tmp;
            guidata(hObject, data);
            EegCallbacks.trackSpectrum(app, data.EEG, 'Remove eyes open');
        end

        % Button pushed function: pushbuttonSave
        function pushbuttonSave_Callback(app, event) %#ok<INUSD>
            data = guidata(app.eeg_workflow);
            tmp = data.EEG;
            h = figSaveModal(tmp, tmp.filename, data.fontsize, data.SETTINGSDIR);
            uiwait(h);
        end

        % Button pushed function: pbView
        function pbView_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                app.pbView.BackgroundColor = [1 .6 .6];
                return
            end
            screensize = get(groot, 'Screensize' );
            tmp = data.EEG;
            eegplot(tmp.data,'srate',tmp.srate,'eloc_file',tmp.chanlocs,'spacing',50,...
                'limits',[tmp.xmin tmp.xmax],'winlength',12,'position',screensize,...
                'events',tmp.event, 'submean', 'on', 'dispchans', min(64,tmp.nbchan));
            guidata(hObject,data);
        end

        % Button pushed function: pushbuttonOverlay
        function pushbuttonOverlay_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                app.pbView.BackgroundColor = [1 .6 .6];
                return
            end
            if ~isfield(data,'EEG') || isempty(data.EEG) || data.EEG.nbchan==0
                EegCallbacks.AddToListbox(app, app.listboxStdout, '*** Warning *** no data available');
                return
            end
            if ~isfield(data,'Stack') || isempty(data.Stack)
                EegCallbacks.AddToListbox(app, app.listboxStdout, '*** Warning *** no comparison data available');
                return
            end
            screensize = get(groot, 'Screensize' );
            tmp1 = data.EEG;
            tmp2 = data.Stack{length(data.Stack)};
            if tmp1.nbchan==tmp2.nbchan && tmp1.pnts==tmp2.pnts
                eegplot(tmp2.data,'srate',tmp2.srate,'eloc_file',tmp2.chanlocs,'spacing',50,...
                    'limits',[tmp2.xmin tmp2.xmax],'winlength',12,'position',screensize,...
                    'events',tmp2.event, 'data2', tmp1.data, 'dispchans', min(tmp1.nbchan, 64));
            else
                EegCallbacks.AddToListbox(app, app.listboxStdout, '*** Warning *** data and comparison data are incompatible to plot together');
                u = union({tmp2.chanlocs.labels}, {tmp1.chanlocs.labels});
                tmp1 = pop_select(tmp1, 'channel', u);
                tmp2 = pop_select(tmp2, 'channel', u);
                if tmp1.srate~=tmp2.srate
                    newrate = min(tmp1.srate, tmp2.srate);
                    tmp1 = pop_resample(tmp1,newrate);
                    tmp2 = pop_resample(tmp2,newrate);
                end
                if tmp1.pnts ~= tmp2.pnts
                    newpnts = min(tmp1.pnts, tmp2.pnts);
                    tmp1 = pop_select(tmp1, 'point', [1 newpnts]);
                    tmp2 = pop_select(tmp2, 'point', [1 newpnts]);
                    warning('EEG data does not match in size!')
                end
                if tmp1.nbchan ~= tmp2.nbchan
                    both = intersect({tmp1.chanlocs.labels},{tmp2.chanlocs.labels});
                    tmp1 = pop_select(tmp1, 'channel', both);
                    tmp2 = pop_select(tmp2, 'channel', both);
                end
                eegplot(tmp2.data,'srate',tmp2.srate,'eloc_file',tmp2.chanlocs,'spacing',50,...
                    'limits',[tmp2.xmin tmp2.xmax],'winlength',12,'position',screensize,...
                    'events',tmp2.event, 'data2', tmp1.data, 'dispchans', min(tmp2.nbchan, 64));
            end
            guidata(hObject,data);
        end

        % Button pushed function: pbCompERP
        function pbCompERP_Callback(app, event) %#ok<INUSD>
            data = guidata(app.eeg_workflow);
            tmp = data.EEG;
            if size(tmp.data,3) == 1
                error('Data must be in epochs')
            end
            if iscell(tmp.epoch(1).eventtype)
                evtlist = arrayfun(@(x)x.eventtype(1),tmp.epoch,'uniformoutput',false);
                try
                    evtlist = [evtlist{:}];
                catch
                end
            else
                evtlist = {tmp.epoch.eventtype};
            end
            uniqueevt = unique(evtlist);
            uniqueevt = uniqueevt(~strcmp(uniqueevt,'boundary'));
            figure('pos',[50,50,800,600]);
            nchans = size(tmp.icawinv,2);
            nrows = floor(sqrt(nchans));
            ncols = ceil(sqrt(nchans));
            if nrows*ncols<nchans
                ncols=ncols+1;
            end
            colors = {[0 0 0]; [.8 .2 .2]; [.3 .3 .8]; [.2 .7 .2]; [1 0 0]; [0 0 .8]; [0 .7 0]};
            times = tmp.xmin:(1/tmp.srate):tmp.xmax;
            if ~isfield(tmp,'icaact') || isempty(tmp.icaact)
                tmp.icaact = icaact(tmp.data, tmp.icaweights*tmp.icasphere);
            end
            epochdata = reshape(tmp.icaact,size(tmp.icaact,1),tmp.pnts,[]);
            ERP = [];
            for e=1:length(uniqueevt)
                ERP(:,:,e) = mean(epochdata(:,:,strcmpi(evtlist,uniqueevt{e})),3); %#ok<AGROW>
            end
            for ch=1:size(epochdata,1)
                subplot(nrows,ncols,ch);
                hold on
                for e=1:length(uniqueevt)
                    plot(times, ERP(ch,:,e),'-','color',colors{e});
                end
                ylim([min(ERP(:)) max(ERP(:))])
                xlim([tmp.xmin tmp.xmax])
                box off
                title(sprintf('comp%d',ch))
                pos = get(gca,'pos');
                axes('pos',[pos(1)+.65*pos(3) pos(2) pos(3)*.4 pos(4)*.4])
                topoplot(tmp.icawinv(:,ch),tmp.chanlocs,'electrodes','off')
                box off
                axis off
                drawnow;
            end
        end

        % Button pushed function: pbDFA
        function pbDFA_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            tmp = data.EEG;
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Plotting DFA of alpha oscillations, EEG channels');
            lo = 7.0;
            hi = 13;
            deselect = contains({tmp.chanlocs.labels}, 'eog', 'ignorecase', true) | ...
                strcmpi({tmp.chanlocs.labels}, 'A1') | ...
                strcmpi({tmp.chanlocs.labels}, 'A2') | ...
                strcmpi({tmp.chanlocs.labels}, 'M1') | ...
                strcmpi({tmp.chanlocs.labels}, 'M2') | ...
                isempty([tmp.chanlocs.X]);
            exp = dfa(abs(hilbert(filter_fir(tmp.data(~deselect,:), tmp.srate, lo, hi, 3.0, true)')), tmp.srate);
            if ~isfield(data, 'figDFA') || isempty(data.figDFA) || ~isstruct(data.figDFA) || ~isfield(data.figDFA, 'fig') || ~isvalid(data.figDFA.fig)
                EegCallbacks.AddToListbox(app, app.listboxStdout, '- creating figure');
                if ~isfield(data, 'figDFA')
                    data.figDFA = struct;
                end
                fig = figure('Name', 'plot DFA topoplot', 'Tag', 'figDFA', 'NumberTitle', 'off', ...
                             'Position', [50, 50, 300, 250]);
                ax = axes('Parent', fig);
                data.figDFA.fig = fig;
                data.figDFA.ax = ax;
            else
                EegCallbacks.AddToListbox(app, app.listboxStdout, '- activating figure');
                fig = data.figDFA.fig;
                figure(fig);
            end
            tmp.nbchan = size(tmp.data,1);
            clf;
            topoplot(exp, tmp.chanlocs(~deselect), 'maplimits', [.5 .8]);
            colorbar;
            guidata(hObject, data);
        end

        % Button pushed function: pbERP
        function pbERP_Callback(app, event) %#ok<INUSD>
            data = guidata(app.eeg_workflow);
            tmp = data.EEG;
            if size(tmp.data,3) == 1
                error('Data must be in epochs')
            end
            if isfield(tmp, 'eventlist')
                uniqueevt = tmp.eventlist;
                evtlist = {};
                for ep=1:length(tmp.epoch)
                    try
                        ndx = find(ismember(tmp.eventlist, tmp.epoch(ep).eventtype));
                    catch
                        ndx = [];
                    end
                    if ~isempty(ndx)
                        evtlist{ep} = tmp.eventlist{ndx}; %#ok<AGROW>
                    else
                        evtlist{ep} = []; %#ok<AGROW>
                    end
                end
            else
                if iscell(tmp.epoch(1).eventtype)
                    evtlist = arrayfun(@(x)x.eventtype(1),tmp.epoch,'uniformoutput',false);
                    evtlist = [evtlist{:}];
                    uniqueevt = unique(evtlist);
                else
                    evtlist = {tmp.epoch.eventtype};
                    uniqueevt = unique(evtlist);
                end
                uniqueevt = uniqueevt(~strcmp(uniqueevt,'boundary'));
            end
            figure('pos',[50,50,800,600]);
            colors = {[0 0 0]; [.8 .2 .2]; [.3 .3 .8]; [.2 .7 .2]; [1 0 0]; [0 0 .8]; [0 .7 0]};
            times = tmp.xmin:(1/tmp.srate):tmp.xmax;
            Xs = [tmp.chanlocs.Y];
            Ys = [tmp.chanlocs.X];
            Zs = [tmp.chanlocs.Z];
            ERP = [];
            for e=1:length(uniqueevt)
                ERP(:,:,e) = mean(tmp.data(:,:,strcmpi(evtlist,uniqueevt{e})),3); %#ok<AGROW>
            end
            for ch=1:size(tmp.data,1)
                pos = [.2+.6*(Xs(ch)-min(Xs))/range(Xs) .2+.6*(Ys(ch)-min(Ys))/range(Ys) .08 .08];
                if Zs(ch)<20
                    pos(1:2) = (pos(1:2)-.5)*(1.1+abs(.3*((Zs(ch)-20)/range([20 min(Zs)]))))+.5;
                end
                axes('position',pos);
                hold on
                for e=1:length(uniqueevt)
                    plot(times, ERP(ch,:,e),'-','color',colors{e});
                end
                ylim([min(ERP(:)) max(ERP(:))])
                xlim([tmp.xmin tmp.xmax])
                box off
                axis off
                text(min(xlim)+range(xlim)*.1,min(ylim)+range(ylim)*.1,tmp.chanlocs(ch).labels)
                drawnow;
            end
        end

        % Button pushed function: pbViewPSD
        function pbViewPSD_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            tmp = data.EEG;
            chanlocs = tmp.chanlocs;
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Plotting PSD');
            lo = P('filter','low');
            hi = P('filter','high');
            [P,fs] = EegCallbacks.calc_PSD(app, tmp);
            if ~isfield(data, 'figPSD') || isempty(data.figPSD) || ~isstruct(data.figPSD) || ~isfield(data.figPSD, 'fig') || ~isvalid(data.figPSD.fig)
                EegCallbacks.AddToListbox(app, app.listboxStdout, '- creating figure');
                if ~isfield(data, 'figPSD')
                    data.figPSD = struct;
                end
                fig = figure('Name', 'plot Power Spectrum', 'Tag', 'figPSD', 'NumberTitle', 'off', ...
                             'Position', [100, 100, 700, 500]);
                frame = uipanel('Parent', fig, 'Title', 'Power spectrum', 'Tag', 'framePSD', ...
                                'Position', [0.005, 0.1, .98, 0.9]);
                ax = axes('Parent', frame, 'Position', [0.07, 0.1, .91, 0.88], 'Tag', 'axesPSD');
                checkbox = uicontrol('Style', 'checkbox', 'Parent', fig, 'Tag', 'cbSummary', ...
                                     'String', 'Summarize into regions', 'Fontsize', data.fontsize, ...
                                     'Units', 'normalized', ...
                                     'Position', [0.05, 0.03, 0.4, 0.05], ...
                                     'Callback', @(src, ev) EegCallbacks.checkbox_callback(app, src, ax, fs, P, chanlocs, data.fontsize+1, lo, hi));
                data.figPSD.fig = fig;
                data.figPSD.ax = ax;
                data.figPSD.checkbox = checkbox;
                set(ax,'xlim',[0 min(hi,45)])
            else
                EegCallbacks.AddToListbox(app, app.listboxStdout, '- activating figure');
                fig = data.figPSD.fig;
                ax = data.figPSD.ax;
                checkbox = data.figPSD.checkbox;
                figure(fig);
            end
            EegCallbacks.checkbox_callback(app, checkbox, ax, fs, P, chanlocs, data.fontsize+1, lo, hi);
            guidata(hObject, data);
        end

        % Button pushed function: pbOverlayPSD
        function pbOverlayPSDButtonPushed(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                app.pbView.BackgroundColor = [1 .6 .6];
                return
            end
            if ~isfield(data,'EEG') || isempty(data.EEG) || data.EEG.nbchan==0
                EegCallbacks.AddToListbox(app, app.listboxStdout, '*** Warning *** no data available');
                return
            end
            if ~isfield(data,'Stack') || isempty(data.Stack)
                EegCallbacks.AddToListbox(app, app.listboxStdout, '*** Warning *** no comparison data available');
                return
            end
            EegCallbacks.AddToListbox(app, app.listboxStdout, 'Plotting PSDs of current and previous');
            tmp = data.EEG;
            tmp2 = data.Stack{length(data.Stack)};
            chanlocs = tmp.chanlocs;
            lo = P('filter','low');
            hi = P('filter','high');
            [P,fs]   = EegCallbacks.calc_PSD(app, tmp);
            [P2,fs2] = EegCallbacks.calc_PSD(app, tmp2);
            % compare the channels both versions have, by label; the previous
            % version usually has more (channels removed since)
            [P, P2, fs, chanlocs, info] = EegCallbacks.matchPSDs(P, fs, tmp.chanlocs, P2, fs2, tmp2.chanlocs);
            if isempty(chanlocs)
                EegCallbacks.AddToListbox(app, app.listboxStdout, '*** Warning *** the current and previous data have no channel labels in common');
                return
            end
            EegCallbacks.AddToListbox(app, app.listboxStdout, info);
            if ~isfield(data, 'figPSD') || isempty(data.figPSD) || ~isstruct(data.figPSD) || ~isfield(data.figPSD, 'fig') || ~isvalid(data.figPSD.fig)
                EegCallbacks.AddToListbox(app, app.listboxStdout, '- creating figure');
                if ~isfield(data, 'figPSD')
                    data.figPSD = struct;
                end
                fig = figure('Name', 'plot Power Spectrum', 'Tag', 'figPSD', 'NumberTitle', 'off', ...
                             'Position', [100, 100, 700, 500]);
                frame = uipanel('Parent', fig, 'Title', 'Power spectrum', 'Tag', 'framePSD', ...
                                'Position', [0.005, 0.1, .98, 0.9]);
                ax = axes('Parent', frame, 'Position', [0.07, 0.1, .91, 0.88], 'Tag', 'axesPSD');
                checkbox = uicontrol('Style', 'checkbox', 'Parent', fig, 'Tag', 'cbSummary', ...
                                     'String', 'Summarize into regions', 'Fontsize', data.fontsize, ...
                                     'Units', 'normalized', ...
                                     'Position', [0.05, 0.03, 0.4, 0.05], ...
                                     'Callback', []);
                data.figPSD.fig = fig;
                data.figPSD.ax = ax;
                data.figPSD.checkbox = checkbox;
                data.figPSD.chanlocs = chanlocs;
                set(ax,'xlim',[0 min(hi,45)])
            else
                EegCallbacks.AddToListbox(app, app.listboxStdout, '- activating figure');
                fig = data.figPSD.fig;
                ax = data.figPSD.ax;
                checkbox = data.figPSD.checkbox;
                try
                    chanlocs = data.figPSD.chanlocs; %#ok<NASGU>
                catch
                    data.figPSD.chanlocs = tmp.chanlocs;
                    chanlocs = tmp.chanlocs; %#ok<NASGU>
                end
                figure(fig);
            end
            % (re)bind the checkbox to THIS comparison; a figure that already
            % existed would otherwise keep toggling the spectra it was made with
            data.figPSD.chanlocs = chanlocs;
            checkbox.Callback = @(src, ev) EegCallbacks.checkbox_callback_overlay(app, src, ev, ax, fs, {P, P2}, data.fontsize+1, lo, hi, chanlocs);
            EegCallbacks.checkbox_callback_overlay(app, checkbox, event, ax, fs, {P, P2}, data.fontsize+1, lo, hi, chanlocs);
            guidata(hObject, data);
        end

        % PSD checkbox callback (bound to a uicontrol inside the PSD figure)
        function checkbox_callback(app, hObject, ax, fs, P, chanlocs, fontsize, lo, hi) %#ok<INUSL>
            ndx = fs>lo & fs<hi;
            if get(hObject,'Value')==0
                cla
                traces = plot(ax, fs(ndx), 10*log10(P(ndx,:)));
                leg = {chanlocs.labels};
            else
                numlabels = {'theta','radius','X','Y','Z','sph_theta','sph_phi','sph_radius'};
                tab = struct2table(chanlocs);
                for lab=1:length(numlabels)
                    if ismember(numlabels{lab}, tab.Properties.VariableNames)
                        if iscell(tab.(numlabels{lab}))
                            values = cellfun(@(x)ifthen(isempty(x), nan, double(x)), tab.(numlabels{lab}));
                            tab.(lab) = values;
                        end
                    end
                end
                relX = tab.X ./ sqrt(tab.X.^2+tab.Y.^2+tab.Z.^2);
                relY = tab.Y ./ sqrt(tab.X.^2+tab.Y.^2+tab.Z.^2);
                ant    = relX>=-1E-5;
                post   = ~ant;
                medial = abs(relY)<.41;
                left   = relY>=.41;
                right  = relY<=.41;
                regP = nan(size(P,1),6);
                regP(:,1) = mean(P(:,ant & left),2);
                regP(:,2) = mean(P(:,ant & medial),2);
                regP(:,3) = mean(P(:,ant & right),2);
                regP(:,4) = mean(P(:,post & left),2);
                regP(:,5) = mean(P(:,post & medial),2);
                regP(:,6) = mean(P(:,post & right),2);
                traces = plot(ax, fs(ndx), 10*log10(regP(ndx,:)));
                leg = {'ant left','ant medial','ant right','post left','post medial','post right'};
            end
            set(ax, 'fontsize', fontsize+2)
            xlabel('frequency (Hz)')
            ylabel('Power({\mu}V^2/Hz)')
            patch([8 12 12 8], repelem(ylim, 2), [.5 .5 .5], ...
                'facealpha', .2, 'edgecolor','none')
            line([1 1], ylim, 'linestyle', ':')
            line([4 4], ylim, 'linestyle', ':')
            line([10 10], ylim, 'linestyle', ':')
            line([19 19], ylim, 'linestyle', ':')
            line([30 30], ylim, 'linestyle', ':')
            legend(traces, leg, 'location','northeast')
        end

        % PSD overlay checkbox callback (bound to a uicontrol inside the figure)
        function checkbox_callback_overlay(app, hObject, event, ax, fs, Ps, fontsize, lo, hi, chanlocs) %#ok<INUSD>
            if nargin < 10
                data = guidata(app.eeg_workflow);
                chanlocs = data.EEG.chanlocs;
            end
            ndx = fs>lo & fs<hi;
            P = Ps{1};
            P2= Ps{2};
            if any(size(P)~=size(P2))
                errordlg2('Only available for data of the same size.','Error')
                return
            end
            cla;
            if get(hObject,'Value')==0
                plot(ax, fs(ndx), 10*log10(P(ndx,:)));
                hold on
                ax.ColorOrderIndex = 1;
                plot(ax, fs(ndx)+.5, 10*log10(P2(ndx,:)),'--');
                legend({chanlocs.labels}, 'location','northeast')
            else
                numlabels = {'theta','radius','X','Y','Z','sph_theta','sph_phi','sph_radius'};
                tab = struct2table(chanlocs);
                for lab=1:length(numlabels)
                    if ismember(numlabels{lab}, tab.Properties.VariableNames)
                        if iscell(tab.(numlabels{lab}))
                            values = cellfun(@(x)ifthen(isempty(x), nan, double(x)), tab.(numlabels{lab}));
                            tab.(lab) = values;
                        end
                    end
                end
                relX = tab.X ./ sqrt(tab.X.^2+tab.Y.^2+tab.Z.^2);
                relY = tab.Y ./ sqrt(tab.X.^2+tab.Y.^2+tab.Z.^2);
                ant    = relX>=-1E-5;
                post   = ~ant;
                medial = abs(relY)<.41;
                left   = relY>=.41;
                right  = relY<=.41;
                regP = nan(size(P,1),6);
                regP(:,1) = mean(P(:,ant & left),2);
                regP(:,2) = mean(P(:,ant & medial),2);
                regP(:,3) = mean(P(:,ant & right),2);
                regP(:,4) = mean(P(:,post & left),2);
                regP(:,5) = mean(P(:,post & medial),2);
                regP(:,6) = mean(P(:,post & right),2);
                regP2 = nan(size(P2,1),6);
                regP2(:,1) = mean(P2(:,ant & left),2);
                regP2(:,2) = mean(P2(:,ant & medial),2);
                regP2(:,3) = mean(P2(:,ant & right),2);
                regP2(:,4) = mean(P2(:,post & left),2);
                regP2(:,5) = mean(P2(:,post & medial),2);
                regP2(:,6) = mean(P2(:,post & right),2);
                plot(ax, fs(ndx)+.5, 10*log10(regP(ndx,:)));
                hold on
                ax.ColorOrderIndex = 1;
                plot(ax, fs(ndx)+.5, 10*log10(regP2(ndx,:)),'--');
                legend('ant left','ant medial','ant right','post left','post medial','post right')
            end
            set(ax, 'fontsize', fontsize+2)
            xlabel('frequency (Hz)')
            ylabel('Power({\mu}V^2/Hz)')
            set(line([1 1], ylim), 'linestyle', '--');
            set(line([4 4], ylim), 'linestyle', '--');
            set(line([8 8], ylim), 'linestyle', '--');
            set(line([13 13], ylim), 'linestyle', '--');
            set(line([30 30], ylim), 'linestyle', '--');
        end

        % Overwrite/skip button callback (bound to a uicontrol in the dialog)
        function buttonCallback(app, f, choice) %#ok<INUSL>
            global skipOverwrite_choice %#ok<GVMIS>
            skipOverwrite_choice = choice;
            uiresume(f);
            delete(f);
        end

        % Button pushed function: pushbuttonEditCode
        function pushbuttonEditCodeButtonPushed(app, event) %#ok<INUSD>
            EegCallbacks.openBigEditor(app);
        end

        % Button pushed function: pushbuttonExecuteCode
        function pushbuttonExecuteCodeButtonPushed(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            EEG = data.EEG;
            if ~iscell(app.textareaExecuteCode.Value)
                cmd = app.textareaExecuteCode.Value;
            else
                cmd = app.textareaExecuteCode.Value{1};
                for i=2:length(app.textareaExecuteCode.Value)
                    cmd = [cmd sprintf('\n%s', app.textareaExecuteCode.Value{i})]; %#ok<AGROW>
                end
            end
            try
                eval(cmd);
                EEG.history = [EEG.history '\n' cmd];
            catch E
                EegCallbacks.AddToListbox(app, app.listboxStdout, 'An error occured on executing the code');
                disp(E.message)
                if EegCallbacks.isBatchEvent(event)
                    rethrow(E);          % batch: abort file, do NOT save
                else
                    msgbox(sprintf('Error executing code:\n%s', E.message));
                    return               % interactive: leave data.EEG untouched, do not save the failed step
                end
            end
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, 'Execute code');
            data.EEG = EEG;
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, 'ExecuteCommand');
            guidata(hObject, data);
        end

        % Button pushed function: pushbuttonIntClean
        function pushbuttonIntCleanButtonPushed(app, event)
            data = guidata(app.eeg_workflow);
            method = EegParams.get(data.params, 'interpolation clean', 'method');
            if strcmpi(method, 'RANSAC')
                how = sprintf('RANSAC prediction, median of %d subsets of %.0f%% of the channels', ...
                    EegParams.get(data.params, 'interpolation clean', 'ransac_draws'), ...
                    100*EegParams.get(data.params, 'interpolation clean', 'ransac_fraction'));
            else
                how = 'leave-one-out prediction from all other channels';
            end
            EegCallbacks.runPeriods(app, event, 'interpolation clean', ...
                sprintf('Removing channels/periods that differ from their interpolation (%s).', how), ...
                @(EEG, chans, ep, P) EegPeriods.maskInterpolation(EEG, chans, ep, ...
                    P('sdcrit'), P('rcrit'), P('method'), P('ransac_draws'), P('ransac_fraction')), ...
                'InterpolationCleaning');
        end

        % Button pushed function: pbBatch_v2
        function pbBatchButtonPushed(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            global skipOverwrite_choice %#ok<GVMIS>
            data = guidata(hObject);


            % new code for CLAUDE to know what the order of the buttons is,
            % the label for the buttons

            
            str = {
                '--------'
                'Flatline'
                'Bad channels'
                'Lookup'
                'Resample'
                'Rereference'
                'Cleanline'
                'Filter'
                'EOG'
                'Alt EOG'
                'EMG'
                'Alt EMG'
                'ASR'
                'ICA'
                '--------'
                'Execute code'
                'Remove ~EEG'
                'Remove chans'
                '--------'
                'Flat periods'
                'Excessive periods'
                'EMG periods'
                'Interpolation periods'
                '--------'
                'Remove 1st s'
                'Remove last s'
                'Remove task'
                'Remove no-task'
                };
            buttons = struct2table(struct('idx', num2cell( (1:numel(str))' ), 'label',str));
            % buttons now holds the order in idx to be used in the case.
            % str is now to be passed to the new figure-creating script,
            % not the fixed script it used to be before where the dropdown items were manually filled. 
            
            if ~isfield(data, 'DEFAULTDIR')
                data.DEFAULTDIR = EegCallbacks.loadDefaultDir(data.SETTINGSDIR);
            end
            % The window returns at once; the run starts from its Start button
            % (a uifigure cannot dispatch its own buttons while this callback
            % is still on the stack).
            EegBatch.dialog(hObject, data.SETTINGSDIR, data.DEFAULTDIR, str, ...
                @(sel) EegCallbacks.runBatch(app, sel));
        end

        % ------------------------------------------------------------------
        % Run the batch described by sel (see EegBatch.dialog). sel.steps holds
        % the indices into the step list built above, in the chosen order.
        function runBatch(app, sel)
            hObject = app.eeg_workflow;
            global skipOverwrite_choice %#ok<GVMIS>
            data = guidata(hObject);
            data.batchfilenames = sel.files;
            data.batchpathname  = sel.folder;
            data.batchprefix    = sel.prefix;
            data.batchoutputdir = sel.outputdir;
            data.batchsteps     = sel.steps;
            guidata(hObject, data);
            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                'Batch: %d files, steps: %s', numel(sel.files), strjoin(sel.labels, ' > ')));
            if isfield(data,'batchfilenames') &&  ~isempty(data.batchfilenames)
                resultTable = [];
                data.batchfilenamesout = {};
                for f=1:length(data.batchfilenames)
                    [~, Fname, Fext] = fileparts(data.batchfilenames{f});
                    FNOut = [data.batchoutputdir '/' sprintf('%s%s%s', data.batchprefix, Fname, Fext)];
                    data.batchfilenamesout{end+1} = FNOut;
                end
                ex = false;
                for f=1:length(data.batchfilenamesout)
                    if exist(data.batchfilenamesout{f},'file')
                        ex = true;
                        break
                    end
                end
                if ex
                    EegCallbacks.skip_or_overwrite_dialog(app);
                end
                for f=1:length(data.batchfilenames)
                    try
                        FNOut = data.batchfilenamesout{f};
                        if exist(FNOut) && strcmpi(skipOverwrite_choice, 'skip')
                            continue
                        end
                        [Fpath, Fname, Fext] = fileparts(data.batchfilenames{f});
                        EEG = EegCallbacks.loadfile(app, Fpath, [Fname Fext], ...
                            app.listboxStdout, ...
                            EegParams.get(data.params,'open','biosig'), ...
                            EegParams.get(data.params,'open','avgref'));
                        if ~isfield(EEG, 'event') || isempty(EEG.event)
                            EEG.event=struct('type', {}, 'latency', {}, 'urevent', {});
                        end
                        data.EEG = EEG;
                        if ~isfield(data, 'Stack') || isempty(data.Stack)
                            data.Stack = {};
                        end
                        if ~isfield(data, 'StackLabel') || isempty(data.StackLabel)
                            data.StackLabel = {};
                        end
                        data.tabLine = struct();
                        data.tabLine.starttime = size(data.EEG.data(:,:),2) / data.EEG.srate;
                        guidata(hObject, data);
                        EegCallbacks.trackSpectrum(app, data.EEG, 'Open', [Fname Fext]);
                        EegCallbacks.ResetHistory(app, sprintf('Open %s', [Fname Fext]));
                        pause(0.01);
                        % Dropdown item -> action map (see original comments).

                        

                        for cb = 1:numel(sel.steps)
                            EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                                '== %s ==', sel.labels{cb}));
                            switch sel.steps(cb)
                                case 2,  EegCallbacks.pushbuttonFlatline_Callback(app, 'batch');
                                case 3,  EegCallbacks.pushbuttonExcessive_Callback(app, 'batch');
                                case 4,  EegCallbacks.pushbuttonChanlocs_Callback(app, 'batch');
                                case 5,  EegCallbacks.pushbuttonResample_Callback(app, 'batch');
                                case 6,  EegCallbacks.pushbuttonRereference_Callback(app, 'batch');
                                case 7,  EegCallbacks.pushbuttonLineNoise_Callback(app, 'batch');
                                case 8,  EegCallbacks.pushbuttonFilter_Callback(app, 'batch');
                                case 9,  EegCallbacks.pushbuttonInitialICA_Callback(app, 'batch');
                                case 10, EegCallbacks.pushbuttonAltEOG_Callback(app, 'batch');
                                case 11, EegCallbacks.pushbuttonEMG_Callback(app, 'batch');
                                case 12, EegCallbacks.pushbuttonAltEMG_Callback(app, 'batch');
                                case 13, EegCallbacks.pushbuttonASR_Callback(app, 'batch');
                                case 14, EegCallbacks.pushbuttonICA_Callback(app, 'batch');
                                case 16, EegCallbacks.pushbuttonExecuteCodeButtonPushed(app, 'batch');
                                case 17, EegCallbacks.pushbuttonRemoveNoEEG_Callback(app, 'batch');
                                case 18, EegCallbacks.pushbuttonRemoveEOG_Callback(app, 'batch');
                                case 20, EegCallbacks.pushbuttonFlatPeriods_Callback(app, 'batch');
                                case 21, EegCallbacks.pushbuttonExcessivePeriodsButtonPushed(app, 'batch');
                                case 22, EegCallbacks.pushbuttonEMGPeriods_Callback(app, 'batch');
                                case 23, EegCallbacks.pushbuttonIntCleanButtonPushed(app, 'batch');
                                case 25, EegCallbacks.pushbuttonRemoveFirstSec_Callback(app, 'batch');
                                case 26, EegCallbacks.pushbuttonLastSecButtonPushed(app, 'batch');
                                case 27, EegCallbacks.pushbuttonRemoveTask_Callback(app, 'batch');
                                case 28, EegCallbacks.pushbuttonRemoveResting_Callback(app, 'batch');
                                otherwise        % separators: nothing to do
                            end
                            data = guidata(hObject);
                        end
                        data.Stack = {};
                        data.StackLabel = {};
                        guidata(hObject, data);
                        [filepath, filename, fileext] = fileparts(data.batchfilenamesout{f});
                        if isempty(filepath)
                            filepath = data.batchpathname;
                        end
                        pop_saveset(data.EEG, 'filename', [filename, fileext], ...
                                              'filepath', filepath, ...
                                              'savemode', 'onefile');
                        data.tabLine.finaltime = size(data.EEG.data(:,:),2) / data.EEG.srate;
                        try
                            data.tabLine.filename = data.batchfilenames{f};
                            if isempty(resultTable)
                                resultTable = data.tabLine;
                            else
                                resultTable = cat(1, resultTable, data.tabLine);
                            end
                            rowData = struct2cell(data.tabLine);
                            fid = fopen(sprintf('%s/Metadata_Batch.txt', data.batchpathname), 'a');
                            if fid>0
                                first = true;
                                for c=1:length(rowData)
                                    if ~first, fprintf(fid, '\t'); end
                                    switch class(rowData{c})
                                        case {'int8','int16','int32','int64','uint8','uint16','uint32','uint64'}
                                            fprintf(fid, '%d', rowData{c});
                                        case 'double'
                                            fprintf(fid, '%.4g', rowData{c});
                                        case 'char'
                                            fprintf(fid, '%s', rowData{c});
                                    end
                                    first = false;
                                end
                                fprintf(fid, '\n');
                            end
                            fclose(fid);
                        catch
                            warning('collecting metadata for this subject failed.');
                        end
                    catch E
                        fprintf(2, 'Error in file %s:\n%s\n', data.batchfilenames{f}, E.getReport());
                        continue
                    end
                end
            end
        end

        % Button down function: UITrackAxes
        % Opens the spectrum tracker in a larger, non-modal viewer window. The
        % window is kept: closing it only hides it, and a next click shows the
        % same window again with the current contents of the tracker.
        function UITrackAxesButtonDown(app, event) %#ok<INUSD>
            src = app.UITrackAxes;
            fig = findall(groot, 'Type', 'figure', 'Tag', 'EegTrackViewer');
            if isempty(fig)
                fig = figure('Name', 'Spectrum tracker', 'NumberTitle', 'off', ...
                    'Tag', 'EegTrackViewer', 'Position', [200 200 820 560], ...
                    'CloseRequestFcn', @(f,~) set(f, 'Visible', 'off'));
                axes('Parent', fig, 'Tag', 'EegTrackViewerAxes');
            end
            fig = fig(1);
            ax = findall(fig, 'Tag', 'EegTrackViewerAxes');

            cla(ax, 'reset');
            set(ax, 'Tag', 'EegTrackViewerAxes');   % 'reset' clears the tag too
            hold(ax, 'on');
            % copy bottom-most first so the stacking (band behind lines) stays
            objs = flipud(allchild(src));
            for k = 1:numel(objs)
                if isgraphics(objs(k), 'line') || isgraphics(objs(k), 'patch')
                    c = copyobj(objs(k), ax);
                    c.HitTest = 'on';                    % datatips in the viewer
                    if isgraphics(c, 'line')
                        c.LineWidth = 1.5 * c.LineWidth;
                    end
                end
            end
            hold(ax, 'off');
            ax.XLim = src.XLim;
            ax.YLim = src.YLim;
            ax.FontSize = 12;
            ax.Box = 'on';
            grid(ax, 'on');
            title(ax, src.Title.String, 'Interpreter', 'none');
            xlabel(ax, 'Frequency (Hz)');
            ylabel(ax, 'Power (dB)');
            if ~isempty(findobj(ax, 'Tag', 'trackMean'))
                legend(ax, 'show', 'Location', 'northeast', 'Interpreter', 'none');
            end
            fig.Visible = 'on';
            figure(fig);
        end

        % Button pushed function: pushbuttonOpen
        function pushbuttonOpen_Callback(app, event) %#ok<INUSD>
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(b,p) EegParams.get(data.params, b, p);
            savecolour = app.pushbuttonOpen.BackgroundColor; % restore after cancel
            app.pushbuttonOpen.BackgroundColor = EegCallbacks.ColBusy;

            % First row is the default filter. R2026 dedupes any row containing
            % '*.*' into the auto "All Files (*)" entry it appends, so keep the
            % default row free of '*.*': a combined list of the readable types.
            FilterSpec = {'*.bdf;*.cnt;*.edf;*.set;*.vhdr', 'Any readable EEG'
                '*.bdf', 'Biosemi'
                '*.cnt', 'ANT Neuro / Neuroscan'
                '*.edf', 'European data format'
                '*.set', 'EEGLAB'
                '*.vhdr', 'BrainVision'
                };
            if ~isfield(data, 'DEFAULTDIR')
                data.DEFAULTDIR = EegCallbacks.loadDefaultDir(data.SETTINGSDIR);
            end
            [FileName, PathName, FilterIndex] = uigetfile(FilterSpec,'Select an EEG file', ...
                [data.DEFAULTDIR filesep]); %#ok<ASGLU>
            EegCallbacks.bringToFront(app.eeg_workflow);

            % remember the folder, but not when the dialog was cancelled
            if ischar(PathName)
                data.DEFAULTDIR = EegCallbacks.saveDefaultDir(data.SETTINGSDIR, PathName);
                guidata(hObject, data);
            end

            if isnumeric(FileName) && FileName==0
                % cancelled: keep whatever was loaded, including its table line
                EegCallbacks.AddToListbox(app, app.listboxStdout, '*** warning *** no file selected');
                app.pushbuttonOpen.BackgroundColor = savecolour;
                return
            else
                data.EEG = EegCallbacks.loadfile(app, PathName, FileName, app.listboxStdout, ...
                    P('open','biosig'), P('open','avgref'));
                data.EEG = eeg_checkset(data.EEG);

                % mask event types to their low 8 bits (was its own button)
                if P('open','mask')
                    data.EEG = EegCallbacks.maskEvents(app, app.listboxStdout, data.EEG);
                end

                % correct nonstandard EOG channel naming
                if P('open','rename')
                    EegCallbacks.AddToListbox(app, app.listboxStdout, '- ranaming EOGV and EOGH');
                    adjust = FindSetNdx({data.EEG.chanlocs.labels},{'EOGH'});
                    if ~isempty(adjust)
                        data.EEG.chanlocs(adjust).labels = 'HEOG';
                    end
                    adjust = FindSetNdx({data.EEG.chanlocs.labels},{'EOGV'});
                    if ~isempty(adjust)
                        data.EEG.chanlocs(adjust).labels = 'VEOG';
                    end
                end

                % reset the stack for Undo operations
                data.Stack = {};
                data.StackLabel = {};

                % A file is open: Open itself is done (idle) and every other
                % workflow button becomes available (ready).
                EegCallbacks.setMainButtonsColor(app, EegCallbacks.ColReady);
                app.pushbuttonOpen.BackgroundColor = EegCallbacks.ColIdle;
                EegCallbacks.trackSpectrum(app, data.EEG, 'Open', FileName);
                EegCallbacks.ResetHistory(app, sprintf('Open %s', FileName));
                EegCallbacks.bringToFront(app.eeg_workflow);
            end

            % start a meta-data table line.
            data.tabLine = struct();
            data.tabLine.starttime = size(data.EEG.data(:,:),2) / data.EEG.srate;

            guidata(hObject,data)
        end


        % ------------------------------------------------------------------
        % Give a window the focus again. On macOS the native file and folder
        % dialogs belong to the MATLAB desktop, so when they close the focus
        % goes to the Command Window rather than back to the app.
        function bringToFront(fig)
            try
                if ~isempty(fig) && isvalid(fig)
                    drawnow;
                    figure(fig);
                end
            catch
            end
        end

        % ------------------------------------------------------------------
        % Settings folder: the system's application-data location.
        %   Windows  %APPDATA%\Matlab_EegAutoFlow
        %   macOS    ~/Library/Application Support/Matlab_EegAutoFlow
        %   Linux    $XDG_CONFIG_HOME/Matlab_EegAutoFlow (default ~/.config)
        % ------------------------------------------------------------------
        % RANSAC bad-channel test of the Bad chans button, as PREP's
        % findNoisyChannels: a channel is bad when it correlates below
        % ransac_r with its RANSAC prediction (eeg_ransac) in more than
        % ransac_maxbad of the recording. Tested and used as predictors: the
        % located channels not named *eog*, without flat ones. Returns the
        % indices of the bad channels in EEG; removing them is the caller's.
        function hit = ransacBadChannels(app, EEG, P)
            say = @(varargin) EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(varargin{:}));
            hit = [];
            rc = EegPeriods.eegChannels(EEG);
            flat = std(double(EEG.data(rc, :)), [], 2)' < 1e-9;
            if any(flat)
                say('  *** warning *** RANSAC skips %d flat channels (%s): run Flatline first', ...
                    sum(flat), strjoin({EEG.chanlocs(rc(flat)).labels}, ' '));
                rc = rc(~flat);
            end
            say('- RANSAC (PREP): r < %.2f in more than %.0f%% of %g s windows; median of %d subsets of %.0f%% of %d channels', ...
                P('ransac_r'), 100*P('ransac_maxbad'), P('ransac_window'), P('ransac_draws'), ...
                100*P('ransac_fraction'), numel(rc));
            try
                R = eeg_ransac(EEG, 'Channels', rc, 'WindowSeconds', P('ransac_window'), ...
                    'Draws', P('ransac_draws'), 'Fraction', P('ransac_fraction'));
            catch E
                say('  *** warning *** RANSAC not run: %s', E.message);
                return
            end
            % PREP's rule: bad windows x window length > max bad time x recording
            badTime = sum(R < P('ransac_r'), 2)' * P('ransac_window') * EEG.srate;
            isBad   = badTime > P('ransac_maxbad') * EEG.pnts;
            if all(isBad)
                say('  *** error *** every channel failed the RANSAC test; none removed');
                return
            end
            hit = rc(isBad);
            if isempty(hit)
                say('  no channels');
            else
                pct = 100 * badTime(isBad) / EEG.pnts;
                say('  %d channels: %s', numel(hit), strjoin(arrayfun(@(k, q) ...
                    sprintf('%s (%.0f%%)', EEG.chanlocs(k).labels, q), hit, pct, 'uni', 0), ' '));
            end
        end

        % ------------------------------------------------------------------
        % Full path of a file in resources/ or datasets/ next to the repository
        % root (this class lives in code/). Data files, unlike code, are not
        % found through the MATLAB path, so every channel-location file and
        % lookup dataset is resolved here. A name that is not in either folder
        % is returned unchanged, so files shipped with EEGLAB still work.
        function FN = resourceFile(name)
            here = fileparts(which('EegCallbacks'));
            FN = name;
            if isempty(here)
                return
            end
            root = fileparts(here);
            for d = {'resources', 'datasets'}
                cand = fullfile(root, d{1}, name);
                if exist(cand, 'file')
                    FN = cand;
                    return
                end
            end
        end

        function d = settingsDir()
            if ispc
                base = getenv('APPDATA');
                if isempty(base)
                    base = fullfile(getenv('USERPROFILE'), 'AppData', 'Roaming');
                end
            elseif ismac
                base = fullfile(getenv('HOME'), 'Library', 'Application Support');
            elseif isunix
                base = getenv('XDG_CONFIG_HOME');
                if isempty(base)
                    base = fullfile(getenv('HOME'), '.config');
                end
            else
                base = pwd;
            end
            d = fullfile(base, 'Matlab_EegAutoFlow');
        end

        % ------------------------------------------------------------------
        % Copy settings files from the folders earlier versions used
        % (~/Application Support on macOS; a literal '~' folder in the current
        % folder on Windows, where MATLAB does not expand '~') into newDir.
        % Files already present in newDir are left alone, and the old folder
        % is not touched: the previous EegWorkflowApp still reads it.
        function [n, from] = migrateSettings(newDir)
            n = 0;
            from = '';
            old = {};
            if ismac
                old{end+1} = fullfile(getenv('HOME'), 'Application Support', 'Matlab_EegAutoFlow');
            elseif ispc
                old{end+1} = fullfile(pwd, '~', 'AppData', 'Matlab_EegAutoFlow');
            elseif isunix
                old{end+1} = fullfile(getenv('HOME'), '.config', 'Matlab_EegAutoFlow');
            end
            for k = 1:numel(old)
                if strcmp(old{k}, newDir) || ~isfolder(old{k})
                    continue
                end
                files = dir(old{k});
                files = files(~[files.isdir]);
                for i = 1:numel(files)
                    target = fullfile(newDir, files(i).name);
                    if ~isfile(target)
                        if ~isfolder(newDir)
                            mkdir(newDir);
                        end
                        copyfile(fullfile(old{k}, files(i).name), target);
                        n = n + 1;
                        from = old{k};
                    end
                end
            end
        end

        % ------------------------------------------------------------------
        % The folder file dialogs start in, kept in
        % <SETTINGSDIR>/.EegWorkflow_DefaultPath.ini. Falls back to the folder
        % the batch window used to remember, then to the current folder.
        function d = loadDefaultDir(SETTINGSDIR)
            d = '';
            for f = {'.EegWorkflow_DefaultPath.ini', '.figRunBatch_DefaultPath.ini'}
                fid = fopen(fullfile(SETTINGSDIR, f{1}), 'r');
                if fid > 0
                    line = fgetl(fid);
                    fclose(fid);
                    if ischar(line)
                        line = strtrim(line);
                        if ~isempty(line) && isfolder(line)
                            d = line;
                            break
                        end
                    end
                end
            end
            if isempty(d)
                d = pwd;
            end
            d = EegCallbacks.stripSep(d);
        end

        % ------------------------------------------------------------------
        function d = saveDefaultDir(SETTINGSDIR, d)
            d = EegCallbacks.stripSep(d);
            if ~isfolder(SETTINGSDIR)
                mkdir(SETTINGSDIR);
            end
            fid = fopen(fullfile(SETTINGSDIR, '.EegWorkflow_DefaultPath.ini'), 'w');
            if fid > 0
                fprintf(fid, '%s', d);
                fclose(fid);
            end
        end

        % ------------------------------------------------------------------
        function d = stripSep(d)
            d = char(d);
            while numel(d) > 1 && any(d(end) == '/\')
                d(end) = [];
            end
        end

    end % methods (Static)

    % ----------------------------------------------------------------------
    methods (Static, Access = private)

        % ------------------------------------------------------------------
        % Batch/interactive mode is carried by the callback's own `event` arg:
        % the batch loops invoke each callback with the sentinel event 'batch',
        % while a real UI click passes a ButtonPushedData object.  This keeps
        % the mode local to each call — it can never be stale.
        function tf = isBatchEvent(event)
            tf = (ischar(event) || (isstring(event) && isscalar(event))) ...
                 && strcmpi(char(event), 'batch');
        end

        % Abort the current step.  In batch this THROWS so the file loop catches
        % it and skips saving; in interactive it shows a dialog and the caller
        % returns normally.
        function abortStep(app, event, msg)
            step = EegCallbacks.callerStepName();
            % messages like 'Filter: ...' already name the step
            msg = regexprep(char(msg), ['^' regexptranslate('escape', step) '\s*:\s*'], '', 'ignorecase');
            try
                EegCallbacks.AddToListbox(app, app.listboxStdout, ...
                    sprintf(' *** %s: not executed - %s ***', step, msg));
            catch
            end
            EegCallbacks.AddHistoryNote(app, sprintf('%s: not executed (%s)', step, msg));
            if EegCallbacks.isBatchEvent(event)
                error('EegWorkflow:abort', '%s', msg);
            else
                msgbox(msg);
            end
        end

        % ------------------------------------------------------------------
        % pop_clean_rawdata with every criterion taken from the settings of
        % button key. A value of 0 switches that criterion off, so the same
        % code serves a combined pass (the usual ASR step), a bad-channels-only
        % pass or a bad-segments-only pass (as DISCOVER-EEG uses).
        function runCleanRawdata(app, event, key, btnName, undoLabel, powerTag)
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            P = @(p) EegParams.get(data.params, key, p);
            set(app.(btnName),'backgroundcolor',EegCallbacks.ColBusy)
            pause(0.005);
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                set(app.(btnName),'backgroundcolor',EegCallbacks.ColIdle)
                return
            end
            onoff = @(v) ifthen(v > 0, v, 'off');

            flat   = P('flatline');
            chanr  = P('chanminr');
            lineSD = P('linenoise');
            burst  = P('burstcriterion');
            window = P('windowcriterion');
            % clean_channels needs a numeric correlation criterion whenever
            % the line-noise test runs; -1 never flags a channel
            chanArg = onoff(chanr);
            if chanr <= 0 && lineSD > 0
                chanArg = -1;
            end
            say = @(varargin) EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(varargin{:}));
            say('Clean data using clean_rawdata');
            say('- flatline %s | channel r %s | line noise %s | max bad time %.2f', ...
                ifthen(flat > 0, sprintf('%g s', flat), 'off'), ifthen(chanr > 0, sprintf('%.2f', chanr), 'off'), ...
                ifthen(lineSD > 0, sprintf('%g SD', lineSD), 'off'), P('maxbadtime'));
            say('- drift highpass %s | burst %s (%s) | window %s', ...
                ifthen(P('highpass'), '0.25-0.75 Hz', 'off'), ifthen(burst > 0, sprintf('%g SD', burst), 'off'), ...
                ifthen(P('burstdelete'), 'removed', 'corrected'), ifthen(window > 0, sprintf('%.2f', window), 'off'));

            [tmp, cmd] = pop_clean_rawdata(data.EEG, ...
                'FlatlineCriterion', onoff(flat), ...
                'ChannelCriterion', chanArg, ...
                'LineNoiseCriterion', onoff(lineSD), ...
                'Highpass', ifthen(P('highpass'), [0.25 0.75], 'off'), ...
                'BurstCriterion', onoff(burst), ...
                'BurstCriterionRefTolerances', [-Inf 3.0], ...
                'WindowCriterion', onoff(window), ...
                'BurstRejection', ifthen(P('burstdelete'), 'on', 'off'), ...
                'Distance', 'Euclidian', ...
                'ChannelCriterionMaxBadTime', P('maxbadtime'));
            tmp.history = [tmp.history newline cmd];

            % short stretches left between the cut-out bursts
            gapSec = P('burstmergegap');
            if gapSec > 0
                [tmp, nShort, secShort] = EegPeriods.removeShortStretches(tmp, round(gapSec * tmp.srate));
                if nShort > 0
                    say('- removed %d stretches shorter than %g s between boundaries (%.1f s)', nShort, gapSec, secShort);
                end
            end

            data.tabLine.cleanDeletedChans = abs(data.EEG.nbchan-tmp.nbchan);
            if (data.tabLine.cleanDeletedChans>0)
                prevdata = pop_select(data.EEG, 'channel', {tmp.chanlocs.labels});
            else
                prevdata = data.EEG;
            end
            data.tabLine.cleanDeletedTime = (prevdata.pnts - tmp.pnts) ./ data.EEG.srate;
            say('- removed %d channels and %.1f s', data.tabLine.cleanDeletedChans, data.tabLine.cleanDeletedTime);

            data = EegCallbacks.pushUndo(app, event, data, data.EEG, EegCallbacks.callerStepName());
            data.EEG = tmp;
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, powerTag);
            guidata(hObject, data);
            set(app.(btnName),'backgroundcolor',EegCallbacks.ColIdle)
            pause(0.005);
        end

        % ------------------------------------------------------------------
        % Average reference, excluding HEOG/VEOG. With interpRemoved the
        % channels removed earlier (EEG.chaninfo.removedchans, those with a
        % location) are interpolated first so the average is taken over the
        % full montage, then removed again (as pop_reref 'interpchan' [] and
        % DISCOVER-EEG). Done here rather than through pop_reref so that the
        % HEOG/VEOG exclusion indexes the channels after interpolation.
        function [EEG, cmd] = averageReference(app, EEG, interpRemoved)
            added = {};
            if interpRemoved
                rc = [];
                if isfield(EEG.chaninfo, 'removedchans') && ~isempty(EEG.chaninfo.removedchans) ...
                        && isfield(EEG.chaninfo.removedchans, 'theta')
                    rc = EEG.chaninfo.removedchans;
                    rc = rc(~cellfun(@isempty, {rc.theta}));
                    rc = rc(~ismember(upper({rc.labels}), upper({EEG.chanlocs.labels})));
                end
                if isempty(rc)
                    EegCallbacks.AddToListbox(app, app.listboxStdout, '- no removed channels with a location to include');
                else
                    EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                        '- including %d removed channels in the average (interpolated): %s', numel(rc), strjoin({rc.labels}, ' ')));
                    EEG = pop_interp(EEG, rc, 'spherical');
                    added = {rc.labels};
                end
            end
            [EEG, cmd] = pop_reref(EEG, [], 'exclude', find(ismember(upper({EEG.chanlocs.labels}), {'HEOG','VEOG'})));
            if ~isempty(added)
                EEG = pop_select(EEG, 'nochannel', find(ismember({EEG.chanlocs.labels}, added)));
                cmd = sprintf('%s %% average over the full montage: %d removed channels interpolated for the average, then removed again', cmd, numel(added));
            end
        end

        % ------------------------------------------------------------------
        % FIR order for a Hamming-windowed sinc with transition bandwidth df
        % (Hz) at sampling rate fs: 3.3*fs/df, rounded up to an even number.
        function n = firOrderHamming(fs, df)
            n = ceil(3.3 * fs / df);
            if mod(n, 2) ~= 0
                n = n + 1;
            end
        end

        % ------------------------------------------------------------------
        % Number of ICA runs from the 'ICA number' dropdown (1 when absent).
        function n = icaRuns(app)
            n = 1;
            try
                if isprop(app,'ICANumberDropDown') && isvalid(app.ICANumberDropDown)
                    n = max(1, round(str2double(app.ICANumberDropDown.Value)));
                end
            catch
            end
            if isnan(n), n = 1; end
        end

        % ------------------------------------------------------------------
        % One ICA decomposition. Run 1 is exactly the single-run behaviour.
        % runica and binica start from a random point on every call; PICARD as
        % called through pop_runica starts from the identity matrix and would
        % give the same answer each time, so later PICARD runs get a random
        % orthogonal starting matrix (picard 'w_init').
        function [EEG, cmd] = decomposeIca(app, EEG, icatype, ncomps, goodchans, run)
            extra = {};
            if run > 1 && strcmpi(icatype, 'picard')
                [Q, ~] = qr(randn(ncomps));
                extra = {'w_init', Q};
            end
            try
                if strcmpi(icatype,'jader')
                    warning('using the JADER ICA algorithm. Ncomps/ PCA not used')
                    [EEG, cmd] = pop_runica(EEG, 'icatype', 'jader', 'chanind', goodchans);
                else
                    [EEG, cmd] = pop_runica(EEG, 'icatype', icatype, 'pca', ncomps, 'chanind', goodchans, extra{:});
                end
            catch E
                EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf( ...
                    ' *** warning *** %s failed (%s); using extended runica', icatype, E.message));
                [EEG, cmd] = pop_runica(EEG, 'icatype', 'runica', 'extended', 1, 'pca', ncomps, 'chanind', goodchans);
            end
        end

        % ------------------------------------------------------------------
        % ICLabel, then flag components by per-class rules (see icRules): a
        % component is flagged when ANY switched-on class has a probability of
        % at least its threshold. Brain is never a reason to remove.
        function [EEG, reject] = flagIcs(EEG, rules)
            EEG = pop_iclabel(EEG, 'default');
            reject = EegCallbacks.icReject(EEG.etc.ic_classification.ICLabel.classes, ...
                EEG.etc.ic_classification.ICLabel.classifications, rules);
        end

        % ------------------------------------------------------------------
        % Per-class removal rules from the ICA settings: rows {class, threshold}
        % for every class whose use_ flag is on. DISCOVER-EEG: only Muscle and
        % Eye, both 0.8.
        function rules = icRules(params)
            P = @(p) EegParams.get(params, 'ICA', p);
            map = {'Muscle',        'muscle'
                   'Eye',           'eye'
                   'Heart',         'heart'
                   'Line Noise',    'linenoise'
                   'Channel Noise', 'channoise'
                   'Other',         'other'};
            rules = cell(0, 2);
            for k = 1:size(map, 1)
                if P(['use_' map{k,2}])
                    rules(end+1, :) = {map{k,1}, P(map{k,2})}; %#ok<AGROW>
                end
            end
        end

        % ------------------------------------------------------------------
        % reject(c) true when for some rule the probability of that class in
        % component c is >= its threshold. classes: ICLabel class names;
        % probs: components x classes.
        function reject = icReject(classes, probs, rules)
            reject = false(size(probs, 1), 1);
            for r = 1:size(rules, 1)
                col = find(strcmpi(classes, rules{r,1}), 1);
                if ~isempty(col)
                    reject = reject | probs(:, col) >= rules{r,2};
                end
            end
        end

        % ------------------------------------------------------------------
        % Component removal by subtraction, used by ICA, EOG and Alt EOG:
        % the back-projection of the flagged components is subtracted from the
        % channels the decomposition was computed on (icachansind). Unlike
        % pop_subcomp this keeps variance outside the ICA subspace (PCA-reduced
        % ICA) and leaves the stored decomposition unchanged. The activations
        % are recomputed from the weights, so a stale icaact cannot be used.
        function EEG = subtractIcs(EEG, reject)
            reject = logical(reject(:));
            if ~any(reject)
                return
            end
            act = (EEG.icaweights * EEG.icasphere) * double(EEG.data(EEG.icachansind, :));
            EEG.data(EEG.icachansind, :) = EEG.data(EEG.icachansind, :) - EEG.icawinv(:, reject) * act(reject, :);
            EEG.icaact = [];
        end

        % ------------------------------------------------------------------
        % Samples clean_rawdata would keep (true) in a bad-segments-only pass,
        % with the settings DISCOVER-EEG uses for its step 6. Only used to
        % compare ICA runs; the data itself are not changed here.
        function mask = badSegmentMask(EEG, burst, window)
            EEG.icaweights = []; EEG.icasphere = []; EEG.icawinv = []; EEG.icaact = []; EEG.icachansind = [];
            out = pop_clean_rawdata(EEG, 'FlatlineCriterion', 'off', 'ChannelCriterion', 'off', ...
                'LineNoiseCriterion', 'off', 'Highpass', 'off', ...
                'BurstCriterion', burst, 'WindowCriterion', ifthen(window > 0, window, 'off'), ...
                'BurstRejection', 'on', 'Distance', 'Euclidian', 'WindowCriterionTolerances', [-Inf 7]);
            if isfield(out, 'etc') && isfield(out.etc, 'clean_sample_mask')
                mask = logical(out.etc.clean_sample_mask(:)');
            else
                mask = true(1, EEG.pnts);
            end
        end

        % ------------------------------------------------------------------
        % Pick the ICA run whose bad-segment mask is closest to the average
        % mask over all runs.
        %
        % Adapted from preprocessing_select_ICA_rep.m in DISCOVER-EEG:
        %   Gil Avila C, Bott FS, Tiemann L, Hohn VD, May ES, Nickel MM,
        %   Zebhauser PT, Gross J, Ploner M. DISCOVER-EEG: an open, fully
        %   automated EEG pipeline for biomarker discovery in clinical
        %   neuroscience. Sci Data 10, 613 (2023). doi:10.1038/s41597-023-02525-0
        %   https://github.com/crisglav/discover-eeg (CC BY 4.0)
        % Changes: takes the runs' masks as a matrix (runs x samples) instead of
        % a cell array of EEG structs, and also returns the distances. The
        % selection rule itself is unchanged.
        function [selRun, distToAvg] = selectIcaRun(masks)
            csmAverage = mean(masks, 1);
            distToAvg = sum(abs(masks - csmAverage), 2);
            [~, selRun] = min(distToAvg);
        end

        % ------------------------------------------------------------------
        % Shared shell of the four 'periods' steps (see EegPeriods). key is
        % the button key in EegWorkflow_parameters.xlsx, maskFcn the method's
        % criterion: @(EEG, chans, ep, P) -> channel x epoch logical, where
        % P(parname) reads a setting of that button.
        function runPeriods(app, event, key, title, maskFcn, powerTag)
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            say = @(varargin) EegCallbacks.AddToListbox(app, app.listboxStdout, sprintf(varargin{:}));
            if ~isfield(data,'EEG') || isempty(data.EEG.data)
                EegCallbacks.abortStep(app, event, 'No data available');
                return
            end
            if data.EEG.trials > 1
                EegCallbacks.abortStep(app, event, 'Periods can only be removed from continuous data.');
                return
            end
            P = @(p) EegParams.get(data.params, key, p);

            say('%s', title);
            EEG   = data.EEG;
            chans = EegPeriods.eegChannels(EEG);
            ep    = EegPeriods.epochs(EEG, P('epochlen'), P('overlap'));
            say('- %d channels tested (located, not *eog*), %d epochs of %.2f s, %.0f%% overlap', ...
                numel(chans), ep.n, ep.len/EEG.srate, 100*P('overlap'));

            M = maskFcn(EEG, chans, ep, P);
            D = EegPeriods.decide(M, ep, EEG.pnts, EEG.srate, P('maxbadtime'), P('mergegap'));
            D = EegPeriods.absorbShortIslands(D, EEG, round(P('mergegap')*EEG.srate));

            first = find(D.badChans & ~D.extraChans);
            extra = find(D.extraChans);
            if ~isempty(first)
                say('- removing %d channels bad in more than %.0f%% of the epochs:', numel(first), 100*D.maxbadtime);
                for b = first(:)'
                    say('   %s (%.0f%%)', EEG.chanlocs(chans(b)).labels, 100*D.fracBad(b));
                end
            end
            if ~isempty(extra)
                say('- removing %d more channels to bring the removed time towards %.0f%%:', numel(extra), 100*D.maxbadtime);
                for b = extra(:)'
                    say('   %s (%.0f%%)', EEG.chanlocs(chans(b)).labels, 100*D.fracBad(b));
                end
            end
            if ~any(D.badChans)
                say('- no channels removed');
            end
            if 1 - D.remain > D.maxbadtime
                say(' *** warning *** still %.0f%% of the data to be removed; removing more channels did not help further', ...
                    100*(1 - D.remain));
            end
            if isempty(D.regions) && ~any(D.badChans)
                say('- nothing to remove');
                return
            end
            if D.remain < 0.05 || all(D.badChans)
                say(' *** error *** this would remove (nearly) all data; nothing removed. Loosen the criterion.');
                EegCallbacks.abortStep(app, event, sprintf( ...
                    'only %.0f%% of the data would remain; nothing removed', 100*D.remain));
                return
            end
            if D.islands > 0
                say('- including %d stretches shorter than %g s left between boundaries/removed periods', ...
                    D.islands, P('mergegap'));
            end
            if isempty(D.regions)
                say('- no periods removed');
            else
                say('- removing %d periods, %.1f s (%.0f%% of the data remains)', ...
                    size(D.regions,1), (1-D.remain)*EEG.pnts/EEG.srate, 100*D.remain);
            end

            EEG = EegPeriods.apply(EEG, chans, D);
            data = EegCallbacks.pushUndo(app, event, data, data.EEG, EegCallbacks.callerStepName());
            data.EEG = EEG;
            guidata(hObject, data);
            data.EEG = EegCallbacks.recordPower(app, data.EEG, powerTag);
            guidata(hObject, data);
        end

        % ------------------------------------------------------------------
        % Label of the button step that is running, as shown on the button, for
        % the undo stack, the history and 'not executed' notes. Taken from the
        % nearest button callback on the call stack, so helpers such as
        % runPeriods get the label of the button that called them.
        function name = callerStepName()
            name = 'Step';
            try
                st = dbstack(1);
                for k = 1:numel(st)
                    fn = regexprep(st(k).name, '^.*\.', '');
                    label = EegCallbacks.buttonLabel(fn);
                    if ~isempty(label)
                        name = label;
                        return
                    end
                end
            catch
            end
        end

        % Callback function name -> button text. Empty for non-callbacks.
        % Unknown button callbacks fall back to a name derived from the
        % function ('pushbuttonFoo_Callback' -> 'Foo').
        function label = buttonLabel(fn)
            persistent map
            if isempty(map)
                map = containers.Map( ...
                    {'pushbuttonOpen_Callback', 'pushbuttonFlatline_Callback', 'pushbuttonExcessive_Callback', ...
                     'pushbuttonChanlocs_Callback', 'pushbuttonResample_Callback', 'pushbuttonRereference_Callback', ...
                     'pushbuttonLineNoise_Callback', 'pushbuttonFilter_Callback', 'pushbuttonInitialICA_Callback', ...
                     'pushbuttonAltEOG_Callback', 'pushbuttonEMG_Callback', 'pushbuttonAltEMG_Callback', ...
                     'pushbuttonASR_Callback', 'pushbuttonICA_Callback', 'pushbuttonExecuteCodeButtonPushed', ...
                     'pushbuttonRemoveNoEEG_Callback', 'pushbuttonRemoveEOG_Callback', 'pushbuttonFlatPeriods_Callback', ...
                     'pushbuttonExcessivePeriodsButtonPushed', 'pushbuttonEMGPeriods_Callback', 'pushbuttonIntCleanButtonPushed', ...
                     'pushbuttonRemoveTask_Callback', 'pushbuttonRemoveResting_Callback', 'pushbuttonRemoveFirstSec_Callback', ...
                     'pushbuttonLastSecButtonPushed', 'pushbuttonLast4minButtonPushed', 'pushbuttonImputeAll_Callback', ...
                     'pushbuttonReview_Callback', 'pbUndo_Callback', 'pbBatchButtonPushed'}, ...
                    {'Open', 'Flatline', 'Bad chans', ...
                     'Lookup', 'Resample', 'Rereference', ...
                     'Line noise', 'Filter', 'EOG', ...
                     'Alt EOG', 'EMG', 'Alt EMG', ...
                     'ASR', 'ICA', 'Execute code', ...
                     'Rm ~EEG ch', 'Remove chans', 'Flat periods', ...
                     'Excessive periods', 'EMG periods', 'Interpolation Clean', ...
                     'Remove task', 'Remove no-task', 'Rm first 1s', ...
                     'Rm last 1s', 'Keep last 4 min', 'Impute All', ...
                     'Review', 'Undo', 'Batch'});
            end
            label = '';
            if isKey(map, fn)
                label = map(fn);
            elseif ~isempty(regexp(fn, '^(pushbutton|pb)\w*(_Callback|ButtonPushed)$', 'once'))
                label = regexprep(regexprep(fn, '(_Callback|ButtonPushed)$', ''), '^(pushbutton|pb)', '');
            end
        end

        % Push the previous EEG onto the undo stack — but ONLY in interactive
        % mode.  In batch the undo history is pure wasted memory (discarded when
        % the next file opens), so it is skipped entirely.
        function data = pushUndo(app, event, data, prevEEG, label)
            EegCallbacks.AddToHistory(app, label);
            if EegCallbacks.isBatchEvent(event), return; end
            n = numel(data.Stack) + 1;
            data.Stack{n}      = prevEEG;
            data.StackLabel{n} = label;
        end

        % ------------------------------------------------------------------
        % History pane: one line per step that changed the data, in the same
        % order as the undo stack (pushUndo adds the line, Undo removes it).
        function AddToHistory(app, str)
            EegCallbacks.addHistoryLine(app, str, 'step');
        end

        % A red note in the history (a step that was not executed). Undo skips
        % notes, so the history stays in step with the undo stack.
        function AddHistoryNote(app, str)
            EegCallbacks.addHistoryLine(app, str, 'note');
        end

        function addHistoryLine(app, str, kind)
            try
                if ~isprop(app, 'listboxHistory') || ~isvalid(app.listboxHistory)
                    return
                end
                lb = app.listboxHistory;
                [items, kinds] = EegCallbacks.historyState(lb);
                lb.Value = {};
                lb.Items = [items, {char(str)}];
                lb.UserData = [kinds, {kind}];
                EegCallbacks.colourHistory(lb);
                scroll(lb, 'bottom');
            catch
            end
        end

        % Remove the last executed step (Undo); notes after it stay.
        function RemoveLastHistory(app)
            try
                if ~isprop(app, 'listboxHistory') || ~isvalid(app.listboxHistory)
                    return
                end
                lb = app.listboxHistory;
                [items, kinds] = EegCallbacks.historyState(lb);
                k = find(strcmp(kinds, 'step'), 1, 'last');
                if isempty(k)
                    return
                end
                items(k) = [];
                kinds(k) = [];
                lb.Value = {};
                lb.Items = items;
                lb.UserData = kinds;
                EegCallbacks.colourHistory(lb);
            catch
            end
        end

        function [items, kinds] = historyState(lb)
            items = cellstr(lb.Items);
            items = items(:)';
            kinds = lb.UserData;
            if ~iscell(kinds) || numel(kinds) ~= numel(items)
                kinds = repmat({'step'}, 1, numel(items));
            end
        end

        function colourHistory(lb)
            try
                removeStyle(lb);
                notes = find(strcmp(lb.UserData, 'note'));
                if ~isempty(notes)
                    addStyle(lb, uistyle('FontColor', [.8 0 0]), 'item', notes);
                end
            catch
            end
        end

        % Start a new history (a file was opened).
        function ResetHistory(app, str)
            try
                if isprop(app, 'listboxHistory') && isvalid(app.listboxHistory)
                    app.listboxHistory.Value = {};
                    app.listboxHistory.Items = {char(str)};
                    app.listboxHistory.UserData = {'open'};
                    removeStyle(app.listboxHistory);
                end
            catch
            end
        end

        % ------------------------------------------------------------------
        function AddToListbox(app, listboxObject, str)
            % App Designer listboxes hold their lines in Items (a cell array of
            % char), not in the GUIDE-era String property.
            old = listboxObject.Items;
            if ischar(old),   old = {old}; end
            if ischar(str),   new = {str};
            elseif isstring(str), new = cellstr(str);
            elseif iscell(str),   new = str;
            else, error('Input must be a char, string, or cell array of char.');
            end
            combined = [old(:); new(:)]';
            if numel(combined) > 100
                combined = combined(end-99:end);
            end
            listboxObject.Value = {};          % Value must stay inside Items
            listboxObject.Items = combined;
            EegCallbacks.colourLines(listboxObject);
            scroll(listboxObject, 'bottom');
        end

        % ------------------------------------------------------------------
        % Colour the output lines by content: '*** error' and 'not executed'
        % red and bold, '*** warning' orange. Styles belong to item positions,
        % so they are rebuilt whenever the lines change.
        function colourLines(lb)
            try
                removeStyle(lb);
                items = cellstr(lb.Items);
                err  = contains(items, '*** error', 'IgnoreCase', true) | ...
                       contains(items, 'not executed', 'IgnoreCase', true);
                warn = contains(items, '*** warning', 'IgnoreCase', true) & ~err;
                if any(err)
                    addStyle(lb, uistyle('FontColor', [.8 0 0], 'FontWeight', 'bold'), 'item', find(err));
                end
                if any(warn)
                    addStyle(lb, uistyle('FontColor', [.85 .45 0]), 'item', find(warn));
                end
            catch
            end
        end

        % ------------------------------------------------------------------
        % Spectrum tracker in app.UITrackAxes: after every step the mean power
        % spectrum over all channels is added, so the effect of each step on
        % the spectrum can be followed. reset (with a title) clears the axes;
        % used when a file is opened.
        %
        %   pfft with hanning(2*srate) and 50% overlap -> 0.5 Hz resolution.
        %   DC is dropped; 1 Hz bins average two points: [1 1.5], [2 2.5], ...
        %   Power is converted to dB per channel; the line is the mean over
        %   channels, the grey band +/- 1 SD over channels (latest step only).
        %   Earlier steps stay as thin lines. Never lets plotting break a step.
        function trackSpectrum(app, EEG, label, resetTitle)
            try
                if ~isprop(app, 'UITrackAxes') || ~isvalid(app.UITrackAxes) || isempty(EEG) || isempty(EEG.data)
                    return
                end
                ax = app.UITrackAxes;
                if nargin >= 4
                    cla(ax);
                    title(ax, resetTitle, 'Interpreter', 'none', 'FontSize', 9);
                end
                [x, m, sd] = EegCallbacks.binnedSpectrum(EEG);
                if isempty(x)
                    return
                end
                % earlier steps: thin lines, no band
                delete(findall(ax, 'Tag', 'trackSD'));      % findall: the band is hidden from legends
                set(findobj(ax, 'Tag', 'trackMean'), 'LineWidth', 0.5);
                hold(ax, 'on');
                % HitTest off: a click on a line or band reaches the axes'
                % ButtonDownFcn, which opens the larger viewer
                fill(ax, [x fliplr(x)], [m+sd fliplr(m-sd)], [.5 .5 .5], ...
                    'FaceAlpha', .25, 'EdgeColor', 'none', 'Tag', 'trackSD', ...
                    'HandleVisibility', 'off', 'HitTest', 'off');
                % remember the undo depth this line belongs to, so Undo can
                % remove exactly the lines of the steps it takes back
                depth = 0;
                try
                    gd = guidata(app.eeg_workflow);
                    if isfield(gd, 'Stack'), depth = numel(gd.Stack); end
                catch
                end
                plot(ax, x, m, 'LineWidth', 1.5, 'Tag', 'trackMean', 'DisplayName', label, ...
                    'HitTest', 'off', 'UserData', struct('depth', depth));
                disableDefaultInteractivity(ax);   % no pan/zoom/datatips on the thumbnail
                hold(ax, 'off');
                ax.XTickMode = 'auto';
                ax.YTickMode = 'auto';
                % after a Filter step: only the passband (stored by the filter)
                xl = [x(1) x(end)];
                if isfield(EEG, 'etc') && isstruct(EEG.etc) && isfield(EEG.etc, 'filterBand')
                    fb = EEG.etc.filterBand;
                    xl = [max(xl(1), fb(1)), min(xl(2), fb(2))];
                    if xl(2) <= xl(1)
                        xl = [x(1) x(end)];
                    end
                end
                ax.XLim = xl;
                % y-range from what is visible, so out-of-band values
                % (e.g. far below the lowpass) do not squash the plot
                lines = findobj(ax, 'Tag', 'trackMean');
                ys = [];
                for L = lines(:)'
                    in = L.XData >= xl(1) & L.XData <= xl(2);
                    ys = [ys, L.YData(in)]; %#ok<AGROW>
                end
                sdIn = x >= xl(1) & x <= xl(2);
                ys = [ys, m(sdIn) + sd(sdIn), m(sdIn) - sd(sdIn)];
                ys = ys(isfinite(ys));
                if numel(ys) > 1 && max(ys) > min(ys)
                    margin = 0.05 * (max(ys) - min(ys));
                    ax.YLim = [min(ys) - margin, max(ys) + margin];
                end
                ax.FontSize = 8;
                xlabel(ax, 'Hz');
                ylabel(ax, 'dB');
                % the SD band behind the lines
                band = findall(ax, 'Tag', 'trackSD');
                if ~isempty(band)
                    uistack(band, 'bottom');
                end
                drawnow limitrate
            catch E
                fprintf(2, 'Spectrum tracker: %s\n', E.message);
            end
        end

        % ------------------------------------------------------------------
        % After Undo: remove the tracker lines of the steps taken back (those
        % above the restored undo depth) and redraw the restored data as the
        % newest line, with its SD band, under the label it had.
        function untrackSpectrum(app, EEG, depth)
            try
                if ~isprop(app, 'UITrackAxes') || ~isvalid(app.UITrackAxes)
                    return
                end
                ax = app.UITrackAxes;
                lines = findobj(ax, 'Tag', 'trackMean');          % newest first
                d = arrayfun(@(L) EegCallbacks.lineDepth(L), lines);
                delete(lines(d > depth));
                lines = lines(d <= depth);
                if isempty(lines)
                    return
                end
                label = lines(1).DisplayName;
                delete(lines(1));                                % replotted with its band
                EegCallbacks.trackSpectrum(app, EEG, label);
            catch E
                fprintf(2, 'Spectrum tracker (undo): %s\n', E.message);
            end
        end

        function d = lineDepth(L)
            d = 0;
            if isstruct(L.UserData) && isfield(L.UserData, 'depth')
                d = L.UserData.depth;
            end
        end

        % ------------------------------------------------------------------
        % 1 Hz binned spectrum in dB: x = bin start (Hz), m = mean over
        % channels, sd = SD over channels.
        function [x, m, sd] = binnedSpectrum(EEG)
            win = hanning(2*EEG.srate);
            [P, fs] = pfft(double(EEG.data(:,:))', EEG.srate, win, .5);   % freq x chan
            first = find(fs >= 1, 1);                % skip DC (and 0.5 Hz)
            K = floor((numel(fs) - first + 1) / 2);  % complete 1 Hz bins
            if K < 1
                x = []; m = []; sd = [];
                return
            end
            r1 = first + 2*(0:K-1);
            Pb = (P(r1, :) + P(r1 + 1, :)) / 2;      % [k, k+0.5]
            Pb(Pb <= 0) = NaN;                       % flat channels: no log of 0
            D  = 10*log10(Pb);
            x  = fs(r1);
            x  = x(:)';
            m  = mean(D, 2, 'omitnan')';
            sd = std(D, 0, 2, 'omitnan')';
        end

        % ------------------------------------------------------------------
        % Pure function: takes EEG, returns EEG with recorded power in
        % EEG.etc.recordPower. Only tracks channels with location info.
        function EEG = recordPower(app, EEG, tag)
            EegCallbacks.trackSpectrum(app, EEG, tag);
            hasLoc = false(1, EEG.nbchan);
            for ch = 1:EEG.nbchan
                try
                    hasLoc(ch) = ~isempty(EEG.chanlocs(ch).sph_theta) && ...
                                 ~isempty(EEG.chanlocs(ch).sph_phi);
                catch
                    hasLoc(ch) = false;
                end
            end
            if ~any(hasLoc), return; end

            if ~isfield(EEG.etc, 'recordPowerChannels')
                EEG.etc.recordPowerChannels = {EEG.chanlocs(hasLoc).labels};
            end
            refChannels = EEG.etc.recordPowerChannels;

            locIdx        = find(hasLoc);
            currentLabels = {EEG.chanlocs(locIdx).labels};
            bands = [1 4; 4 7.5; 7.5 13; 13 19; 19 30];
            [P, fs] = pfft(EEG.data(locIdx, :)', EEG.srate, hanning(EEG.srate * 2));

            S = struct;
            S.Tag      = tag;
            S.bands    = bands;
            S.channels = refChannels;
            S.bandpower = nan(length(bands), length(refChannels));
            for f = 1:size(bands, 1)
                ndx     = fs >= bands(f,1) & fs < bands(f,2);
                bandpow = 10*log10(sum(P(ndx, :)) * (fs(2) - fs(1)));
                for c = 1:length(currentLabels)
                    refIdx = strcmp(currentLabels{c}, refChannels);
                    if any(refIdx)
                        S.bandpower(f, refIdx) = bandpow(c);
                    end
                end
            end

            if ~isfield(EEG.etc, 'recordPower')
                EEG.etc.recordPower = {};
            end
            EEG.etc.recordPower{end+1} = S;
        end

        % ------------------------------------------------------------------
        function setFontSize(app, hObject, fs)
            try
                chlist = hObject.Children;
            catch
                return
            end
            for ch=1:length(chlist)
                if isprop(chlist(ch), 'fontsize')
                    try, set(chlist(ch), 'fontsize', fs); catch, beep; end
                end
                if isprop(chlist(1), 'children')
                    if ~isempty(get(chlist(ch), 'children'))
                        EegCallbacks.setFontSize(app, chlist(ch), fs);
                    end
                end
            end
        end

        % ------------------------------------------------------------------
        % ICA algorithm, taken from the ICA type dropdown in the main window.
        % pop_runica wants the lowercase name.
        function icatype = icaType(app)
            icatype = 'picard';
            try
                if isprop(app,'ICAtypeDropDown') && isvalid(app.ICAtypeDropDown)
                    icatype = lower(app.ICAtypeDropDown.Value);
                end
            catch
            end
        end

        % ------------------------------------------------------------------
        % GetUIControlData / SetUIControlData: taken back from the previous
        % EegWorkflowApp. They save and restore the value controls placed
        % directly in the main window (dropdowns, checkboxes, edit fields,
        % sliders, spinners, text areas) as key/value lines in
        % <SETTINGSDIR>/<figure name>.ini. Everything that is a parameter lives in
        % EegParams instead. One change: widgets are keyed by Tag, or by their
        % property name in app when the Tag is empty (App Designer leaves it
        % empty), so that untagged controls such as ICAtypeDropDown and
        % ICANumberDropDown do not all share the key ''.
        function key = widgetKey(app, h)
            key = get(h, 'tag');
            if ~isempty(key)
                return
            end
            props = properties(app);
            for k = 1:numel(props)
                try
                    v = app.(props{k});
                    if isscalar(v) && isgraphics(v) && v == h
                        key = props{k};
                        return
                    end
                catch
                end
            end
        end

        % ------------------------------------------------------------------
        function strlist = GetUIControlData(app, hObject)
            strlist = struct;
            ch = get(hObject,'ch');
            count = 0;
            for c=1:length(ch)
                skip = false;
                tag = EegCallbacks.widgetKey(app, ch(c));
                switch get(ch(c), 'type')
                    case 'uiedit'
                        val = sprintf('%s',get(ch(c),'string'));
                    case {'uicheckbox','uislider','uispinner'}
                        val = sprintf('%.4f',get(ch(c),'value'));
                    case {'uidropdown'}
                        if isprop(ch(c),'valueindex')
                            val = get(ch(c), 'valueindex');
                        else
                            str = get(ch(c), 'value');
                            items = get(ch(c), 'items');
                            val = find(strcmpi(items, str));
                            if isempty(val)
                                val = 0;
                            elseif length(val)>1
                                val = val(1);
                            end
                        end
                        val = sprintf('%.4f', val);
                    case {'uitextarea'}
                        val = get(ch(c), 'value');
                        if iscell(val)
                            tmp = val{1};
                            for i=2:length(val)
                                tmp = [tmp '\n' val{i}]; %#ok<AGROW>
                            end
                            val = tmp;
                        end
                    otherwise
                        skip=true;
                end
                if ~skip
                    count = count + 1;
                    strlist(count).key = tag;
                    strlist(count).val = val;
                end
            end
            strlist = struct2table(strlist, 'AsArray', true);
        end

        % ------------------------------------------------------------------
        function SetUIControlData(app, hObject, strlist)
            ch = get(hObject,'ch');
            T = table;
            T.type = get(ch,'type');
            T.tag = arrayfun(@(h) EegCallbacks.widgetKey(app, h), ch, 'UniformOutput', false);
            for c=1:size(T,1)
                try
                    tag = T.tag{c};
                    fprintf('%s  %s\n',T.type{c}, T.tag{c})
                    switch T.type{c}
                        case 'uiedit'
                            ndx = find(strcmpi(strlist.key, tag));
                            if length(ndx)==1 %#ok<ISCL>
                                if isnumeric(strlist.val)
                                    set(ch(c),'string', sprintf('%.4f', strlist.val(ndx)));
                                else
                                    set(ch(c),'string', sprintf('%s', strlist.val{ndx}));
                                end
                            end
                            pause(0.005)
                        case {'uicheckbox'}
                            ndx = find(strcmpi(strlist.key, tag));
                            if length(ndx)==1 %#ok<ISCL>
                                if isnumeric(strlist.val)
                                    set(ch(c),'value', strlist.val(ndx));
                                else
                                    set(ch(c),'value', str2num(strlist.val{ndx})); %#ok<ST2NM>
                                end
                            end
                            pause(0.005);
                        case {'uidropdown'}
                            ndx = find(strcmpi(strlist.key, tag));
                            if length(ndx)==1 %#ok<ISCL>
                                fld = ifthen(isprop(ch(c),'value'), 'value', 'valueindex');
                                if isnumeric(strlist.val)
                                    set(ch(c), fld, strlist.val(ndx));
                                elseif strcmpi(fld, 'valueindex')
                                    val = strlist.val{ndx};
                                    set(ch(c), fld, str2num(val)); %#ok<ST2NM>
                                else
                                    val = strlist.val{ndx};
                                    if ~isempty(str2num(val)) %#ok<ST2NM>
                                        items = get(ch(c),'items');
                                        set(ch(c), fld, items{str2num(val)}); %#ok<ST2NM>
                                    else
                                        set(ch(c), fld, val);
                                    end
                                end
                            end
                            pause(0.005);
                        case {'uislider','uispinner'}
                            ndx = find(strcmpi(strlist.key, tag));
                            if length(ndx)==1 %#ok<ISCL>
                                if isnumeric(strlist.val)
                                    set(ch(c),'value', strlist.val(ndx));
                                else
                                    set(ch(c),'value', str2num(strlist.val{ndx})); %#ok<ST2NM>
                                end
                            end
                            pause(0.005);
                        case 'uitextarea'
                            ndx = find(strcmpi(strlist.key, tag));
                            if length(ndx)==1 %#ok<ISCL>
                                tmp = strsplit(sprintf(strlist.val{ndx}),'\n');
                                set(ch(c),'value', tmp);
                            end
                            pause(0.005)
                    end
                catch
                    pause(.05)
                end
            end
        end

        % ------------------------------------------------------------------
        % The EEG properties panel is not part of the current layout. Kept so
        % that re-adding a listboxEegProperties component makes it work again.
        function listboxEegProperties_Update(app, hObject)
            if ~isprop(app, 'listboxEegProperties')
                return
            end
            data = guidata(hObject);
            fields = fieldnames(data.EEG);
            tmpstr = {};
            for f=1:length(fields)
                if ischar(data.EEG.(fields{f}))
                    tmpstr = [tmpstr ; {sprintf('%-12s "%s"',fields{f}, data.EEG.(fields{f}))}]; %#ok<AGROW>
                elseif isinteger(data.EEG.(fields{f})) && length(data.EEG.(fields{f}))==1
                    tmpstr = [tmpstr ; {sprintf('%-12s %d',fields{f}, data.EEG.(fields{f}))}]; %#ok<AGROW>
                elseif isnumeric(data.EEG.(fields{f})) && length(data.EEG.(fields{f}))==1
                    tmpstr = [tmpstr ; {sprintf('%-12s %f',fields{f}, data.EEG.(fields{f}))}]; %#ok<AGROW>
                else
                    tmpstr = [tmpstr ; {sprintf('%-12s <%s>',fields{f}, class(data.EEG.(fields{f})))}]; %#ok<AGROW>
                end
            end
            app.listboxEegProperties.Items = tmpstr(:)';
        end

        % ------------------------------------------------------------------
        % Restrict two spectra (freq x chan) to the channels both have, matched
        % by label (case-insensitive, in the order of the first), and put the
        % second on the frequency axis of the first when they differ (e.g.
        % after resampling; frequencies outside its range become NaN).
        function [P, P2, fs, chanlocs, info] = matchPSDs(P, fs, chanlocs, P2, fs2, chanlocs2)
            lab1 = upper({chanlocs.labels});
            lab2 = upper({chanlocs2.labels});
            [~, i1, i2] = intersect(lab1, lab2, 'stable');
            info = sprintf('- comparing %d channels present in both (current %d, previous %d)', ...
                numel(i1), numel(lab1), numel(lab2));
            P  = P(:, i1);
            P2 = P2(:, i2);
            chanlocs = chanlocs(i1);
            if isempty(i1)
                return
            end
            if numel(fs) ~= numel(fs2) || any(abs(fs(:) - fs2(:)) > 1e-9)
                P2 = interp1(fs2(:), P2, fs(:), 'linear', NaN);
                info = sprintf('%s; previous spectrum interpolated to the current frequencies', info);
            end
        end

        % ------------------------------------------------------------------
        function [P, fs] = calc_PSD(app, EEG) %#ok<INUSL>
            if size(EEG.data,3)>1
                win = ones(1,size(EEG.data,2));
                [P, fs] = pfft(EEG.data(:,:)', EEG.srate, win, 0);
            else
                win = hanning(EEG.srate*2);
                [P, fs] = pfft(EEG.data(:,:)', EEG.srate, win, .5);
            end
        end

        % ------------------------------------------------------------------
        function skip_or_overwrite_dialog(app)
            fig = uifigure('Name', 'Action Required', 'Position', [500 500 300 150]);
            uilabel(fig, ...
                'Text', 'Some output files exists. What would you like to do?', ...
                'Position', [25 80 250 40], ...
                'HorizontalAlignment', 'center');
            uibutton(fig, 'push', ...
                'Text', 'Skip all', ...
                'Position', [50 30 80 30], ...
                'ButtonPushedFcn', @(btn,event) EegCallbacks.buttonCallback(app, fig, 'skip'));
            uibutton(fig, 'push', ...
                'Text', 'Overwrite', ...
                'Position', [170 30 80 30], ...
                'ButtonPushedFcn', @(btn,event) EegCallbacks.buttonCallback(app, fig, 'overwrite'));
            uiwait(fig);
        end

        % ------------------------------------------------------------------
        function openBigEditor(app)
            d = uifigure('Name','Edit text', ...
                         'Position',[0 0 520 400], ...
                         'WindowStyle','modal');
            movegui(d,'center');
            ta = uitextarea(d, ...
                'Position',[20 60 480 320], ...
                'Value', app.textareaExecuteCode.Value);
            uibutton(d,'push','Text','OK (Enter)', ...
                'Position',[320 18 85 30], ...
                'ButtonPushedFcn', @(~,~) commit());
            uibutton(d,'push','Text','Cancel (Esc)', ...
                'Position',[410 18 90 30], ...
                'ButtonPushedFcn', @(~,~) cancel());
            d.WindowKeyPressFcn = @onKey;
            d.CloseRequestFcn   = @(~,~) cancel();
            try focus(ta); catch, end
            uiwait(d);
            function onKey(~,evt)
                switch evt.Key
                    case 'return'
                        if isempty(evt.Modifier)
                            commit();
                        end
                    case 'escape'
                        cancel();
                end
            end
            function commit()
                lines = ta.Value;
                while numel(lines) > 1 && (isempty(lines{end}) || all(isspace(lines{end})))
                    lines(end) = [];
                end
                if isa(app.textareaExecuteCode,'matlab.ui.control.TextArea')
                    app.textareaExecuteCode.Value = lines;
                else
                    app.textareaExecuteCode.Value = strjoin(string(lines), ' ');
                end
                uiresume(d); delete(d);
            end
            function cancel()
                uiresume(d); delete(d);
            end
        end

        % ------------------------------------------------------------------
        function EEG = loadfile(app, PathName, FileName, listboxStdout, ...
                useBiosig, dataAreAvgRef)
            if PathName(end)~='/'
                PathName = [PathName '/'];
            end
            zz = strsplit(FileName,'.');
            switch zz{end}
                case 'cnt'
                    try
                        cLoadANTNeuro = true;
                        [EEG, cmd] = pop_loadeep_v4([PathName FileName], 'triggerfile', 'on');
                    catch E
                        if strcmpi(E.message, 'Error getting samples')
                            cLoadANTNeuro = false;
                            [EEG, cmd] = pop_loadcnt([PathName FileName] , 'dataformat', 'auto', 'memmapfile', '');
                        else
                            rethrow(E)
                        end
                    end
                    EEG.filename = [PathName FileName];
                    if cLoadANTNeuro
                        tmp = EEG;
                        evtcnt = 0;
                        for ev=1:length(tmp.event)
                            if tmp.event(ev).latency>1
                                evtcnt = evtcnt + 1;
                                if strcmp(tmp.event(ev).type,'__')
                                    typ = 'boundary';
                                else
                                    typ = strtrim(tmp.event(ev).type);
                                end
                                tmp.event(evtcnt).type = typ;
                                tmp.event(evtcnt).latency = EEG.event(ev).latency;
                                tmp.event(evtcnt).duration = EEG.event(ev).duration;
                            end
                        end
                        if length(tmp.event)>10 && ~isempty(str2num(tmp.event(1).type)) ... %#ok<ST2NM>
                                && tmp.event(1).latency>100*tmp.srate ...
                                && sum(strcmpi('boundary',{tmp.event.type}))==0
                            sig = tmp.data(:,(tmp.event(1).latency-tmp.srate*4):tmp.event(1).latency);
                            z = abs(zscore(mean(abs(diff(sig')')))');
                            ndx = tmp.event(1).latency-tmp.srate*4 + min(find(z>10)) + 0; %#ok<MXFND>
                            dummy = tmp.event(1);
                            dummy.type = 'boundary';
                            dummy.latency = ndx;
                            dummy.duration = 0;
                            tmp.event = cat(1,dummy,tmp.event(:));
                        end
                        EEG = tmp;
                        EegCallbacks.AddToListbox(app, listboxStdout, 'Read ANT CNT file')
                    else
                        EegCallbacks.AddToListbox(app, listboxStdout, 'Read Neuroscan CNT file')
                    end
                case 'set'
                    [EEG, cmd] = pop_loadset([PathName FileName]);
                    EEG.filename = [PathName FileName];
                    EegCallbacks.AddToListbox(app, listboxStdout, 'Read EEGLAB file')
                case 'bdf'
                    if useBiosig
                        EegCallbacks.AddToListbox(app, listboxStdout, 'Read BDF file with pop_biosig');
                        [EEG, cmd] = pop_biosig([PathName FileName], 'bdfeventmode',1);
                        EegCallbacks.AddToListbox(app, listboxStdout, ' - data read');
                    else
                        EegCallbacks.AddToListbox(app, listboxStdout, 'Read BDF file with pop_readbdf');
                        EegCallbacks.AddToListbox(app, listboxStdout, ' - read file header');
                        tmp = sopen([PathName FileName]);
                        EegCallbacks.AddToListbox(app, listboxStdout, sprintf(' - %d channels', tmp.NS));
                        [EEG, cmd] = pop_readbdf([PathName FileName], [], tmp.NS);
                        EegCallbacks.AddToListbox(app, listboxStdout, ' - data read');
                    end
                    EEG.filename = FileName;
                    EegCallbacks.AddToListbox(app, listboxStdout, ' - *** NOTE Data are raw. Please rereference in the next steps.');
                    EegCallbacks.AddToListbox(app, listboxStdout, ' - *** NOTE Renaming EXG* to EXT*.');
                    for ch=1:EEG.nbchan
                        if strncmp(EEG.chanlocs(ch).labels, "EXG", 3)
                            EEG.chanlocs(ch).labels(3) = "T";
                        end
                    end
                case 'edf'
                    [EEG, cmd] = pop_biosig([PathName FileName]);
                    EEG.filename = FileName;
                    EegCallbacks.AddToListbox(app, listboxStdout, 'Read EDF file.');
                case 'vhdr'
                    [EEG, cmd] = pop_loadbv(PathName, FileName, [], []);
                    EEG.filename = FileName;
                    EegCallbacks.AddToListbox(app, listboxStdout, 'Read BrainVision file.');
            end
            if dataAreAvgRef
                EegCallbacks.AddToListbox(app, listboxStdout, 'Assuming data are in average reference.');
                EEG.ref = 'average';
                test = mean(EEG.data(:,100:199));
                if any(test>.5)
                    warning('The data do not average out to near zero. check AvgRef setting.')
                    EegCallbacks.AddToListbox(app, listboxStdout, '*** The data do not average out to near zero.\n*** check AvgRef setting.');
                end
            end
            EEG.history = [EEG.history char(uint8(10)) cmd];
        end

        % ------------------------------------------------------------------
        % The workflow buttons of the left-hand column, top to bottom, EXCEPT
        % Open, which is coloured on its own (see ColIdle/ColReady/ColBusy).
        % Buttons in the other columns are not part of the scheme.
        function names = mainButtonNames(~)
            names = {'pushbuttonFlatline','pushbuttonExcessive','pushbuttonChanlocs', ...
                     'pushbuttonResample','pushbuttonRereference','pushbuttonFilter', ...
                     'pushbuttonInitialICA','pushbuttonAltEOG','pushbuttonEMG', ...
                     'pushbuttonAltEMG','pushbuttonASR','pushbuttonICA','pushbuttonLineNoise'};
        end

        % ------------------------------------------------------------------
        function setMainButtonsColor(app, color)
            names = EegCallbacks.mainButtonNames(app);
            for b = 1:numel(names)
                if isprop(app, names{b}) && isvalid(app.(names{b}))
                    app.(names{b}).BackgroundColor = color;
                end
            end
        end

        % ------------------------------------------------------------------
        % Startup state: nothing is loaded, so only Open can be used.
        function resetButtonColors(app)
            EegCallbacks.setMainButtonsColor(app, EegCallbacks.ColIdle);
            if isprop(app, 'pushbuttonOpen') && isvalid(app.pushbuttonOpen)
                app.pushbuttonOpen.BackgroundColor = EegCallbacks.ColReady;
            end
        end

    end % methods (Static, Access = private)

end % classdef
