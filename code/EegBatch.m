classdef EegBatch
% EegBatch  The batch window, built in code instead of from a .fig.
%
%   EegBatch.dialog(parentFig, SETTINGSDIR, DEFAULTDIR, steps, onStart)
%
%   steps    cell array of step labels, in the order of the caller's switch:
%            the label at index k is step k. Separator labels ('--------')
%            are offered like any other item but do nothing.
%   onStart  called when Start is pressed, with one struct:
%              .files      selected files, full paths
%              .folder     folder of the first file
%              .prefix     prefix of the output file names
%              .outputdir  folder the results are written to
%              .steps      chosen step indices, in the order of the slots
%              .labels     their labels
%
%   The window returns at once and does NOT use uiwait: it is opened from an
%   App Designer callback, and a uifigure cannot dispatch its own buttons
%   while the calling callback is still on the stack. So the batch run starts
%   from Start, not after the window closes.
%
%   What the user set is remembered in <SETTINGSDIR>/EegBatch.ini (prefix,
%   folders, and the chosen steps BY LABEL, so that reordering
%   the step list later does not scramble a saved recipe) and the file list
%   in <SETTINGSDIR>/.EegBatch_Filenames.ini.

    properties (Constant)
        NSLOTS  = 16;                     % number of step slots
        SELECT  = '<select>';             % first item of every slot
    end

    methods (Static)

        % ------------------------------------------------------------------
        function dialog(parentFig, SETTINGSDIR, DEFAULTDIR, steps, onStart)
            old = findall(groot, 'Type', 'figure', 'Tag', 'EegBatchDialog');
            if ~isempty(old)
                figure(old(1));
                return
            end
            steps = cellstr(steps(:));
            items = [{EegBatch.SELECT}; steps];

            st = EegBatch.loadState(SETTINGSDIR);
            if isempty(st.folder) || ~isfolder(st.folder)
                st.folder = DEFAULTDIR;
            end
            if isempty(st.outputdir)
                st.outputdir = st.folder;
            end

            % height from the layout: padding, 4 form rows, the step slots with
            % their spacing, the button row, and the spacing between the rows
            pad = 10; formH = 28; slotH = 26; slotGap = 4; btnH = 34; rowGap = 6;
            figW = 430;
            figH = 2*pad + 4*formH + EegBatch.NSLOTS*slotH + (EegBatch.NSLOTS-1)*slotGap ...
                   + btnH + 5*rowGap + 8;
            scr  = get(groot, 'ScreenSize');
            figH = min(figH, scr(4) - 120);     % small screens: the slots scroll
            pos  = [100 100 figW figH];
            try
                pp = parentFig.Position;
                pos(1) = pp(1) + pp(3) + 20;
                pos(2) = pp(2) + max(0, pp(4) - figH);
            catch
            end
            pos(1) = min(max(pos(1), 20), max(20, scr(3) - figW - 20));
            pos(2) = min(max(pos(2), 40), max(40, scr(4) - figH - 60));

            fig = uifigure('Name', 'Run batch', 'Position', pos, 'Tag', 'EegBatchDialog');
            outer = uigridlayout(fig, [6 1]);
            outer.RowHeight  = {formH, formH, formH, formH, '1x', btnH};
            outer.RowSpacing = rowGap;
            outer.Padding    = [pad pad pad pad];

            h = struct();
            % ---- files -------------------------------------------------------
            r = EegBatch.formRow(outer, 1, 'Files');
            h.count = uilabel(r, 'Text', sprintf('%d files', numel(st.files)), 'Tag', 'count');
            h.count.Layout.Column = 2;
            b = uibutton(r, 'Text', 'Select...', 'ButtonPushedFcn', @(~,~) pickFiles());
            b.Layout.Column = 3;

            r = EegBatch.formRow(outer, 2, 'Folder');
            h.folder = uieditfield(r, 'Value', st.folder, 'Editable', 'off', 'Tag', 'folder');
            h.folder.Layout.Column = 2;

            r = EegBatch.formRow(outer, 3, 'Prefix');
            h.prefix = uieditfield(r, 'Value', st.prefix, 'Tag', 'prefix');
            h.prefix.Layout.Column = 2;
            h.prefix.Tooltip = 'Put in front of every output file name';

            r = EegBatch.formRow(outer, 4, 'Output folder');
            h.outdir = uieditfield(r, 'Value', st.outputdir, 'Tag', 'outputdir');
            h.outdir.Layout.Column = 2;
            b = uibutton(r, 'Text', 'Select...', 'ButtonPushedFcn', @(~,~) pickOutput());
            b.Layout.Column = 3;

            % ---- steps -------------------------------------------------------
            g = uigridlayout(outer, [EegBatch.NSLOTS 2]);
            g.Layout.Row  = 5;
            g.RowHeight   = repmat({slotH}, 1, EegBatch.NSLOTS);
            g.ColumnWidth = {60, '1x'};
            g.RowSpacing  = slotGap;
            g.Padding     = [0 0 0 0];
            g.Scrollable  = 'on';
            h.slot = gobjects(1, EegBatch.NSLOTS);
            for k = 1:EegBatch.NSLOTS
                lbl = uilabel(g, 'Text', sprintf('%d', k), 'HorizontalAlignment', 'right');
                lbl.Layout.Row = k; lbl.Layout.Column = 1;
                value = EegBatch.SELECT;
                if k <= numel(st.labels) && ~isempty(st.labels{k}) && any(strcmp(items, st.labels{k}))
                    value = st.labels{k};
                end
                h.slot(k) = uidropdown(g, 'Items', items, 'Value', value, 'Tag', sprintf('step%02d', k));
                h.slot(k).Layout.Row = k; h.slot(k).Layout.Column = 2;
            end

            % ---- buttons -----------------------------------------------------
            bg = uigridlayout(outer, [1 3]);
            bg.Layout.Row  = 6;
            bg.ColumnWidth = {'1x', 90, 90};
            bg.Padding     = [0 0 0 0];
            uilabel(bg, 'Text', '');
            uibutton(bg, 'Text', 'Cancel', 'ButtonPushedFcn', @(~,~) closeMe());
            uibutton(bg, 'Text', 'Start',  'ButtonPushedFcn', @(~,~) startMe());
            fig.CloseRequestFcn = @(~,~) closeMe();

            files = st.files;

            % ---- nested callbacks --------------------------------------------
            function pickFiles()
                start = h.folder.Value;
                if ~isfolder(start), start = DEFAULTDIR; end
                % '*' plus the REFilter below = every readable type; the
                % specific entries narrow it down. Folders are always listed
                % (redirs off), so the list stays navigable by double click.
                picked = uipickfiles('filterspec', start, ...
                    'type', {'*','All readable EEG'; '*.bdf','Biosemi'; '*.cnt','Neuroscan or ANT'; ...
                             '*.edf','European data format'; '*.set','EEGLAB'; '*.vhdr','BrainVision'}, ...
                    'REFilter', '(\.bdf|\.cnt|\.edf|\.set|\.vhdr)$', 'redirs', false, ...
                    'prompt', 'Select files for processing', 'output', 'cell');
                EegCallbacks.bringToFront(fig);
                if isequal(picked, 0) || isempty(picked)
                    return                      % cancelled or nothing picked
                end
                files = picked(:)';
                folder = fileparts(files{1});
                h.folder.Value = folder;
                h.count.Text = sprintf('%d files', numel(files));
                if isempty(h.outdir.Value)
                    h.outdir.Value = folder;
                end
                EegBatch.rememberFolder(parentFig, SETTINGSDIR, folder);
            end

            function pickOutput()
                start = h.outdir.Value;
                if ~isfolder(start), start = h.folder.Value; end
                if ~isfolder(start), start = DEFAULTDIR; end
                out = uigetdir(start, 'Pick an output folder');
                EegCallbacks.bringToFront(fig);
                if ischar(out)                  % 0 when cancelled
                    h.outdir.Value = out;
                end
            end

            function sel = collect()
                sel.files     = files;
                sel.folder    = h.folder.Value;
                sel.prefix    = h.prefix.Value;
                sel.outputdir = h.outdir.Value;
                sel.steps     = [];
                sel.labels    = {};
                for s = 1:EegBatch.NSLOTS
                    v = find(strcmp(items, h.slot(s).Value), 1) - 1;   % 0 = <select>
                    if v > 0
                        sel.steps(end+1)  = v; %#ok<AGROW>
                        sel.labels{end+1} = items{v+1}; %#ok<AGROW>
                    end
                end
            end

            function startMe()
                sel = collect();
                if isempty(sel.files)
                    uialert(fig, 'No files selected.', 'Batch');
                    return
                end
                if isempty(sel.outputdir) || ~isfolder(sel.outputdir)
                    uialert(fig, 'The output folder does not exist.', 'Batch');
                    return
                end
                if isempty(sel.steps)
                    uialert(fig, 'No steps selected.', 'Batch');
                    return
                end
                EegBatch.saveState(SETTINGSDIR, sel, h.slot, items);
                delete(fig);
                if isa(onStart, 'function_handle')
                    onStart(sel);
                end
            end

            function closeMe()
                EegBatch.saveState(SETTINGSDIR, collect(), h.slot, items);
                delete(fig);
            end
        end

        % ------------------------------------------------------------------
        % One 'label | field | button' row of the form part.
        function r = formRow(parent, row, label)
            r = uigridlayout(parent, [1 3]);
            r.Layout.Row  = row;
            r.ColumnWidth = {90, '1x', 80};
            r.Padding     = [0 0 0 0];
            r.ColumnSpacing = 6;
            lbl = uilabel(r, 'Text', label, 'HorizontalAlignment', 'right');
            lbl.Layout.Column = 1;
        end

        % ------------------------------------------------------------------
        function FN = iniFile(SETTINGSDIR)
            FN = fullfile(SETTINGSDIR, 'EegBatch.ini');
        end

        function FN = filesFile(SETTINGSDIR)
            FN = fullfile(SETTINGSDIR, '.EegBatch_Filenames.ini');
        end

        % ------------------------------------------------------------------
        % What was set last time. Unknown keys and labels that no longer exist
        % are ignored, so the step list can change between versions.
        function st = loadState(SETTINGSDIR)
            st = struct('prefix', 'c_', 'folder', '', ...
                        'outputdir', '', 'labels', {{}}, 'files', {{}});
            try
                FN = EegBatch.iniFile(SETTINGSDIR);
                if exist(FN, 'file')
                    opts = detectImportOptions(FN, 'TextType', 'string', 'filetype', 'text', 'Delimiter', '\t');
                    opts.DataLines = [2 Inf];
                    opts.VariableTypes(:) = {'string'};
                    T = readtable(FN, opts);
                    T.Properties.VariableNames = lower(T.Properties.VariableNames);
                    val = @(k) EegBatch.lookup(T, k);
                    for f = {'prefix','folder','outputdir'}
                        v = val(f{1});
                        if ~isempty(v), st.(f{1}) = v; end
                    end
                    st.labels = cell(1, EegBatch.NSLOTS);
                    for k = 1:EegBatch.NSLOTS
                        st.labels{k} = val(sprintf('step%02d', k));
                    end
                end
            catch
                warning('EegBatch:readFailed', 'Could not read the batch settings. Using defaults.');
            end
            fid = fopen(EegBatch.filesFile(SETTINGSDIR), 'r');
            if fid > 0
                line = fgetl(fid);
                while ischar(line)
                    if ~isempty(strtrim(line))
                        st.files{end+1} = strtrim(line);
                    end
                    line = fgetl(fid);
                end
                fclose(fid);
            end
        end

        % ------------------------------------------------------------------
        function saveState(SETTINGSDIR, sel, slots, items) %#ok<INUSD>
            try
                if ~isfolder(SETTINGSDIR)
                    mkdir(SETTINGSDIR);
                end
                key = {'prefix'; 'folder'; 'outputdir'};
                val = {sel.prefix; sel.folder; sel.outputdir};
                for k = 1:numel(slots)
                    key{end+1} = sprintf('step%02d', k); %#ok<AGROW>
                    v = slots(k).Value;
                    if strcmp(v, EegBatch.SELECT)
                        v = '';
                    end
                    val{end+1} = v; %#ok<AGROW>
                end
                writetable(table(key, val, 'VariableNames', {'key','val'}), ...
                    EegBatch.iniFile(SETTINGSDIR), 'Delimiter', '\t', 'FileType', 'text');
                fid = fopen(EegBatch.filesFile(SETTINGSDIR), 'w');
                if fid > 0
                    for f = 1:numel(sel.files)
                        fprintf(fid, '%s\n', sel.files{f});
                    end
                    fclose(fid);
                end
            catch
                warning('EegBatch:writeFailed', 'Could not save the batch settings.');
            end
        end

        % ------------------------------------------------------------------
        function v = lookup(T, key)
            v = '';
            ndx = find(strcmpi(T.key, key), 1);
            if ~isempty(ndx) && ~ismissing(T.val(ndx))
                v = char(T.val(ndx));
            end
        end

        % ------------------------------------------------------------------
        % Keep the folder as the one file dialogs start in, here and in the
        % main window.
        function rememberFolder(parentFig, SETTINGSDIR, folder)
            try
                d = EegCallbacks.saveDefaultDir(SETTINGSDIR, folder);
                pdata = guidata(parentFig);
                pdata.DEFAULTDIR = d;
                guidata(parentFig, pdata);
            catch
            end
        end

    end
end
