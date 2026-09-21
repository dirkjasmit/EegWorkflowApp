classdef EegParams
% EegParams  Definition-driven settings store for the EEG workflow app.
%
%   Every editable setting is declared in EegWorkflow_parameters.xlsx, one
%   row per (button, parname):
%
%       button          key of the workflow button the setting belongs to
%       parname         parameter name, unique within a button
%       pardescription  label shown in the edit dialog
%       partype         'checkbox' | 'value' | 'dropdown'
%       paritems        'value'    -> 'min=..,max=..,by=..'
%                       'dropdown' -> '"item,item,item"'
%                       'checkbox' -> (empty)
%       tooltip         tooltip shown on the control
%       initialpardefault
%                       the value used the FIRST time the app runs, i.e. only
%                       until a settings file exists.  When the column is
%                       absent the fallback is the minimum ('value'), the
%                       first item ('dropdown') or false ('checkbox').
%
%   Settings live in two files with the same rows:
%
%     EegWorkflow_parameters.xlsx      in resources/; declares the
%                                      parameters and their initial values
%     <SETTINGSDIR>/EegWorkflow_parameters.ini
%                                      tab separated button/parname/value, in
%                                      the OS-specific settings directory;
%                                      written on every OK and on close, and
%                                      it is what is read back on start-up.
%
%   At load time a 'value' column is appended that holds the current setting
%   as a string.  The whole table lives in guidata as data.params.
%
%   Typical use inside a callback:
%
%       data = guidata(hObject);
%       P    = @(b,p) EegParams.get(data.params, b, p);
%       lo   = P('filter','low');          % double
%       useB = P('open','biosig');         % logical
%       ref  = P('rereference','refchans');% char

    methods (Static)

        % ------------------------------------------------------------------
        % Full path of the spreadsheet that declares the parameters.
        function FN = definitionFile()
            % this class lives in code/, the spreadsheet in resources/
            here = fileparts(which('EegParams'));
            if isempty(here)
                here = pwd;
            else
                here = fullfile(fileparts(here), 'resources');
            end
            FN = fullfile(here, 'EegWorkflow_parameters.xlsx');
        end

        % ------------------------------------------------------------------
        % Full path of the file holding the saved values.
        function FN = valueFile(SETTINGSDIR)
            FN = fullfile(SETTINGSDIR, 'EegWorkflow_parameters.ini');
        end

        % ------------------------------------------------------------------
        % Read the definitions, add the 'value' column and apply saved values.
        function T = load(SETTINGSDIR)
            FN = EegParams.definitionFile();
            opts = detectImportOptions(FN, 'FileType', 'spreadsheet', ...
                                           'VariableNamingRule', 'preserve');
            % Never let MATLAB guess where the table starts: an empty cell
            % (e.g. a missing initialpardefault) can make it pick a later
            % row, silently dropping every parameter above it.
            opts.VariableNamesRange = 'A1';
            opts.DataRange          = 'A2';
            opts = setvartype(opts, opts.VariableNames, 'string');
            T = readtable(FN, opts);

            T.Properties.VariableNames = lower(T.Properties.VariableNames);
            required = {'button','parname','pardescription','partype','paritems'};
            for r = 1:numel(required)
                if ~ismember(required{r}, T.Properties.VariableNames)
                    error('EegParams:missingColumn', ...
                          'Column ''%s'' is missing from %s', required{r}, FN);
                end
            end
            if ~ismember('tooltip', T.Properties.VariableNames)
                T.tooltip = repmat("", height(T), 1);
            end

            % blanks read back as <missing>; make them plain empty strings
            for v = 1:width(T)
                col = T.(v);
                col(ismissing(col)) = "";
                T.(v) = strtrim(col);
            end

            % drop rows without a button or a parameter name
            T(T.button=="" | T.parname=="", :) = [];

            % initial values: only ever used the first time the app runs, i.e.
            % until the settings file in SETTINGSDIR exists.
            T.value = repmat("", height(T), 1);
            hasDefault = ismember('initialpardefault', T.Properties.VariableNames);
            for i = 1:height(T)
                if hasDefault && T.initialpardefault(i) ~= ""
                    T.value(i) = T.initialpardefault(i);
                else
                    T.value(i) = EegParams.defaultValue(T.partype(i), T.paritems(i));
                end
            end

            T = EegParams.applySaved(T, SETTINGSDIR);
        end

        % ------------------------------------------------------------------
        % Overwrite the 'value' column with whatever was saved last session.
        % Unknown or stale (button,parname) pairs in the file are ignored.
        function T = applySaved(T, SETTINGSDIR)
            FN = EegParams.valueFile(SETTINGSDIR);
            if ~exist(FN, 'file')
                return
            end
            try
                opts = detectImportOptions(FN, 'FileType', 'text', ...
                                               'Delimiter', '\t', ...
                                               'VariableNamingRule', 'preserve');
                opts.VariableNamesLine = 1;
                opts.DataLines         = [2 Inf];
                opts = setvartype(opts, opts.VariableNames, 'string');
                S = readtable(FN, opts);
                S.Properties.VariableNames = lower(S.Properties.VariableNames);
                for i = 1:height(S)
                    ndx = EegParams.rowIndex(T, S.button(i), S.parname(i));
                    if ~isempty(ndx)
                        val = S.value(i);
                        if ismissing(val)
                            val = "";
                        end
                        T.value(ndx) = val;
                    end
                end
            catch
                warning('EegParams:readFailed', ...
                        'Could not read %s. Falling back to defaults.', FN);
            end
            T = EegParams.validate(T);
        end

        % ------------------------------------------------------------------
        % Write the current values.  Called from the figure CloseRequestFcn.
        function save(T, SETTINGSDIR)
            if ~exist(SETTINGSDIR, 'dir')
                mkdir(SETTINGSDIR);
            end
            S = table(T.button, T.parname, T.value, ...
                      'VariableNames', {'button','parname','value'});
            writetable(S, EegParams.valueFile(SETTINGSDIR), ...
                       'Delimiter', '\t', 'FileType', 'text');
        end

        % ------------------------------------------------------------------
        % Clamp/repair values that no longer fit their definition (e.g. after
        % the spreadsheet changed the limits or the dropdown items).
        function T = validate(T)
            for i = 1:height(T)
                switch lower(T.partype(i))
                    case 'checkbox'
                        if ~ismember(lower(T.value(i)), ["0","1","true","false"])
                            T.value(i) = "0";
                        end
                    case 'value'
                        [mn, mx] = EegParams.valueLimits(T.paritems(i));
                        v = str2double(T.value(i));
                        if isnan(v)
                            v = mn;
                        end
                        T.value(i) = string(min(max(v, mn), mx));
                    case 'dropdown'
                        items = EegParams.dropdownItems(T.paritems(i));
                        if isempty(items)
                            continue
                        end
                        if ~any(strcmpi(items, T.value(i)))
                            T.value(i) = string(items{1});
                        end
                end
            end
        end

        % ------------------------------------------------------------------
        % Row index of one setting, or [] when it is not declared.
        function ndx = rowIndex(T, button, parname)
            ndx = find(strcmpi(T.button, button) & strcmpi(T.parname, parname), 1);
        end

        % ------------------------------------------------------------------
        % Current value of one setting, converted to its natural MATLAB type:
        %   checkbox -> logical, value -> double, dropdown -> char
        function v = get(T, button, parname)
            ndx = EegParams.rowIndex(T, button, parname);
            if isempty(ndx)
                error('EegParams:unknownParameter', ...
                      'No parameter ''%s'' declared for button ''%s''.', parname, button);
            end
            switch lower(T.partype(ndx))
                case 'checkbox'
                    s = lower(T.value(ndx));
                    v = (s == "1") || (s == "true");
                case 'value'
                    v = str2double(T.value(ndx));
                otherwise
                    v = char(T.value(ndx));
            end
        end

        % ------------------------------------------------------------------
        % Index of the current item of a dropdown (1-based), 0 when unknown.
        % Handy where the old code switched on popupmenu.Value.
        function n = getIndex(T, button, parname)
            ndx = EegParams.rowIndex(T, button, parname);
            n = 0;
            if isempty(ndx)
                return
            end
            items = EegParams.dropdownItems(T.paritems(ndx));
            hit = find(strcmpi(items, T.value(ndx)), 1);
            if ~isempty(hit)
                n = hit;
            end
        end

        % ------------------------------------------------------------------
        % Copy just one button's values from Tedited into T. Used when a
        % settings dialog is accepted: several dialogs can be open at once,
        % so writing back a whole captured table would revert whatever
        % another dialog changed in the meantime.
        function T = mergeButton(T, Tedited, buttonKey)
            rows = find(strcmpi(Tedited.button, buttonKey));
            for r = rows(:)'
                ndx = EegParams.rowIndex(T, Tedited.button(r), Tedited.parname(r));
                if ~isempty(ndx)
                    T.value(ndx) = Tedited.value(r);
                end
            end
        end

        % ------------------------------------------------------------------
        function T = set(T, button, parname, value)
            ndx = EegParams.rowIndex(T, button, parname);
            if isempty(ndx)
                error('EegParams:unknownParameter', ...
                      'No parameter ''%s'' declared for button ''%s''.', parname, button);
            end
            if islogical(value)
                value = double(value);
            end
            if isnumeric(value)
                value = string(value);
            end
            T.value(ndx) = string(value);
        end

        % ------------------------------------------------------------------
        function s = defaultValue(partype, paritems)
            switch lower(partype)
                case 'checkbox'
                    s = "0";
                case 'value'
                    mn = EegParams.valueLimits(paritems);
                    s = string(mn);
                case 'dropdown'
                    items = EegParams.dropdownItems(paritems);
                    if isempty(items)
                        s = "";
                    else
                        s = string(items{1});
                    end
                otherwise
                    s = "";
            end
        end

        % ------------------------------------------------------------------
        % Parse 'min=0.2,max=10.0,by=0.2'.  Missing entries get safe defaults.
        function [mn, mx, step] = valueLimits(paritems)
            mn = 0; mx = 100; step = 1;
            s = char(paritems);
            t = regexp(s, 'min\s*=\s*(-?[\d.eE+-]+)', 'tokens', 'once');
            if ~isempty(t), mn = str2double(t{1}); end
            t = regexp(s, 'max\s*=\s*(-?[\d.eE+-]+)', 'tokens', 'once');
            if ~isempty(t), mx = str2double(t{1}); end
            t = regexp(s, 'by\s*=\s*(-?[\d.eE+-]+)', 'tokens', 'once');
            if ~isempty(t), step = str2double(t{1}); end
            if isnan(mn), mn = 0;   end
            if isnan(mx), mx = 100; end
            if isnan(step) || step<=0, step = 1; end
            if mx < mn
                [mn, mx] = deal(mx, mn);
            end
        end

        % ------------------------------------------------------------------
        % Parse '"a,b,c"' (the surrounding quotes are optional) into a cellstr.
        function items = dropdownItems(paritems)
            s = strtrim(char(paritems));
            s = regexprep(s, '^"|"$', '');          % strip wrapping quotes
            if isempty(s)
                items = {};
                return
            end
            items = strtrim(strsplit(s, ','));
            items = items(~cellfun(@isempty, items));
        end

        % ------------------------------------------------------------------
        % The '<key>=<value>' string shown in a listbox.
        function s = displayItem(T, i)
            switch lower(T.partype(i))
                case 'checkbox'
                    val = "false";
                    if EegParams.get(T, T.button(i), T.parname(i))
                        val = "true";
                    end
                case 'value'
                    val = string(sprintf('%g', str2double(T.value(i))));
                otherwise
                    val = T.value(i);
            end
            s = char(T.parname(i) + "=" + val);
        end

        % ------------------------------------------------------------------
        % Listbox component name -> spreadsheet button key -> dialog title.
        % Add a row here when a button with settings is added to the window.
        function m = listboxMap()
            m = { 'ListBoxOpen',        'open',             'Open'
                  'ListBoxFlatline',    'flatline',         'Flatline'
                  'ListBoxExcessive',   'excessive signal', 'Excessive signal'
                  'ListBoxLookup',      'lookup',           'Lookup'
                  'ListBoxResample',    'resample',         'Resample'
                  'ListBoxRereference', 'rereference',      'Rereference'
                  'ListBoxFilter',      'filter',           'Filter'
                  'ListBoxCleanline',   'cleanline',        'Line noise'
                  'ListBoxEOG',         'eog',              'EOG'
                  'ListBoxAltEOG',      'alt eog',          'Alt EOG'
                  'ListBoxEMG',         'emg',              'EMG'
                  'ListBoxAltEMG',      'alt emg',          'Alt EMG'
                  'ListBoxASR',         'asr',              'ASR'
                  'ListBoxICA',         'ICA',              'ICA'
                  'ListBoxFlatPeriods',          'flat periods',        'Flat periods'
                  'ListBoxExcessivePeriods',     'excessive periods',   'Excessive periods'
                  'ListBoxEMGPeriods',           'EMG periods',         'EMG periods'
                  'ListBoxInterpolationPeriods', 'interpolation clean', 'Interpolation clean' };
        end

        % ------------------------------------------------------------------
        % Fill every listbox with the '<key>=<value>' lines of its button.
        function refresh(app, T)
            m = EegParams.listboxMap();
            for k = 1:size(m,1)
                EegParams.refreshOne(app, T, m{k,1}, m{k,2});
            end
        end

        % ------------------------------------------------------------------
        function refreshOne(app, T, lbName, buttonKey)
            if ~isprop(app, lbName) || ~isvalid(app.(lbName))
                return
            end
            rows = find(strcmpi(T.button, buttonKey));
            items = cell(1, numel(rows));
            for r = 1:numel(rows)
                items{r} = EegParams.displayItem(T, rows(r));
            end
            lb = app.(lbName);
            lb.Items = items;
            lb.ItemsData = num2cell(rows(:)');
            lb.Value = {};                     % nothing selected: display only
            lb.Tooltip = 'Click to edit these settings';
        end

        % ------------------------------------------------------------------
        % Editor for the settings of one button.
        %
        % ASYNCHRONOUS on purpose: it builds the window and returns at once,
        % rather than blocking on uiwait. The caller is itself inside an App
        % Designer callback, and a uiwait there keeps that callback on the
        % stack, so the dialog's own OK/Cancel presses never get dispatched
        % and the window appears dead. onAccept(Tnew) is invoked when OK is
        % pressed; Cancel and the window's X just close it.
        function editDialog(parentFig, T, buttonKey, dlgTitle, onAccept)
            rows = find(strcmpi(T.button, buttonKey));
            if isempty(rows)
                return
            end
            n = numel(rows);

            % never leave two editors open for the same button
            old = findall(groot, 'Type', 'figure', 'Tag', ['EegParamsDialog_' buttonKey]);
            if ~isempty(old)
                figure(old(1));
                return
            end

            rowH   = 28;                        % one parameter row
            rowGap = 4;                         % spacing between parameter rows
            btnH   = 34;                        % OK / Cancel row
            pad    = 10;                        % window padding
            figW   = 420;
            figH   = 2*pad + n*rowH + (n-1)*rowGap + pad + btnH;
            % never taller than the screen; the parameters scroll instead
            scr = get(groot, 'ScreenSize');
            figH = min(figH, scr(4) - 120);
            pos  = [100 100 figW figH];
            try
                pp = parentFig.Position;
                pos(1) = pp(1) + round((pp(3)-figW)/2);
                pos(2) = pp(2) + round((pp(4)-figH)/2);
            catch
            end
            pos(1) = min(max(pos(1), 20), max(20, scr(3) - figW - 20));
            pos(2) = min(max(pos(2), 40), max(40, scr(4) - figH - 60));

            fig = uifigure('Name', dlgTitle, 'Position', pos, ...
                           'Tag', ['EegParamsDialog_' buttonKey], ...
                           'WindowStyle', 'alwaysontop', 'Resize', 'on');

            % parameters in a scrollable grid, buttons in a fixed row below it
            outer = uigridlayout(fig, [2 1]);
            outer.RowHeight  = {'1x', btnH};
            outer.RowSpacing = pad;
            outer.Padding    = [pad pad pad pad];

            g = uigridlayout(outer, [n 2]);
            g.Layout.Row  = 1;
            g.RowHeight   = repmat({rowH}, 1, n);
            g.ColumnWidth = {180, '1x'};
            g.RowSpacing  = rowGap;
            g.Padding     = [0 0 0 0];
            g.Scrollable  = 'on';

            ctrls = gobjects(1, n);
            for r = 1:n
                i    = rows(r);
                desc = char(T.pardescription(i));
                if isempty(desc)
                    desc = char(T.parname(i));
                end
                tip = char(T.tooltip(i));

                switch lower(T.partype(i))
                    case 'checkbox'
                        c = uicheckbox(g, 'Text', desc, ...
                                          'Value', EegParams.get(T, T.button(i), T.parname(i)));
                        c.Layout.Row    = r;
                        c.Layout.Column = [1 2];

                    case 'value'
                        lbl = uilabel(g, 'Text', desc, 'HorizontalAlignment', 'right');
                        lbl.Layout.Row = r; lbl.Layout.Column = 1;
                        [mn, mx, st] = EegParams.valueLimits(T.paritems(i));
                        v = min(max(EegParams.get(T, T.button(i), T.parname(i)), mn), mx);
                        % A spinner with Editable off is exactly a read-only
                        % numeric box driven by its up/down arrows.
                        c = uispinner(g, 'Limits', [mn mx], 'Step', st, ...
                                         'Value', v, 'Editable', 'off', ...
                                         'ValueDisplayFormat', EegParams.numFormat(st));
                        c.Layout.Row = r; c.Layout.Column = 2;
                        if ~isempty(tip), lbl.Tooltip = tip; end

                    case 'dropdown'
                        lbl = uilabel(g, 'Text', desc, 'HorizontalAlignment', 'right');
                        lbl.Layout.Row = r; lbl.Layout.Column = 1;
                        items = EegParams.dropdownItems(T.paritems(i));
                        if isempty(items)
                            items = {char(T.value(i))};
                        end
                        cur = char(T.value(i));
                        if ~any(strcmpi(items, cur))
                            cur = items{1};
                        end
                        c = uidropdown(g, 'Items', items, 'Value', cur);
                        c.Layout.Row = r; c.Layout.Column = 2;
                        if ~isempty(tip), lbl.Tooltip = tip; end

                    otherwise
                        lbl = uilabel(g, 'Text', desc, 'HorizontalAlignment', 'right');
                        lbl.Layout.Row = r; lbl.Layout.Column = 1;
                        c = uilabel(g, 'Text', char(T.value(i)));
                        c.Layout.Row = r; c.Layout.Column = 2;
                end
                if ~isempty(tip)
                    c.Tooltip = tip;
                end
                ctrls(r) = c;
            end

            bg = uigridlayout(outer, [1 3]);
            bg.Layout.Row    = 2;
            bg.Layout.Column = 1;
            bg.ColumnWidth   = {'1x', 90, 90};
            bg.Padding       = [0 0 0 0];
            uilabel(bg, 'Text', '');
            uibutton(bg, 'Text', 'Cancel', 'ButtonPushedFcn', @onCancel);
            uibutton(bg, 'Text', 'OK',     'ButtonPushedFcn', @onOk);
            fig.CloseRequestFcn = @onCancel;

            % no uiwait: the window lives on and its callbacks do the work

            function onOk(~,~)
                try
                    for k = 1:n
                        j = rows(k);
                        if ~isvalid(ctrls(k))
                            continue
                        end
                        switch lower(T.partype(j))
                            case 'checkbox'
                                T.value(j) = string(double(ctrls(k).Value));
                            case {'value','dropdown'}
                                T.value(j) = string(ctrls(k).Value);
                        end
                    end
                    if isa(onAccept, 'function_handle')
                        onAccept(T);
                    end
                catch ME
                    EegParams.closeIfValid(fig);
                    rethrow(ME);
                end
                EegParams.closeIfValid(fig);
            end

            function onCancel(~,~)
                EegParams.closeIfValid(fig);
            end
        end

        % ------------------------------------------------------------------
        function fmt = numFormat(step)
            dec = 0;
            if step > 0 && step < 1
                dec = max(0, ceil(-log10(step)));
            end
            if dec == 0
                fmt = '%d';
            else
                fmt = ['%.' num2str(dec) 'f'];
            end
        end

        % ------------------------------------------------------------------
        function closeIfValid(fig)
            if isvalid(fig)
                delete(fig);
            end
        end

    end
end
