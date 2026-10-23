classdef EegPreset
% EegPreset  Pipeline presets of the batch window, kept as files.
%
%   P     = EegPreset.list(SETTINGSDIR)        every preset found, in order
%   P     = EegPreset.load(file)               one preset
%   EegPreset.apply(app, P, say)               write its settings and code,
%                                              logging through say(fmt, ...)
%   EegPreset.saveCurrent(app, name, labels)   the current state as a preset
%
% A preset is a JSON file. The shipped ones live in resources/presets next to
% this class's repository; a user's own in <SETTINGSDIR>/presets, and those
% are the only ones Save writes to (and the only ones that can be deleted).
%
%   {
%     "name":     "DISCOVER-EEG",            shown on the button
%     "order":    3,                         button order, lower is first
%     "requires": "PrepPipeline",            EEGLAB plugin needed (optional)
%     "steps":    ["Lookup", "Resample"],    batch step labels, in order
%     "settings": { "resample": { "srate": "250" } },
%                                            button -> parameter -> value, as
%                                            in EegWorkflow_parameters.xlsx
%     "controls": { "ICAtypeDropDown": "runica" },
%                                            main-window controls to set
%     "code":     ["EEG = ...;"],            fills the Execute code box
%     "batch":    { "prefix": "disc_" },     batch-window fields
%     "check":    ["Lookup", "Line frequency"],
%                                            what the user must still set
%     "notes":    ["..."]                    logged in the output pane
%   }
%
% Every field except name is optional: a preset that only lists steps is
% fine, and so is one that only writes settings. Settings that no longer
% exist in the spreadsheet are reported and skipped, so an old preset keeps
% working.
%
% Written for EegWorkflowApp. D.J.A. Smit & Claude (Anthropic), 2026.

    properties (Constant)
        FOLDER = 'presets';        % under the repository root and SETTINGSDIR
        EXT    = '.json';
    end

    methods (Static)

        % ------------------------------------------------------------------
        % Folders holding presets: the shipped ones first, then the user's.
        function d = folders(SETTINGSDIR)
            d = {};
            here = fileparts(which('EegPreset'));          % code/
            if ~isempty(here)
                d{end+1} = fullfile(fileparts(here), 'resources', EegPreset.FOLDER);
            end
            if nargin >= 1 && ~isempty(SETTINGSDIR)
                d{end+1} = fullfile(SETTINGSDIR, EegPreset.FOLDER);
            end
        end

        % ------------------------------------------------------------------
        % Every preset that can be found, shipped ones first, then the user's,
        % each sorted by 'order' and then by name. A preset in the user's
        % folder replaces a shipped one with the same name.
        function P = list(SETTINGSDIR)
            if nargin < 1, SETTINGSDIR = ''; end
            P = struct('name', {}, 'order', {}, 'file', {}, 'user', {});
            dirs = EegPreset.folders(SETTINGSDIR);
            for k = 1:numel(dirs)
                f = dir(fullfile(dirs{k}, ['*' EegPreset.EXT]));
                for i = 1:numel(f)
                    try
                        one = EegPreset.load(fullfile(dirs{k}, f(i).name));
                    catch E
                        warning('EegPreset:unreadable', '%s: %s', f(i).name, E.message);
                        continue
                    end
                    one.file = fullfile(dirs{k}, f(i).name);
                    one.user = k > 1;                      % the user's folder
                    old = find(strcmpi({P.name}, one.name), 1);
                    if isempty(old)
                        P(end+1) = EegPreset.row(one); %#ok<AGROW>
                    else
                        P(old) = EegPreset.row(one);       % the user's wins
                    end
                end
            end
            if isempty(P), return, end
            [~, ord] = sortrows([[P.order]', (1:numel(P))']);
            P = P(ord);
        end

        % ------------------------------------------------------------------
        function r = row(one)
            r = struct('name', one.name, 'order', one.order, 'file', one.file, 'user', one.user);
        end

        % ------------------------------------------------------------------
        % Read one preset file and fill in what it leaves out.
        function P = load(file)
            raw = jsondecode(fileread(file));
            P = struct();
            P.name     = EegPreset.getField(raw, 'name', '');
            if isempty(P.name)
                [~, P.name] = fileparts(file);
            end
            P.order    = EegPreset.getField(raw, 'order', 99);
            P.requires = EegPreset.getField(raw, 'requires', '');
            P.steps    = EegPreset.cellOf(EegPreset.getField(raw, 'steps', {}));
            P.settings = EegPreset.getField(raw, 'settings', struct());
            P.controls = EegPreset.getField(raw, 'controls', struct());
            P.code     = EegPreset.cellOf(EegPreset.getField(raw, 'code', {}));
            P.batch    = EegPreset.getField(raw, 'batch', struct());
            P.check    = EegPreset.cellOf(EegPreset.getField(raw, 'check', {}));
            P.notes    = EegPreset.cellOf(EegPreset.getField(raw, 'notes', {}));
            P.file     = file;
            P.user     = false;
        end

        % ------------------------------------------------------------------
        % Settings of a preset as an n x 3 cell array {button, parameter,
        % value}. JSON field names cannot hold spaces, so a button key like
        % 'excessive signal' is written with an underscore ('excessive_signal')
        % or a dash; both are turned back here, and the real key is whichever
        % the spreadsheet knows.
        function rows = settingRows(P, T)
            rows = cell(0, 3);
            if ~isstruct(P.settings), return, end
            buttons = fieldnames(P.settings);
            for b = 1:numel(buttons)
                key  = EegPreset.buttonKey(buttons{b}, T);
                pars = P.settings.(buttons{b});
                if ~isstruct(pars), continue, end
                names = fieldnames(pars);
                for p = 1:numel(names)
                    v = pars.(names{p});
                    if isnumeric(v) || islogical(v)
                        v = num2str(double(v));
                    end
                    rows(end+1, :) = {key, names{p}, char(string(v))}; %#ok<AGROW>
                end
            end
        end

        % ------------------------------------------------------------------
        % A button key as a JSON field name: 'excessive signal' cannot hold a
        % space, so it is written 'excessive_signal'.
        function f = fieldName(button)
            f = matlab.lang.makeValidName(regexprep(char(button), '\s+', '_'));
        end

        % ------------------------------------------------------------------
        % The button key as the spreadsheet spells it. The JSON side may use
        % underscores, dashes or camel case for the same button, so both sides
        % are compared with those removed and in lower case.
        function key = buttonKey(name, T)
            key = name;
            if nargin < 2 || isempty(T), return, end
            known = unique(cellstr(T.button));
            flat  = @(s) lower(regexprep(char(s), '[^A-Za-z0-9]', ''));
            hit = find(strcmp(cellfun(flat, known, 'UniformOutput', false), flat(name)), 1);
            if ~isempty(hit)
                key = known{hit};
            end
        end

        % ------------------------------------------------------------------
        % Write the current state as a preset in <SETTINGSDIR>/presets: every
        % setting, the chosen steps, the main-window dropdowns, the Execute
        % code box and the batch prefix. Returns the file written.
        function file = saveCurrent(app, name, labels, prefix, SETTINGSDIR)
            name = strtrim(name);
            if isempty(name)
                error('EegPreset:noName', 'Give the preset a name.');
            end
            d = fullfile(SETTINGSDIR, EegPreset.FOLDER);
            if ~exist(d, 'dir'), mkdir(d); end
            data = guidata(app.eeg_workflow);
            T = data.params;

            S = struct('name', name, 'order', 50);
            S.steps = labels(:)';
            % every setting, grouped per button
            S.settings = struct();
            for i = 1:height(T)
                b = char(T.button(i));
                f = EegPreset.fieldName(b);                % JSON-safe field
                S.settings.(f).(char(T.parname(i))) = char(T.value(i));
            end
            S.controls = struct();
            for c = {'ICAtypeDropDown', 'ICANumberDropDown', 'dropdownInterpolateToChanlocs'}
                try
                    if isprop(app, c{1}) && isvalid(app.(c{1}))
                        S.controls.(c{1}) = app.(c{1}).Value;
                    end
                catch
                end
            end
            try
                v = app.textareaExecuteCode.Value;
                if ~isempty(v) && ~all(cellfun(@isempty, cellstr(v)))
                    S.code = cellstr(v)';
                end
            catch
            end
            if nargin >= 4 && ~isempty(prefix)
                S.batch = struct('prefix', prefix);
            end
            S.notes = {sprintf('saved from the main window on %s', datestr(now, 'yyyy-mm-dd HH:MM'))}; %#ok<TNOW1,DATST>

            file = fullfile(d, [matlab.lang.makeValidName(name) EegPreset.EXT]);
            fid = fopen(file, 'w');
            if fid < 0
                error('EegPreset:cannotWrite', 'Cannot write %s.', file);
            end
            fprintf(fid, '%s', jsonencode(S, 'PrettyPrint', true));
            fclose(fid);
        end

        % ------------------------------------------------------------------
        % Apply a preset: check the plugin it needs, write its settings, set
        % the main-window controls, fill the Execute code box, and log what
        % happened. Returns the steps for the batch window and the lines to
        % show in the 'check the settings' popup.
        function [steps, check] = apply(app, P, say)
            hObject = app.eeg_workflow;
            data = guidata(hObject);
            if nargin < 3 || ~isa(say, 'function_handle')
                say = @(varargin) fprintf([varargin{1} '\n'], varargin{2:end});
            end

            if ~isempty(P.requires)
                EegPreset.requirePlugin(P.requires, P.name, say);
            end

            rows = EegPreset.settingRows(P, data.params);
            miss = {};
            for k = 1:size(rows, 1)
                if isempty(EegParams.rowIndex(data.params, rows{k,1}, rows{k,2}))
                    miss{end+1} = sprintf('%s.%s', rows{k,1}, rows{k,2}); %#ok<AGROW>
                    continue
                end
                data.params = EegParams.set(data.params, rows{k,1}, rows{k,2}, rows{k,3});
            end
            data.params = EegParams.validate(data.params);
            EegParams.save(data.params, data.SETTINGSDIR);
            EegParams.refresh(app, data.params);

            fn = fieldnames(P.controls);
            for k = 1:numel(fn)
                try
                    ctl = app.(fn{k});
                    if any(strcmp(ctl.Items, P.controls.(fn{k})))
                        ctl.Value = P.controls.(fn{k});
                    end
                catch
                end
            end
            if ~isempty(P.code)
                try
                    app.textareaExecuteCode.Value = P.code(:);
                catch
                    say(' *** warning *** could not fill the Execute code box');
                end
            end
            guidata(hObject, data);

            say('Pipeline preset %s: %d settings written, steps %s', P.name, ...
                size(rows, 1) - numel(miss), strjoin(P.steps, ' > '));
            for k = 1:numel(P.notes)
                say('- %s', P.notes{k});
            end
            if ~isempty(miss)
                say(' *** warning *** %d settings of this preset no longer exist and were skipped: %s', ...
                    numel(miss), strjoin(miss, ' '));
            end
            check = P.check;
            if isempty(check)
                check = {'Lookup: the montage of your system, and whether a reference channel is added'
                         'Line frequency: 50 Hz (Europe) or 60 Hz (Americas)'
                         'Open: how your files are read'};
            end
            say('- set yourself, by recording system: %s', strjoin(check, '; '));
            steps = P.steps;
        end

        % ------------------------------------------------------------------
        % An EEGLAB plugin a preset needs: fail at the click, not halfway
        % through a batch run, and put it on the path (its functions often
        % live in subfolders EEGLAB does not add).
        function requirePlugin(plugin, name, say)
            d = fileparts(which(['eegplugin_' lower(plugin)]));
            if isempty(d)
                d = fileparts(which(['eegplugin_' plugin]));
            end
            if isempty(d)
                error('EegPreset:missingPlugin', ...
                    ['The %s preset needs the %s plugin, which is not installed. ' ...
                     'Install it in EEGLAB (File > Manage EEGLAB extensions).'], name, plugin);
            end
            addpath(genpath(d));
            if isa(say, 'function_handle')
                say('- %s plugin found (%s) and added to the path', plugin, d);
            end
        end

        % ------------------------------------------------------------------
        function v = getField(s, name, dflt)
            v = dflt;
            if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
                v = s.(name);
            end
        end

        % ------------------------------------------------------------------
        % jsondecode gives a cell array for mixed strings and a string array
        % (or one char row) otherwise.
        function c = cellOf(v)
            if isempty(v)
                c = {};
            elseif iscell(v)
                c = cellfun(@char, v(:)', 'UniformOutput', false);
            else
                c = cellstr(v)';
            end
        end

    end
end
