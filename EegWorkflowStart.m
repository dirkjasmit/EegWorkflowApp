function app = EegWorkflowStart()
% EEGWORKFLOWSTART  Start the EEG Workflow app.
%
%   EegWorkflowStart          starts the app
%   app = EegWorkflowStart    also returns the app object
%
% Adds the folders of this repository to the MATLAB path and opens
% EegWorkflow.mlapp. Run this instead of opening the .mlapp directly, so that
% the classes in code/ and the functions, GUIs and resources in the other
% folders are found wherever the repository is installed. EEGLAB itself must
% already be on the path (or the app asks for its folder on first use).

here = fileparts(mfilename('fullpath'));
subs = {'code', 'functions', 'gui', 'external', 'resources'};
for k = 1:numel(subs)
    d = fullfile(here, subs{k});
    if isfolder(d) && ~contains([path pathsep], [d pathsep])
        addpath(d);
    end
end

h = EegWorkflow;
if nargout > 0
    app = h;
end
end
