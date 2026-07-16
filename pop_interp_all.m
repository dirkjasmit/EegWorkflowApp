% POP_INTERP_ALL - replace every channel with its leave-one-out spherical-
%                  spline interpolation (each channel reconstructed from all
%                  other channels), computed in one shared pass.
%
% This is the same per-channel result as looping POP_INTERP over every
% channel, but far faster (the spline matrix and setup are shared instead of
% recomputed nchan times). The returned dataset has its data REPLACED by the
% interpolated version.
%
% Usage:
%   >> [EEG, com] = pop_interp_all(EEG);            % pops up a method dialog
%   >> [EEG, com] = pop_interp_all(EEG, method);
%
% Inputs:
%   EEG    - EEGLAB dataset (must contain channel locations)
%   method - 'spherical' (default) or 'sphericalKang'
%
% Outputs:
%   EEG    - dataset whose data is replaced by the leave-one-out interpolation
%            (channels without locations are left unchanged)
%   com    - command string for the EEGLAB history
%
% See also: EEG_INTERP_ALL, EEG_INTERP, POP_INTERP

% Author: drop-in companion to pop_interp.m (2026)

function [EEG, com] = pop_interp_all(EEG, method)

com = '';
if nargin < 1, help pop_interp_all; return; end
if isempty(EEG.data), disp('(pop_interp_all) cannot process an empty dataset'); return; end
if isempty(EEG.chanlocs) || isempty([EEG.chanlocs.X])
    error('(pop_interp_all) requires channel locations');
end

% method (GUI only when not supplied)
% -----------------------------------
if nargin < 2
    uilist = { ...
        { 'style' 'text'      'string' 'Interpolation method:' } ...
        { 'style' 'popupmenu' 'string' 'Spherical|spherical (Kang et al.)' 'tag' 'method' } };
    geom = { [2 1] };
    [results, ~, ~, restag] = inputgui('uilist', uilist, 'geometry', geom, ...
        'title', 'Interpolate all channels (leave-one-out) -- pop_interp_all()', ...
        'helpcom', 'pophelp(''pop_interp_all'')');
    if isempty(results), return; end
    if restag.method == 2, method = 'sphericalKang'; else method = 'spherical'; end
elseif isempty(method)
    method = 'spherical';
end

% leave-one-out interpolation of every channel, and replace the data
% ------------------------------------------------------------------
% eeg_interp_all computes in double; eeg_checkset then applies EEGLAB's data
% precision option (single by default), exactly as a per-channel POP_INTERP /
% EEG_INTERP loop would. Result matches such a loop to the stored precision.
EEG.data = eeg_interp_all(EEG, method);
EEG      = eeg_checkset(EEG);

% history string
% --------------
com = sprintf('EEG = pop_interp_all(EEG, ''%s'');', method);

return;
