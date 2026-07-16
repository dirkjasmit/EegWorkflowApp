% EEG_INTERP_ALL - leave-one-out spherical-spline interpolation of every channel.
%
% For each channel, reconstructs it from ALL other channels using the same
% spherical-spline math as EEG_INTERP, but computed in a single shared pass
% instead of nchan separate EEG_INTERP calls. This is the operation typically
% used for bad-channel detection (compare each channel to its interpolation).
%
% Usage:
%   >> [interpdata, W] = eeg_interp_all(EEG, method, params);
%
% Inputs:
%   EEG    - EEGLAB dataset (needs channel locations)
%   method - 'spherical' (default) or 'sphericalKang'. Ignored if params given.
%   params - [lambda m n] spherical-spline parameters (see EEG_INTERP).
%
% Outputs:
%   interpdata - [nchan x pnts x trials] data where row i is channel i
%                reconstructed from every other (located) channel. Channels
%                without locations are passed through unchanged.
%   W          - [nchan x nchan] leave-one-out interpolation operator, so that
%                interpdata(:,:) = W * EEG.data(:,:). W(i,i) = 0.
%
% Notes:
%   The per-channel result is numerically identical (to ~1e-13) to calling
%   EEG_INTERP(EEG, i, method) for each channel i and keeping channel i.
%
% See also: EEG_INTERP, POP_INTERP

% Author: drop-in companion to eeg_interp.m (2026)

function [interpdata, W] = eeg_interp_all(EEG, method, params)

if nargin < 1, help eeg_interp_all; return; end
if nargin < 2 || isempty(method), method = 'spherical'; end
if nargin < 3 || isempty(params)
    switch lower(method)
        case 'spherical',     params = [0    4  7];
        case 'sphericalkang', params = [1e-8 3 50];
        case 'sphericalcrd',  params = [1e-5 4 500];
        otherwise, error('(eeg_interp_all) unsupported method ''%s''', method);
    end
end
lambda = params(1);

% located channels only (mirrors eeg_interp's nonemptychans handling)
% -------------------------------------------------------------------
chanlocs = EEG.chanlocs;
located  = find(~cellfun('isempty', { chanlocs.theta }) & ~cellfun('isempty', { chanlocs.X }));
n = numel(located);
if n < 3
    error('(eeg_interp_all) need at least 3 located channels');
end

% unit-sphere electrode coordinates
% ---------------------------------
xe = [ chanlocs(located).X ];
ye = [ chanlocs(located).Y ];
ze = [ chanlocs(located).Z ];
rad = sqrt(xe.^2 + ye.^2 + ze.^2);
xe = xe./rad; ye = ye./rad; ze = ze./rad;

% full spline matrix, computed ONCE
% ---------------------------------
G = computeg_local(xe, ye, ze, xe, ye, ze, params);   % n x n

% build the leave-one-out operator W over located channels
% --------------------------------------------------------
Wloc = zeros(n, n);
allidx = 1:n;
for i = 1:n
    idx    = allidx(allidx ~= i);            % the n-1 other channels
    Mi     = [G(idx,idx) + lambda*eye(n-1); ones(1,n-1)];   % n x (n-1)
    u      = G(i,idx) * pinv(Mi);            % 1 x n  (Gsph_i * pinv(Mi))
    u      = u(1:end-1);                     % drop the constraint column
    % fold in eeg_interp's per-window mean-centering: allres = u*(v - mean) + mean
    Wloc(i, idx) = u + (1 - sum(u))/(n-1);
end

% assemble full-size operator (pass through non-located channels) and apply
% -------------------------------------------------------------------------
nchan = EEG.nbchan;
W = zeros(nchan, nchan);
for a = 1:nchan, W(a,a) = 1; end            % default: identity (pass-through)
W(located, :) = 0;
W(located, located) = Wloc;

data2d     = double(reshape(EEG.data, nchan, EEG.pnts*EEG.trials));
interpdata = W * data2d;
interpdata = reshape(interpdata, nchan, EEG.pnts, EEG.trials);

% ------------------------------------------------------------------------
% G function (order-0 Legendre via three-term recurrence; matches eeg_interp)
% ------------------------------------------------------------------------
function g = computeg_local(x,y,z,xelec,yelec,zelec, params)
x = x(:); y = y(:); z = z(:);
xelec = xelec(:).'; yelec = yelec(:).'; zelec = zelec(:).';
EI = 1 - sqrt( (x - xelec).^2 + (y - yelec).^2 + (z - zelec).^2 );
m = params(2); maxn = params(3);
g = zeros(size(EI)); Pnm1 = ones(size(EI)); Pn = EI;
for nn = 1:maxn
    g = g + ((2*nn+1)/(nn^m*(nn+1)^m)) * Pn;
    Pnp1 = ((2*nn+1).*EI.*Pn - nn.*Pnm1) / (nn+1);
    Pnm1 = Pn; Pn = Pnp1;
end
g = g/(4*pi);
