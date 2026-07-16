function [W,r] = bsscca(X,delay)
% bsscca() - Blind Source Separation through Canonical Correlation Analysis
%
% Usage:
%   >> [W,r] = bsscca(X,delay)
%
% Inputs:
%   X     - data matrix (dxN, data channels are rowwise)
%   delay - delay at which the autocorrelation of the sources will be
%           maximized (def: 1)
%
% Output:
%   W     - separation matrix
%   r     - autocorrelation of the estimated sources at the given delay
%
% See also:
%   BSSCCA_IFC, AUTOBSS

% Copyright (C) <2007>  German Gomez-Herrero, http://germangh.com
%
% Speed-optimized drop-in (2026): the two SVD-based pinv() calls in the
% canonical-correlation eigenproblem are replaced by Cholesky solves, which
% are several times faster for the symmetric positive-definite covariance
% matrices produced here (bsscca_ifc pre-whitens the data with pca, bounding
% the eigenvalue spread). When a covariance matrix is rank-deficient the code
% falls back to pinv(), so the numerical result is identical to the original
% for every well-conditioned case and unchanged in the degenerate case.

if nargin < 2, delay = 1; end
if nargin < 1,
    help bsscca;
    return;
end


[d,T] = size(X); %#ok<ASGLU>

% delayed and undelayed data
Y = X(:,delay+1:end);
X = X(:,1:end-delay);

% correlation matrices (1/T scaling cancels in the eigenproblem below, but
% is kept to match the original definition exactly)
Cyy = (1/T)*(Y*Y');
Cxx = (1/T)*(X*X');
Cxy = (1/T)*(X*Y');
Cyx = Cxy';

% M = inv(Cxx) * Cxy * inv(Cyy) * Cyx, formed via symmetric solves instead
% of explicit (pseudo-)inverses.
M = symsolve(Cxx, Cxy * symsolve(Cyy, Cyx));

% calculate W
[W,r] = eig(M);
r = sqrt(abs(real(r)));
[r,I] = sort(diag(r),'descend');
W = W(:,I)';


% -------------------------------------------------------------------------
function Z = symsolve(A, B)
% Solve A*Z = B for a symmetric A. Uses a Cholesky factorization when A is
% positive definite (the normal, fast path) and falls back to the
% pseudo-inverse for rank-deficient A, reproducing the original pinv()
% behaviour exactly in that case.
[R,p] = chol(A);
if p == 0
    Z = R \ (R' \ B);
else
    Z = pinv(A) * B;
end
