function p = pnorm(z,onesided)
% calculate the cumulative normal distribution for z values
% z: vector or matrix of z-values
% onesided: true/false for one or twosided probabilities (def false so
%           two-sided testing)

if nargin<2
    onesided = false;
end

if onesided
    p = normcdf(z,0,1);
else
	p = normcdf(-abs(z),0,1)*2;
end



