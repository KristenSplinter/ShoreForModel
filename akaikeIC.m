function [AIC]=akaikeIC(xm,x,k)
% Akaike's Information Criteria determine's the relative appropriate of a
% model. See Akaike, 1974 or Kuriyama, 2012.

% xm - measured data
% x - model data
% k - the number of free parameters
n = length(xm);
resid = xm-x;
AIC = n.*(log10(2*pi)+1)+n.*log10(var(resid))+2.*k;