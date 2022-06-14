function R2 = calcR2(beta,Ho,T);
%
% function to calculate the elevation of the 2% exceedence level for runup,
% R2, based on empircal parameterizations of Stockdon et al. 2006.
% Includes both wave-induced setup and swash.
% Inputs:
% beta - local beach slope. Can be defined as either the area over
% significant swash activity or shoreline to dune toe.
% Ho - offshore wave height
% Lo - offshore wave length = gT^2/2pi
%
% Kristen Splinter, 2009.
g = 9.81;
Lo = g.*T.^2./(2.*pi);
R2 = 1.1.*(0.35.*beta.*sqrt(Ho.*Lo) + ...
    sqrt(Ho.*Lo.*(0.563.*beta.^2 + 0.004))./2);