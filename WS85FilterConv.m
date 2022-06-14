function [omegaFiltered]=WS85FilterConv(omega, D, phi, dt);

% code to calculate the omegaMean value based on Wright and Short 1985
% paper.
% inputs 
% ~~~~~~
% omega = time series of dimensionless fall velocity
% dt = time step in Omega, seconds
% D = number of days used in back filtering - set this to 2*phi
% phi = number of days when beach memory is 10%
% The method utilises convolution theorem to apply the filter. It is much
% faster than using loops! However, the methodology does not give good
% results for the last phi data points which must be calculated using the
% slow looping method
% Mark davidson 11/7/12

disp('Computing equilibrium omega value ....WS85 convolution')

dt = dt./(3600.*24);
D = round(D./dt);
phi = round(phi./dt);
N=length(omega);
meanOmega=mean(omega);
omega=omega-meanOmega;

% Define filter back-to-front for convolution
ii=0:(D-1);
padding=zeros(D-1,1);
filterCoeff=10.^(-abs(ii)/phi);
filterCoeff=[padding(:);filterCoeff(:)];
window=filterCoeff/sum(filterCoeff);
Nw=length(window);

% perform convolution
omegaFiltered = conv(omega,window,'same');


% Finally add on mean
omegaFiltered=omegaFiltered+meanOmega;