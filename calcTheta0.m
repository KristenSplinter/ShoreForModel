function theta0 = calcTheta10(theta,T,h1)

%calcualtes 2nd (inshore) angle based on snell's law.
% results in degrees
g = 9.81;
%h2 = 10;%.*ones(size(T)); %location where offshore wave was taken
k = dispsol2(h1, 1./T);
kh = k.*h1;
%k2 = dispsol2(h2, 1./T);
%kh2 = k2.*h2;
Co = 1./2.*g.*T./pi;
C1 = g.*T./(2.*pi).*tanh(kh);
theta0 = asind(sind(theta)./C1.*Co);