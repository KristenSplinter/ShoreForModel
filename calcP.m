function P = calcP(H,T,h)
%function to calculate wave power at an arbitrary depth
%P = ECn, where E = wave Energy, Cn=Cg = group velocity

g=9.81;
rho = 1025;
E = 1./16.*rho.*g.*H.^2;
if nargin==3
    k = dispsol2(h,1./T);
    kh=k.*h;
    C = sqrt(g./k.*tanh(kh));
    n = 1./2.*(1+kh.*((1-(tanh(kh)).^2)./tanh(kh)));

else
    %assume deep water
    C= 1./2.*g.*T./pi;
    n=1./2;
end
    P = E.*C.*n;
