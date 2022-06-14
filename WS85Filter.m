function omegaFiltered=WS85Filter(omega, D, phi, dt);

%code to calculate the omegaMean value based on Wright and Short 1985
%paper. For Narrabeen, WS85 found optimum fit with phi=5, D=30.
%inputs 
%omega = time series of dimensionless fall velocity
% dt = time step in Omega, seconds
%D = number of days used in back filtering
%phi = number of days when beach memory is 10%

dt = dt./3600./24;
D = round(D./dt);
phi = round(phi./dt);
omega=flipud(omega(:));
omegaFiltered=nan(size(omega));
for jj=1:length(omega)-D
    for ii=1:D
        T1(ii)=10^(-ii./phi);
        T2(ii) = omega(jj+ii).*T1(ii);
    end
    omegaFiltered(jj) = (sum(T1)).^(-1).*sum(T2);
end

omegaFiltered=flipud(omegaFiltered);

