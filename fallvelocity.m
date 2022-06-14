function [w]=fallvelocity(D,Tw)
% D = Grain size (m)
% Tw = Temp in degrees C
% w returned in m/s

D=D*100;

ROWs=2.75;	% Density of sand (Mark 2.75, Kristen, 2.65?)
g=981;		% Gravity cm/s^2

T   =[5 10 15 20 25];
v   =[0.0157 0.0135 0.0119 0.0105 0.0095];
ROW =[1.028 1.027 1.026 1.025 1.024]; 

vw=interp1(T,v,Tw);
ROWw=interp1(T,ROW,Tw);

A = ((ROWs-ROWw)*g*(D.^3))./(ROWw*(vw.^2));

if A < 39
      w=((ROWs-ROWw)*g*(D.^2))./(18*ROWw*vw);
else
	if A <10.^4   
   	  w=((((ROWs-ROWw)*g./ROWw).^0.7)*(D.^1.1))./(6*(vw.^0.4));
	else
      w=sqrt(((ROWs-ROWw)*g*D)./(0.91*ROWw));
   end
end

w=w./100; % convert to SI (m/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%