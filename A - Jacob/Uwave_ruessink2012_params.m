function [Aw,Sw,Uw]=Uwave_ruessink2012_params(Hmo,k,omega,h)
%
% [Aw,Sw,Uw]=Uwave_ruessink2012_params(Hmo,k,omega,h)
%
% Helper function for Uwave_ruessink2012().  Calculates asymmetry (Aw),
% skewness (Sw), and amplitude (Uw).  This function can also be useful if
% you just want the parameters and don't want u(t).
%
% Hmo   : significant wave height, m
% k     : wavenumber, rad/m
% omega : wave frequency, rad/s
% h     : water depth, m
% 

% Ursell number, eqn (6)
aw=Hmo/2;
Ur=3/4*aw.*k./(k.*h).^3;

% wave velocity magnitude, linear theory, stated in text
Hrms=Hmo/1.4;
Uw = omega/2.*Hrms./sinh(k.*h);

% non-linearity param B, eqn (9).  Using parameter values quoted in text
p1=0;
p2= .857;  % +/- .016
p3=-.471;  % +/- .025
p4= .297;  % +/- .021
ee=(p3-log10(Ur))./p4;  % dec 17, 2021: This should use log10(), confirmed vs fig 2
dens=1+exp(ee);
B0 = p1 + (p2-p1)./dens;

% non-linearity phase param psi, eqn (10)
p5=.815;  % +/- .055
p6=.672;  % +/- .073
psi0 = -pi/2 + pi/2*tanh(p5./Ur.^p6);

% The above baseline B0 and psi0 give a baseline value for asymmetry and
% skewness (Aw and Sw)
Aw = B0.*sin(psi0);
Sw = B0.*cos(psi0);

% apply masking to points with zero wave height, avoids NaN
imask=find(Hmo==0);
Aw(imask)=0;
Sw(imask)=0;
Uw(imask)=0;
