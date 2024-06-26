function u=Uwave_ruessink2012(phs,Aw,Sw,Uw)
%
% u=Uwave_ruessink2012(phs,Aw,Sw,Uw)
%
% Intra-wave velocity time series parameterization, following Ruessink et
% al. (2012), "On the parameterization of the free stream non-linear wave
% orbital motion in nearshore morphodynamic models".
%
% INPUTS:
%
% phs : wave phase vs. time at which to evaluate u, as an array with limits [0,2*pi]
% Aw,Sw,Uw : wave shape parameters, get these from Uwave_ruessink2012_params()
%

% convert to B and psi
B = sqrt(Aw.^2 + Sw.^2);
psi = atan(Aw./Sw);

% % TEST: override with Doering & Bowen (1995) fit, which predicts larger
% % asymmetry and skewness
% B = max(0, 0.8 + 0.62*log10(Ur));  % confirmed is log10() from DB95 fig 5
% psi = -pi/2 + pi/2*tanh(0.73./Ur);

% phi parameter, eqn (12)
phi = -psi-pi/2;

% % v1: r parameter, direct inversion
% rguess=0.5;
% b=@(r)r./(1+sqrt(1-r.^2));
% rfun=@(r)B-3*b(r)./sqrt(2*(1-b(r).^2));
% nx=length(Ur);
% opt=optimset('Display','off');
% if(nx>1)
%   r=fsolve(rfun,rguess*ones(nx,1),opt);
% else
%   r=fzero(rfun,rguess,opt);
% end

% v2: r parameter, 2-step inversion.  First solve for b as a function of B,
% then solve for r as a function of b.  Doing this in 2 steps makes it
% easier to write TL code, though the old one-step version (v1, see above) may
% be slightly more efficient.  Results are identical, tested.
nx=length(Uw);
bguess=0.5;
rguess=0.5;
Bfun=@(b)3*b./sqrt(2*(1-b.^2));
bfun=@(r)r./(1+sqrt(1-r.^2));
opt=optimset('Display','off');
if(nx>1)
  b=fsolve(@(b)B-Bfun(b),bguess*ones(nx,1),opt);
  r=fsolve(@(r)b-bfun(r),rguess*ones(nx,1),opt);
else
  b=fzero(@(b)B-Bfun(b),bguess,opt);
  r=fzero(@(r)b-bfun(r),rguess,opt);
end

% Abreu et al. (2010) time series formula, eqn (4)
f = sqrt(1-r.^2);
f1 = sin(phs) + r.*sin(phi)./(1+sqrt(1-r.^2));
f2 = 1 - r.*cos(phs+phi);
u = Uw.*f.*f1./f2;
