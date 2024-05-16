%% ripple model equations
% required defined factors: 
    % ripple length (lambda)
    % ripple height (eta)
    % time-averaged absolute shields stress (named shields_aa here)
        % defined in category C1
    % d50, d90: grain size (grain size, static values assumed)
    % mobility number psi 
%% defined factors for ripple dimensions
psimax = (max(u))^2 / ((s-1) * g * d50); % maximum mobility number, 
% used to calculate multipliers nlambda and neta in lambda and eta eqns 

% to get values of multipliers mlamda and meta necessary for lambda/eta:
    % lambda multipliers
if d50 <= 0.22,
    mlambda = 0.73
elseif d50 <= 0.3,
    mlambda = .73 + ((0.27*(d50-.22))/(0.3-0.22));
else 
    mlambda = 1;
end 

    % eta multipliers
if d50 <= 0.22,
    meta = .55;
elseif d50 <= 0.3,
    meta = 0.55 + ((0.45 * (d50 - 0.22))/ (0.3-0.22));
else 
    meta = 1;
end 

% to get values of multipliers nlambda and neta, which are important for
% smooth transition between flat bed and ripple regimes:
    % both are going to be called "n" because they are equal (eqn B.5) just
    % used in their respective equations for ripple dimension 

if psimax <= 190,
    n = 1;
elseif psimax <= 240,
    n = 0.5 * (1 + cos(pi * ((psimax - 190)/(240-190))));
else 
    n = 0;
end 
%% actual ripple length and height equations (B.1 and B.2)

% eta (ripple height) equation, 'ahat' given in eqn 9 of category A 
eta = ahat * meta * n * (0.275 - ((0.022 * psimax)^0.42));
% lambda (ripple length) equation
lambda = ahat * mlamda * n * (1.97 - ((0.44 * psimax)^0.21));

%% factor mu (VDA13 eqn A.2)
% mu = fine sand adjustment for sheet flow
% without factor adjustment, underestimation of net transport rate

if d50 <= 0.15, % d50 static input value entered in millimeters
    mu = 6; 
elseif d50 <= 0.2, % 
    mu = 6 - ((5 * (d50 - 0.15))/(0.2 - 0.15)); 
else 
    mu = 1;
end 

%% current-related bed roughness (eqn A.1)
% ksdelta = equation for current-related bed roughness from Ribberink 1998
% requires d90 value 
rough = d50 * (mu + 6(shields_aa - 1)); % expression in the "max" fxn 
ksdelta = max((3*d90), rough) + ((0.4 * eta^2)/lambda);

%% wave-related bed roughness (eqn. A.5)

ksw = max(d50, rough) + ((0.4 * eta^2)/lambda);


