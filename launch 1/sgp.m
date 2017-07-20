function [R,dR] = sgp(modTLE,dT)
%This code is a matlab implementation of the SGP Fortran code in the
% Startrack report #3.
%
% The matlab code was derived mainly by porting and editing the
% FORTRAN routines in NORAD’s Spacetrack report #3.
% According to a statement in that report, the document is free
% of copyrights and open to unlimited public distribution.
%
% This matlab code is distributed under the same conditions
%
%
% Created 15/10-2001 by
    % Klaus Krogsgaard & Torsten Lorentzen
% - - - - - - - - - - - - - - - - - - - - - - - - -
%
% SGP model for predicting orbit position.
%

% - - - - - - - - - - - - - - - - - - - - - - - - -
% Constants given by the modified NORAD TLE (all at epoch)
% - - - - - - - - - - - - - - - - - - - - - - - - -


t0 = modTLE(1);
dn0 = modTLE(2); % first time derivative of the mean motion
ddn0 = modTLE(3); % second time derivative of mean motion
Bstar = modTLE(4); % Drag parameter
i0 = modTLE(5); % inclination
Ohm0 = modTLE(6); % right ascension of the ascending node
e0 = modTLE(7); % eccentricity
w0 = modTLE(8); % argument of perigee
M0 = modTLE(9); % mean anormaly
n0 = modTLE(10); % mean motion
revNo = modTLE(11); % Number of revolutions before epoch
% - - - - - - - - - - - - - - - - - - - - - - - - -
% Other constants
% - - - - - - - - - - - - - - - - - - - - - - - - -
GM = 398603 * 10^9;
er = 6378.135;
ke = 0.0743669161; %converted to er/min ^3/2
% aE = 6378160; % equatorial radius of the Earth
aE = 1;% give result in Earth radii
min_per_day = 1440;
sec_per_day = 86400;
km_per_er = er;
J2 = 5.413080*10^-4 * 2; % second gravitational zonal harmonic of the Earth
J3 = -0.253881*10^-5; % third gravitational zonal harmonic of the Earth
J4 = -0.62098875*10^-6 * 8/3; % fourth gravitational zonal harmonic of the Earth
% - - - - - - - - - - - - - - - - - - - - - - - - -
% Local constants
% - - - - - - - - - - - - - - - - - - - - - - - - -
a1 = (ke/n0)^(2/3);
delta1 = 3/4 * J2 * (aE/a1)^2 * (3*((cos(i0))^2)-1)/((1-e0^2)^(3/2));
a0 = a1 * (1 - 1/3*delta1 - delta1^2 - 134/81*delta1^3);
p0 = a0 * (1 - e0^2);
q0 = a0 * (1-e0);
L0 = M0 + w0 + Ohm0;
dOhm = - (3/2) * J2 * (aE/p0)^2 * n0 * cos(i0);
dw = (3/4) * J2 * (aE/p0)^2 * n0 * (5*(cos(i0))^2-1);
% - - - - - - - - - - - - - - - - - - - - - - - - -
% Secular effects of drag and gravitation
% - - - - - - - - - - - - - - - - - - - - - - - - -
a = a0 * ( n0 ./ (n0 + 2*(dn0/2)*(dT) + 3*(ddn0/6)*dT.^2)).^(2/3); % vector
e = zeros(1,length(dT)); % vector
for i=1:length(a)
if ( a(i)>q0 )
e(i) = 1 - q0/a(i);
else
e(i) = 10^-6;
end;
end;
p = a .* (1-e.^2); % vector
OhmS0 = Ohm0 + dOhm * dT; % vector
wS0 = w0 + dw * (dT); % vector
Ls = mod(L0 + (n0 + dw + dOhm)*(dT) + dn0/2 * (dT).^2 + ddn0/6 * (dT).^3,2*pi); % vector
% - - - - - - - - - - - - - - - - - - - - - - - - -
% Long period periodic effects
% - - - - - - - - - - - - - - - - - - - - - - - - -
axNSL = e .* cos(wS0); % vector
ayNSL = e .* sin(wS0) - 1/2 * J3/J2 * aE./p * sin(i0); % vector
L = mod(Ls - 1/4 * J3/J2 * aE./p .* axNSL * sin(i0) * (3+5*cos(i0))/(1+cos(i0)),2*pi); % vector
% - - - - - - - - - - - - - - - - - - - - - - - - -
% Iteration for short period periodics below...
% - - - - - - - - - - - - - - - - - - - - - - - - -
tol = 10^-12;
U = mod(L - OhmS0,2*pi); % vector
Ew1 = U; % vector
Ew2 = Ew1; % vector
dEw = (U - ayNSL.*cos(Ew1) + axNSL.*sin(Ew1) - Ew1) ./ (-ayNSL.*sin(Ew1) - axNSL.*cos(Ew1) + 1); % vector
for i=1:length(dEw)
if abs(dEw(i)) > 1
Ew2(i) = Ew1(i) + sign(dEw(i)); %vector
else
Ew2(i) = Ew1(i) + dEw(i);
end;
end;
for i=1:length(Ew1)
while abs(dEw(i))>tol
Ew1(i) = Ew2(i);
dEw(i) = (U(i) - ayNSL(i).*cos(Ew1(i)) + axNSL(i).*sin(Ew1(i)) - Ew1(i)) ./ ...
(-ayNSL(i).*sin(Ew1(i)) - axNSL(i).*cos(Ew1(i)) + 1);
if abs(dEw(i)) > 1
Ew2(i) = Ew1(i) + sign(dEw(i)); %vector
else
Ew2(i) = Ew1(i) + dEw(i);
end;
end;
end;
ecosE = axNSL.*cos(Ew2) + ayNSL.*sin(Ew2); % vector
esinE = axNSL.*sin(Ew2) - ayNSL.*cos(Ew2); % vector
SQeL = axNSL.^2 + ayNSL.^2; % vector
pL = a.*(1 - SQeL); % vector
r = a.*(1 - ecosE); % vector
dr = ke * sqrt(a)./r .* esinE; % vector
rdv = ke * sqrt(pL)./r; % vector
sinu = a./r .* (sin(Ew2) - ayNSL - axNSL.*esinE./(1+sqrt(1-SQeL))); % vector
cosu = a./r .* (cos(Ew2) - axNSL + ayNSL.*esinE./(1+sqrt(1-SQeL))); % vectro
u = mod(atan2( sinu,cosu ),2*pi); % vector
% - - - - - - - - - - - - - - - - - - - - - - - - -
% Short term perturbations
% - - - - - - - - - - - - - - - - - - - - - - - - -
rk = r + 1/4 * J2 * aE^2./pL * (sin(i0))^2 .* cos(2*u); % vector
uk = u - 1/8 * J2 * (aE./pL).^2 * (7 * (cos(i0))^2 - 1) .* sin(2*u); % vector
Ohmk = OhmS0 + 3/4 * J2 * (aE./pL).^2 * cos(i0) .* sin(2*u); % vector
ik = i0 + 3/4 * J2 * (aE./pL).^2 * sin(i0) * cos(i0) .* cos(2*u); % vector
% - - - - - - - - - - - - - - - - - - - - - - - - -
% Unit orientation vectors
% - - - - - - - - - - - - - - - - - - - - - - - - -
R = zeros(length(uk),3); % vector
dR = zeros(length(uk),3); % vector
for i=1:length(uk)
    M_vec = [ -sin(Ohmk(i))*cos(ik(i)) cos(Ohmk(i))*cos(ik(i)) sin(ik(i)) ]; % vector
N_vec = [ cos(Ohmk(i)) sin(Ohmk(i)) 0 ]; % vector
U_vec = M_vec * sin(uk(i)) + N_vec * cos(uk(i)); % vector
V_vec = M_vec * cos(uk(i)) - N_vec * sin(uk(i)); % vector
% - - - - - - - - - - - - - - - - - - - - - - - - -
% Position and velocity
% - - - - - - - - - - - - - - - - - - - - - - - - -
R(i,:) = rk(i) * U_vec; % vector
dR(i,:) = dr(i) * U_vec + rdv(i) * V_vec; % vector
end;
% - - - - - - - - - - - - - - - - - - - - - - - - -
% Transforming position to Cartesian coordinates in meters and velocity to meters/sec
% - - - - - - - - - - - - - - - - - - - - - - - - -
R = R * km_per_er /aE;
R = R * 1000;
dR = dR * km_per_er / aE * min_per_day / sec_per_day;
dR = dR * 1000;