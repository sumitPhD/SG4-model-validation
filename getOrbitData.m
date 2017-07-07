    n0 = 2*pi*meanMo/1440; % Mean motion (rad/min)
    e0 = orbEcc; % Eccentricity (0.0<=e0>=1.0)
    i0 = pi*orbInc/180; % Inclination (rad)
    M0 = pi*meanAno/180; % Mean anomaly (rad)
    w0 = pi*argPer/180; % Argument of perigee (rad)
    Ohm0 = pi*rghtAsc/180; % Right ascension of the ascending node (rad)
    dn0 = 2*2*pi*SGPdragp1/(1440^2); % First time derivative of mean motion(rad/min^2)
    ddn0 = 6*2*pi*SGPdragp2/(1440^3); % Second time derivative of mean motion (rad/min^3)
    Bstar = SGP4dragp; % SGP4 type drag coefficient
    % Mean year (400 year period) in days:
    meanYearDays = (400*365 + 4 * (100/4 - 1) + 1) / 400;
    % Time of epoch (since y2k)
    t0 = (EpYear*meanYearDays + EpJD + EpTime)*1440;
    modTLE = [t0 dn0 ddn0 Bstar i0 Ohm0 e0 w0 M0 n0 revNo];