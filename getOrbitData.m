%%SGP 4 model's test case given in (test case is run in 360 minutes(6hr) of interval)
%%http://www.celestrak.com/NORAD/documentation/spacetrk.pdf page 80
%1 88888U 80275.98708465 .00073094 13844-3 66816-4 0 8
%2 88888 72.8435 115.9689 0086731 52.6988 110.5714 16.05824518 105
%important terms
%for reference http://www.stltracker.com/resources/tle
% https://www.celestrak.com/NORAD/documentation/tle-fmt.asp
%also page 62  Chap orbital mechanics of book HandBook of Space Technology
%by Wilfried Ray et al
%1 88888U 80(epoach year)275.98708465(epoach day of year) .00073094(1st derivative of mean motion) 13844-3(second derivative of mean motion)
%66816-4(B-star drag) 0 8(checksum)
%2 88888 72.8435(inclination-degree) 115.9689(right ascension of ascending node) 0086731(eccentricity-decimal point assumed)
%52.6988(argument of perigee- deg) 110.5714(mean anamoly-deg)
%16.05824518(mean motion) 10(revolution number at epoach)5(checksum)%%

%%
% output of the sgp model here is in position (m) and velocity (m/s)
% in the test data it is in position (km) and velocity (km/s)
%%
meanMo = 16.05824518;
orbEcc= 0.0086731;
orbInc= 72.8435;
meanAno= 110.5714;
argPer = 52.6988;
rghtAsc= 115.9689;
SGPdragp1 =.00073094;
SGPdragp2 =0.13844*10^(-3);
SGP4dragp =0.66816*10^(-4);
%Epyear EpJD and EpTime used for calculation of t0
% BUT   t0 is not used in sgp.m 
%sgp is meant to run when tle data equal to epoach
EpYear =  0; %1980
EpJD   = 0; %275.98..
EpTime  = 0 ;
revNo = 10;


%%  
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
    
    %%
dT = 0:0.1:10000*60; % in seconds (360 minutes)
[X, V] = sgp(modTLE, dT/60); 
SGP_test_case = [dT/60; X'; V'];
% plot(dT, X)
save SGP_test_case.mat 
%% Summary and Conclusion
%Test case has be taken from SGP4 model from the referece 
%%http://www.celestrak.com/NORAD/documentation/spacetrk.pdf page 80
% gave into sgp.m and validated