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
%%TLE data of pratham
% 1 99999U 16005A   16((epoach year)) 270.24769735(epoach day of year)  .00000368(1st derivative of mean motion)  
%00000-0(second derivative of mean motion)  70534-4(B-star drag) 0  1235 (checksum)
%2 99999  98.2056(inclination-degree) 327.8012(right ascension of ascending node) 0027951(eccentricity-decimal point assumed) 
% 267.5222(argument of perigee- deg)   9.4278(mean anamoly-deg) 14.63691895(mean motion)    1(revolution number at epoach)3(checksum)


% output of the sgp model here is in position (m) and velocity (m/s)
% in the test data it is in position (km) and velocity (km/s)
%% TLE DATA of pratham dated 19 July 2017 from n2yo
% 1 41783U 16059A   17198.94132223 +.00000076(1st time deri mean motion) +00000-0 +24125-4 0  9997
% 2 41783 098.1648(incl) 259.9939(rsc) 0034865(ecc) 065.5965(arg prg) 294.8869(mean anaomaly) 14.62844671043094
%%
%%get latest tle data from https://www.space-track.org/#/tle 
%%
% commented data belongs to sgp test data
% non commented is of pratham
% meanMo =    14.62844671;  %16.05824518;
% orbEcc=0.0034865;% 0.0027951;%0.0086731;
% orbInc=  098.1648;%72.8435;
% meanAno=  294.8869;% 9.4278;%110.5714;   %% change
% argPer =065.5965;%267.5222; %52.6988;     %% change
% rghtAsc= 259.9939;%327.8012;%115.9689;     %% change
% SGPdragp1 =0.00000076;%.00073094;  % from ---First time derivative of mean motion(rad/min^2)
% SGPdragp2 =0;%0.13844*10^(-3);     % Second time derivative of mean motion (rad/min^3
% SGP4dragp = 0.24125*10^(-4);%0.66816*10^(-4);
% %Epyear EpJD and EpTime used for calculation of t0
% % BUT   t0 is not used in sgp.m 
% %sgp is meant to run when tle data equal to epoach
% EpYear =  0; %1980
% EpJD   = 0; %275.98..
% EpTime  = 0 ;
% revNo = 9;%10;
% %%launch 2
% % %  n0 = 2*pi*meanMo/1440; % Mean motion (rad/min) 
% % %     e0 = 0.0045; % Eccentricity (0.0<=e0>=1.0)   2)
% % %     i0 = pi*orbInc/180; % Inclination (rad)      3)
% % %     M0 = 5.0563; % Mean anomaly (rad)    4)
% % %     w0 = 1.2269; % Argument of perigee (rad)   5)
% % %     Ohm0 = 4.8003;
% 
% 
% %%  
% 
%     n0 = 2*pi*meanMo/1440; % Mean motion (rad/min) 
%     e0 = orbEcc; % Eccentricity (0.0<=e0>=1.0)   2)
%     i0 = pi*orbInc/180; % Inclination (rad)      3)
%     M0 = pi*meanAno/180; % Mean anomaly (rad)    4)
%     w0 = pi*argPer/180; % Argument of perigee (rad)   5)
%     Ohm0 = pi*rghtAsc/180; % Right ascension of the ascending node (rad) 6)
%     dn0 = 2*2*pi*SGPdragp1/(1440^2); % First time derivative of mean motion(rad/min^2)
%     ddn0 = 6*2*pi*SGPdragp2/(1440^3); % Second time derivative of mean motion (rad/min^3)
%     Bstar = SGP4dragp; % SGP4 type drag coefficient
%     % Mean year (400 year period) in days:
%     meanYearDays = (400*365 + 4 * (100/4 - 1) + 1) / 400;
%     % Time of epoch (since y2k)
%     t0 = (EpYear*meanYearDays + EpJD + EpTime)*1440;
%     modTLE = [t0 dn0 ddn0 Bstar i0 Ohm0 e0 w0 M0 n0 revNo];
%     
%     %%
%%
%%Launch 2
orbInc=  098.2089;
meanMo =    14.62844671;
 SGPdragp1 =0.00000076;%.00073094;  % from ---First time derivative of mean motion(rad/min^2)
SGPdragp2 =0;%0.13844*10^(-3);     % Second time derivative of mean motion (rad/min^3
SGP4dragp = 0.24125*10^(-4);%0.66816*10^(-4); Bstar
n0 = 2*pi*meanMo/1440; % Mean motion (rad/min) 
    e0 = 0.0045; % Eccentricity (0.0<=e0>=1.0)   2)
    i0 = pi*orbInc/180; % Inclination (rad)      3)
    M0 = 5.0563; % Mean anomaly (rad)    4)
    w0 = 1.2269; % Argument of perigee (rad)   5)
    Ohm0 = 4.8003; %rsc
    dn0 = 2*2*pi*SGPdragp1/(1440^2); % First time derivative of mean motion(rad/min^2)
    ddn0 = 6*2*pi*SGPdragp2/(1440^3); % Second time derivative of mean motion (rad/min^3)
    Bstar = SGP4dragp; % SGP4 type drag coefficient
    t0 =0;
    revNo = 9;%10;
    %%trying to change just rsc node
    load('D:\Dropbox\Pratham\pratham propogator\SGP4 Validation\launch 1\modTLE_launch1.mat')
modTLE = modTLE_launch1;
modTLE(6) =4.8003 ;

%  modTLE = [t0 dn0 ddn0 Bstar i0 Ohm0 e0 w0 M0 n0 revNo];
modTLElaunch2 = modTLE;
save modTLElaunch2;
%%
dT = 0:1:4*60*60; % in seconds (24 hour)
[X, V] = sgp(modTLE, dT/60); 
SGP_test_case_launch2 = [dT/60; X'; V'];
% plot(dT, X)
save SGP_test_case_launch2.mat 
%% Summary and Conclusion
%Test case has be taken from SGP4 model from the referece 
%%http://www.celestrak.com/NORAD/documentation/spacetrk.pdf page 80
% gave into sgp.m and validated