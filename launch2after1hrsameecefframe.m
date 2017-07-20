%at t = 0 with given pratham tle data, below is the launch position 
% and velocity vector
xp_eci_launch1 =SGP_test_case(2:4,1);
xv_eci_launch1 = SGP_test_case(5:7,1);

stperut = 1.00273790935;    % siderial time = stperut * universal time  
W_EARTH_ROT = 2*pi/(24*60*60);   % rotation angular velocity of earth, SI
%  1.0e+06 *
% 
%    -1.2253
%    -6.9448
%     0.0000
%    -0.0011
%     0.0002
%     0.0074

% convert this into eci frame after 1 hr

%  ut_sec = (today - equinox)*24*60*60 + t*60; % universal time vector in sec
% ut_sec = (198.94132223-78.9375)*24*60*60 +t*60;  %ut_sec in sec
ut_sec = (198.94132223-78.9375)*24*60*60 ;

% delete 60*60 in final
st_sec = stperut*ut_sec;    % sidereal time vector in sec

phi = st_sec*W_EARTH_ROT;            % sidereal time vector in rad

TEI = [ cos(phi) sin(phi) 0;
       -sin(phi) cos(phi) 0;
        0        0        1]; % 
    xp_ecef_launch1 = TEI*xp_eci_launch1;
    xv_ecef_launch1 = TEI*xv_eci_launch1;
    
    
    ut_sec = (198.94132223-78.9375)*24*60*60 +60*60; %(launching after 1 hr)

% delete 60*60 in final
st_sec = stperut*ut_sec;    % sidereal time vector in sec

phi = st_sec*W_EARTH_ROT;            % sidereal time vector in rad
%ecef to eci conversion 
TIE = [ cos(phi) -sin(phi) 0;
       sin(phi) cos(phi) 0;
        0        0        1]; % 
    %%% ecef for launch 1 and 2 are same
    
    xp_eci_launch2 = TIE*xp_ecef_launch1;
    xv_eci_launch2 = TIE*xv_ecef_launch1;
[a,e,i,O,o,nu] = rv2orb(xp_eci_launch2,xv_eci_launch2,6.673e-11*5.9742e24);
orbit_ele_launch2 = [a,e,i,O,o,nu] ;
    save orbit_ele_launch2;
    
