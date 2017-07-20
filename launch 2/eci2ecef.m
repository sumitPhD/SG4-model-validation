function TEI = eci2ecef(today,equinox, stperut,W_EARTH_ROT,t)
% This block gives DCM for converting from ECI to ECEF frames
% inputs: today's date, equinox, sidereal time per unit time, angular velocity of earth's rotation, simulation time
% output: DCM

% 23hr 56m 4.09s per day


% ut_sec = (today - equinox)*24*60*60 + t*60; % universal time vector in sec
ut_sec = (198.94132223-78.9375)*24*60*60 +t*60;  %ut_sec in sec
% delete 60*60 in final
st_sec = stperut*ut_sec;    % sidereal time vector in sec

phi = st_sec*W_EARTH_ROT;            % sidereal time vector in rad

TEI = [ cos(phi) sin(phi) 0;
       -sin(phi) cos(phi) 0;
        0        0        1];

