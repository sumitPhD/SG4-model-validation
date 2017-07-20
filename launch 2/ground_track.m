%%% This script is to determine ground track
%% 
% First orbital element is given to sgp by getorbitdata.m then position and velocity vector is calculated
% now here for each of the position value (output of sgp is in eci frame). Here eci2 ecef converted and then
% by using ecef2lla (matlab inbuilt function) lat long is calculated.
% sgp here is giving output in m and m/s

today = 0;
equinox = 0;
stperut = 1.00273790935;    % siderial time = stperut * universal time  
W_EARTH_ROT = 2*pi/(24*60*60);   % rotation angular velocity of earth, SI
T = SGP_test_case_launch2(1,:); % minutes
x = SGP_test_case_launch2(2:4,:);
N = length(x);
LLA_test_case_launch2 = zeros(4,N);


for i =1:N % keeping the gap high as the (it was taking too much time , not needed that much accuracy)
    TEI = eci2ecef(today,equinox, stperut,W_EARTH_ROT, T(i)+60); % 60*60 = 1hr for the second launch %%%%VVVI
    X_ECEF = (TEI*x(:,i))';
    LLA = ecef2lla(X_ECEF);
    LAT = LLA(1);  % in deg
    LONG = LLA(2); % in deg
    ALT = LLA(3)/1000;  % altitude in Km
    LLA_test_case_launch2(1,i) = T(i);
    LLA_test_case_launch2(2:4,i) = [LAT; LONG; ALT];
% plot(LLA_test_case(3,i),LLA_test_case(2,:))
% hold on
end
save LLA_test_case_launch2 
plot(LLA_test_case_launch2(3,:),LLA_test_case_launch2(2,:),'red')
xlabel('Longitude');
ylabel('Lattitude');
title('Ground track for pratham for 3 revolution');
%%% plotting the longitude on x axis and lattitude in y axis
