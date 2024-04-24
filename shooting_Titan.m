% Script #1
% Gas explosion model

% This script automates the process of finding the gas/solid mass fraction required to calculate gas amount. 

% Steps for running:
% 1) Input cap thickness r21 in meters.
% 2) Input the launch distance in meters.
% 3) Input guess for gas/solid mass fraction (N) for each gas pressure. Ideas for initial gas/solid mass fractions (N) guesses can be found in "Full model outputs.xlsx". 
% 4) Update inputs and constants in 'func_explode' (Script #2). Don't forget this step!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:

%cap thickness (m)
r21=300;

% desired distance (launch distance or radar-bright halo extent)
distance_desired=4100;

% Guess gas/solid mass fraction (N) for each gas pressure. 
N_=[.1 .1 .1 .1 .1 .1 .1 .1 .1]; 

% Input gas pressure (Pa)
Pgz_=[0.5e6 1e6 2.5e6 5e6 7.5e6 10e6 15e6 20e6 30e6];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%timer
tic
%wait time for convergence
wait=10; 
%gain for proportion control
gain=1e-6; %*note: for 30 the gain of 1e-4 works. for 300 a smaller gain is better
%tolerance for final distance error
error_tolerance=10; % get within 10 m of desired distance
output=[];
for m=1:length(Pgz_)
    N=N_(m);
    Pgz=Pgz_(m);
while (1==1)
    %solve my explosion function
    dist=func_explode(N,r21,Pgz); %the max distance my function spits out
    %this is the distance
    distance=dist(3);
    error=distance-distance_desired;
    if (toc >= wait || N < 0) %then it took too long
        disp('could not converge in time');
        break; %end loop
    elseif (abs(error) <= error_tolerance)
        disp(['N: ', num2str(N),'.']);
        break; %end loop
    end
    N=N-error*gain;
end
output=[output; dist];
end
