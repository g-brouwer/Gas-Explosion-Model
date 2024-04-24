function outy = func_explode(N,r21,Pgz)
% Script #2
% Gas explosion model

% Steps for running:

% 1) Input your block properties, gas properties, and time step.
% 2) Run shooting_Titan (Script #1).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input ejetced block properties:

%radius of spherical projectile (m)
R=1;

%projectile density (kg/m3)
rhos=917;

%drag coefficient
p.Cd=0.9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input gas properties:

% gas temperature in gas reservoir (K)
Tg=128;

%Specific gas constant
SGCG=518.3;  %518.3; %Specific gas constant of methane [J/kg K] %297 for N2

% Gas density
rhogz=Pgz/(Tg*SGCG); %initial gas density. 

% Ratio of specific heats (methane: 1.3, nitrogen: 1.4)
gamma=1.3; %1.4; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
% Input time step for the finite difference solution.

dt=.1; %0.001 (for 30). 0.1 works for 300. Its a matter of dt compared to t0 (total time). You have to play with this. Different values work best for different cap
% thicknesses.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Titan atmosphere properties:
tempz=94; %atmospheric temperature (K)
p.rhoa=5.3; %Titan atmosphere density (kg/m3) (Lorenz 2007)
p.eta=6.25e-6; %kinematic viscosity of the atmosphere (Yu et al., 2017)
Pa=0.145e6; %Titan's atmospheric pressure
rhoa=5.3; %air density Titan

% Other constants
g=1.35; % gravity
theta=pi/4; %launch angle 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate gas chamber height

%Calculate roots of polynomial with N & r21
r1=gas_chamber(r21, rhogz, rhos, N); 
r2=r21+r1; %gas chamber + caprock thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite difference solution of gas expansion phase

%storage
time=zeros();
u=zeros();
%initialize
r(1)=r2;
time(1)=0;
up=0;
u(1)=0;
xi(1)=0;
yi(1)=0;

    for i=1:1e8
        time(i+1)=time(i)+dt;        
    if (u >= 0)
        rdd(i)=(3*r(i)^2)*(Pgz*(r1/r(i))^(3*gamma)-Pa)/(rhos*(r2^3-r1^3)+rhoa*((r(i)+r21)^3-r2^3));
        du(i)=dt*rdd(i);
        u(i+1)=u(i)+du(i);
        ubar(i+1)=0.5*(u(i+1)+up);
        dr(i)=u(i)*dt;

        r(i+1)=r(i)+dr(i); 
        xi(i+1)=u(i+1)*cos(theta)*time(i+1);
        yi(i+1)=u(i+1)*sin(theta)*time(i+1);
        else
            break
        end
    end

umaxfind=find(u == max(u)); %max velocity of expanding envelope
umax=u(umaxfind);
t0=time(umaxfind);
p.t0=t0;

xi0=xi(umaxfind);
yi0=yi(umaxfind);
p.xi0=xi0;
p.yi0=yi0;
time(end)=[]; %remove last time value

%when u=0, that is the t' used to calculate tau
u_0=find(u < 0);
t_u_0=time(u_0);
p.tau=t_u_0-p.t0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now the projectile is launched into the moving gas. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the equation of motion for a projectile with drag.

%initial projectile velocity
p.v0=umax;
%intial projectile position
p.r0=sqrt(xi0^2+yi0^2); 
p.g=g;
%projectile density
p.rhos=rhos; 
%projectile surface area
p.A=4*pi*R^2;
%projectile mass
p.m=(4/3)*pi*R^3*p.rhos; 
%diameter (for drag function)
p.d=2*R;   
%45 degree envelope of expanding gas envelope
p.theta=theta; 

%initial conditions and timespan
x0=[xi0 yi0];
v0=p.v0*[cos(p.theta) sin(p.theta)];
%time interval that ode45 will integrate the system of ODE's (160/1600
%seconds)
% each row in the solution array corresponds to a value returned in column vector t
tspan=[0 1600];
                
%ODE45 accuracy settings
opts.reltol=1e-6;
opts.abstol=1e-6;

%Use matlabs built in ode45 solver 
%inputs are found in @gb_ode_func
[t,z]=ode45(@gb_ode_func, tspan, [x0 v0], opts, p);

% Position and velocity results
%"z" is the output of the ode45 solution
x = z(:,1);
y = z(:,2);
vx = z(:,3);
vy = z(:,4);
p.v=sqrt(vx.^2+vy.^2);

vel=p.v;

%"outy" will output the varriables I'm interested in in a convenient table
%that I can copy into a spreadsheet
outy=[Pgz, N, max(x), r1, umax,t0,xi0,yi0];

end