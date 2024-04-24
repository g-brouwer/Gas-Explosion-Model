function gb = gb_ode_func(t,z,p)

%z(1) & z(2) are position in x & y
%velocity in x direction
dz1 = z(3); 
%velocity in y direction
dz2 = z(4); 

%position
p.ri0=sqrt(p.xi0^2+p.yi0^2);
% Gas velocity decays as:
p.u_air=p.v0*((p.ri0/sqrt(z(1)^2+z(2)^2))^2)*exp(-t/p.tau);
% x & y components
p.u_air_x=p.u_air*cos(p.theta);
p.u_air_y=p.u_air*sin(p.theta);

%account for the drag varriation due to the air velocity relative to the projectile
if (p.u_air_x >= 0 && p.u_air_y >= 0) %if the air is moving
    vx=z(3)-p.u_air_x;
    vy=z(4)-p.u_air_y;
else %the air has stopped expanding
    vx=z(3);
    vy=z(4);
end

%velocity magnitude
p.v=sqrt(vx^2+vy^2);

%the equations of motion in x &  y direction to be solved by ode45
dz3=-(0.5*p.Cd*p.rhoa*p.A/p.m)*vx*p.v;
dz4=-p.g-(0.5*p.Cd*p.rhoa*p.A/p.m)*vy*p.v;

% Zero motion if projectile reaches ground
if (z(2) <= 0 && t > 0)
    dz1 = 0;
    dz2 = 0;
    dz3 = 0;
    dz4 = 0;
end

%outputs that ultimately feed into ode45 (inputted in func_explode)
gb = [dz1; dz2; dz3; dz4];

end