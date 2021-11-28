function dY = ode_motion_equation(~,Y,mu)

% ode_motion_equation.m - ODE system for the equations of motion of the obj
%
% PROTOTYPE:
%   dY = ode_motion_equation(~,Y,mu)
%
% INPUT:
%   Y   [6x1] Cartesian state of the body (rx,ry,rz,vx,vy,vz) [km, km/s]
%   mu  [1]   Gravitational parameter                         [km^3/s^2]
%
% OUTPUT:
%   dy  [6x1] Derivative of the state                         [km/s, km/s^2]
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


r = Y(1:3); %3x1
v = Y(4:6); %3x1
                                   

% State vector 6x1
dY=[v
    r.*(-mu/(norm(r))^3)];
    

end