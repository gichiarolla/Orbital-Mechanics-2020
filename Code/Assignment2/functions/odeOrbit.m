function dY = odeOrbit(~,Y,mu,J2,Re)

% odeOrbit.m - ODE system for the equations of motion of the obj,
%              with peturbations
%
% PROTOTYPE:
%   dY = odeOrbit(~,Y,mu,J2,Re)
%
% INPUT:
%   Y   [6x1] Cartesian state of the body (rx,ry,rz,vx,vy,vz) [km, km/s]
%   mu  [1]   Gravitational parameter                         [km^3/s^2]
%   J2        Second zonal harmonic                      
%   RE        Earth's radius
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


r = Y(1:3);
v = Y(4:6);


% second zonal harmonic
a_J2 = ((1.5*J2*mu*Re^2)/(norm(r)^4))*[(r(1)/norm(r))*(5*(r(3)/norm(r))^2 - 1)
                                        (r(2)/norm(r))*(5*(r(3)/norm(r))^2 - 1)
                                        (r(3)/norm(r))*(5*(r(3)/norm(r))^2 - 3)];
                                        

% State vector
dY=[v
    r*(-mu/(norm(r))^3)+ a_J2];
    

end