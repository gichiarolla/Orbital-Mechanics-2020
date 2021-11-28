function [a_sec] = secular_J2(a,e,i,om_E,m,k,mu,J2,RE)

% secular_J2.m - computes the new semi-major axis due to J2 influence. 
%
% PROTOTYPE:
%   [a_sec] = secular_J2(a,e,i,om_E,m,k,mu,J2,RE)
%
% INPUT:
%   a                   Semi-major axis                       [km]
%   e                   Eccentricity                          [-]
%   i                   Inclination                           [rad]
%   om_E                Rotation velocity of Earth            [rad/s]
%   m                   Rotations of the planet
%   k                   revolutions of the satellite
%   mu                  Planetary gravitational constant      [km^3/s^2]
%   J2                  Second zonal harmonic                      
%   RE                  Earth' radius			      [km]
%
% OUTPUT:
%   a_sec               New semi-major axis                   [km]
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


dOM = cos(i);
dom = (2.5*(sin(i)^2)-2);
dM0 = (1-1.5*sin(i)^2);

f = @(a) m*(sqrt(mu/(a^3)) + (-1.5*(sqrt(mu)*J2*RE^2)/(((1-e^2)^2)*a^3.5))*dom + ...
    (1.5*(sqrt(mu)*J2*RE^2)/(((1-e^2)^1.5)*a^3.5))*dM0) ...
    -k*(om_E - (-1.5*(sqrt(mu)*J2*RE^2)/(((1-e^2)^2)*a^3.5))*dOM);

a_sec = fzero(f,a);

return
