function [a_new] = rep_orbit(om_E,m,k,mu)

% rep_orbit.m - computes the new semi-major axis in order to get 
%               a repeated ground track. 
%
% PROTOTYPE:
%   [a_new] = rep_orbit(om_E,m,k,mu)
%
% INPUT:
%   om_E                Rotation velocity of Earth            [rad/s]
%   m                   Rotations of the planet
%   k                   revolutions of the satellite
%   mu                  Planetary gravitational constant      [km^3/s^2]
%
% OUTPUT:
%   a_new               New semi-major axis                    [km]
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


n=om_E*(k/m);
a_new = (mu/(n^2))^(1/3);

return