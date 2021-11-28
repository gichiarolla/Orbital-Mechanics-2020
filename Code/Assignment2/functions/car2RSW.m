function [RotM] = car2RSW(i,RAAN,omega,f)

% car2RSW.m - computes the rotation matrix to change the reference frame
%             from the inertial equatorial one to another containing 
%             vectors radial, tangential and out of plane.
%
% PROTOTYPE:
%   [RotM] = car2RSW(i,RAAN,omega,f) is a function 
%
% INPUT:
%   i       [1x1]   Inclination         \                             [rad]
%   RAAN    [1x1]   Right Ascension of the Ascending Node             [rad]
%   omega   [1x1]   Argument of Periapsis                             [rad]
%   f       [1x1]   True Anomaly                                      [rad]
%
% OUTPUT:
%   RotM    [3x3]   Rotation Matrix to change frome the inertial equatorial
%                   to the radial transversal out of plane RF [-] 
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version

u = omega + f;
RotM = [-sin(RAAN)*cos(i)*sin(u)+cos(RAAN)*cos(u)       cos(RAAN)*cos(i)*sin(u)+sin(RAAN)*cos(u)        sin(i)*sin(u);
        -sin(RAAN)*cos(i)*cos(u)-cos(RAAN)*sin(u)       cos(RAAN)*cos(i)*cos(u)-sin(RAAN)*sin(u)        sin(i)*cos(u);
        sin(RAAN)*sin(i)                                -cos(RAAN)*sin(i)                               cos(i)];
end