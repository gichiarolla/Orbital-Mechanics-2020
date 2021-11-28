function [a,e,i,RAAN,omega,theta] = car2kep(rvect,vvect,kP)

% car2kep.m - transforms the radius and velocity vectors, written in 
%             cartesian coordinates in the 
%             planetocentric inertial RF into keplerian orbital parameters.
%
% PROTOTYPE:
%   [a,e,i,RAAN,omega,theta,M] = car2kep(rvect,vvect,kP) 
%
% INPUT:
%   rvect   [3x1]   Radius vector in the inertial planetocentric RF    [km]
%   vvect   [3x1]   Velocity vector in the inertial planetocentric RF  [km]
%   kP      [1x1]   Planetary gravitational constant                   [km^3/s^2]
%
% OUTPUT:
%   a       [1x1]   Semi-major axis                                    [km]
%   e       [1x1]   Eccentricity                                       [-]
%   i       [1x1]   Inclination                                        [rad]
%   RAAN    [1x1]   Right Ascension of the Ascending Node              [rad]
%   omega   [1x1]   Argument of Periapsis                              [rad]
%   theta   [1x1]   True Anomaly                                       [rad]
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


% Earth is assumed as a planet if the planetary constant isn't between the
% inputs
if nargin == 2
    kP = 398600.44;
end

% Velocity and radius magnitudes are computed in order to find specific
% orbital energy and then the Semi Major Axis
r = norm(rvect);
v = norm(vvect);
E = v^2/2 - kP/r;
a = -kP/(2*E);

% Computation of the specific angular momentum vector and its magnitude 
% to compute e, i, RAAN
hvect = cross(rvect,vvect);
h = norm(hvect,2);
% Eccentricity
evect = cross(vvect,hvect)/kP - rvect/r;
e = norm(evect,2);
% Inclination
i = acos(hvect(1,3)/h);
% Line of Nodes and its magnitude
Nvect = cross([0,0,1],hvect);
N = norm(Nvect,2);
% Right Ascension of the Ascending Node, assumed null when inclination is
% zero
if i == 0
    RAAN = 0;
else
    if Nvect(1,2) >= 0
        RAAN = acos(Nvect(1,1)/N);
    else
        RAAN = 2*pi - acos(Nvect(1,1)/N);
    end
end

% Argument of periapsis
if e == 0
    omega = 0;
else
    if evect(1,3) >= 0
        omega = acos(Nvect*evect'/(N*e));
    else
        omega = 2*pi - acos(Nvect*evect'/(N*e));
    end
end

% True Anomaly
v_rad = rvect*vvect'/r;
if v_rad >= 0
    theta = acos(evect*rvect'/(e*r));
else
    theta = 2*pi - acos(evect*rvect'/(e*r));
end

end