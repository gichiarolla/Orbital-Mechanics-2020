function [rvect,vvect] = kep2car(a,e,i,RAAN,omega,theta,kP)

% kep2car.m - converts keplerian parameters in radius and velocity written
%             in cartesian coordinates in the inertial planetocentric RF.
%
% PROTOTYPE:
%   [rvect,vvect] = kep2car(a,e,i,Omega,omega,theta,mu)
%
% INPUT:
%   a       [1x1]   Semi-major axis                              [km]
%   e       [1x1]   Eccentricity                                 [-]
%   i       [1x1]   Inclination                                  [rad]
%   RAAN    [1x1]   Right Ascension of the Ascending Node        [rad]
%   omega   [1x1]   Argument of Periapsis                        [rad]
%   theta   [1x1]   True Anomaly                                 [rad]
%   kP      [1x1]   Planetary gravitational constant             [km^3/s^2]
%
% OUTPUT:
%   rvect   [3x1]   Radius vector in the inertial planetocentric RF    [km]
%   vvect   [3x1]   Velocity vector in the inertial planetocentric RF  [km]
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version

% Computation of the radius from the keplerian elements
p = a*(1-e^2);
r = p/(1+e*cos(theta));

% Radius components in the perifocal RF
re = r*cos(theta);
rp = r*sin(theta);

% Velocity components in the perifocal RF
h = sqrt(kP*p);
ve = -kP/h*sin(theta);
vp = kP/h*(e+cos(theta));

% Rotation matrices from the inertial planetocentric RF to the perifocal
% RF, they will be transposed to obtain the rotation from perifocal to
% inertial planetocentric
% Rotation of RAAN around K
R_RAAN = [cos(RAAN)   sin(RAAN)  0
           -sin(RAAN)  cos(RAAN)  0
           0            0           1];
% Rotation of i around N
R_i = [1        0           0
       0        cos(i)      sin(i)
       0        -sin(i)     cos(i)];
% Rotation of omega around h   
R_omega = [cos(omega)   sin(omega)  0
           -sin(omega)  cos(omega)  0
           0            0           1];
% Rotation matrix and its transpose
R = R_omega*R_i*R_RAAN;
Rot_inv = R';

% Rotation of r and v as vectors in the perifocal RF to the planetocentric
% inertial RF
rpf = [re;  rp; 0];
vpf = [ve;  vp; 0];
rvect = Rot_inv*rpf;
vvect = Rot_inv*vpf;

end