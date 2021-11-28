function [a_drag_RSW,a_drag_XYZ] = a_drag(t,kep,kP,rP,CD,AtoMratio,rotP)

% a_drag.m - Defines the perturbation for acceleration due to aerodynamic drag in the
%            radial-transversal out of plane reference frame, having the keplerian
%            elements as an input
%
% PROTOTYPE:
% [a_drag_RSW,a_drag_XYZ] = drag_perturbationRSW(t,kep,kP,rP,CD,AtoMratio,rotP)
%
% INPUT:
%   t         (nx1)   Time used as an independent variable            [s]
%   kep       (1x6)   Keplerian elements as a variable dependent on
%                     time                                            [km,-,rad,rad,rad,rad]
%   kP        (1x1)   Planetary gravitational constant                [km^3/s^2]
%   rP        (1x1)   Planetary radius                                [km]
%   J2        (1x1)   Second zonal harmonic                           [-]
%   CD        (1x1)   Drag coefficient                                [-]
%   AtoMratio (1x1)   Area to mass ratio for the spacecraft           [km^2/kg]
%   rotP      (1x1)   Rotational speed of the planet                  [rad/s]
%
% OUTPUT:
%   a_drag_RSW (3x1)  Acceleration due to aerodynamic drag in the     [km^2/s^2] 
%                     radial-transversal-out of plane RF
%   a_drag_XYZ (3x1)  Acceleration due to aerodynamic drag in the     [km/s^2]
%                     planetocentric inertial cartesian RF
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version

% Conversion to cartesian coordinates
[r,v] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),kP);
% Computation of the magnitude of the radius
rmag = norm(r);

% Computation of the altitude, as an input variable for a function 
% computing air density at said altitude
h = rmag-rP;
rho = densityAltitude(t,h);
% Comversion in kilometers
rho_km = rho*1e9;

% Computation of the relative velocity
omegaE = [0;0;rotP];
v_rel = v-cross(omegaE,r);
v_relnorm = norm(v_rel);

% Computation of the acceleration in cartesian coordinates
a_drag_XYZ = -0.5*AtoMratio*CD*rho_km*v_relnorm^2*v_rel/v_relnorm;
% Computation of the acceleration in RSW coordinates
RotM = car2RSW(kep(3),kep(4),kep(5),kep(6));
a_drag_RSW = RotM*a_drag_XYZ;

end