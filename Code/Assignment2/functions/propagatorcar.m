function [tcar,rcar] = propagatorcar(tspan,condition0,kP,rP,J2,CD,AtoMratio,rotP)

% propagatorcar.m - computes the time evolution of the cartesian elements 
%                   for an orbit in which both the second zonal harmonic 
%                   and the drag due to atmosphere are important perturbations. 
%
% PROTOTYPE:
%   [tgauss,kepgauss] = propagatorgauss(tspan,kep0,kP,rP,J2,CD,AtoMratio,rotP)
%
% INPUT:
%   tspan       (nx1)   Vector of time at which the keplerian elements  [s]
%                       will be evaluated
%   condition0  (6x1)   Vector of initial cartesian elements  [km,km,km,km/s,km/s,km/s]
%   kP          (1x1)   Planetary gravitational constant                [km^3/s^2]
%   rP          (1x1)   Planetary radius                                [km]
%   J2          (1x1)   Second zonal harmonic                           [-]
%   CD          (1x1)   Drag coefficient                                [-]
%   AtoMratio   (1x1)   Area to mass ratio for the spacecraft           [km^2/kg]
%   rotP        (1x1)   Rotation rate of the planet                     [rad/s]
%
% OUTPUT:
%   tcar        (nx1)   Vector of time at which the cartesian elements  [s]
%                       have been evaluated
%   rcar        (nx6)   Matrix of the cartesian elements evaluated at
%                       every given time instant	      [km,km,km,km/s,km/s,km/s]
%                       
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


options = odeset('RelTol',1e-13,'AbsTol',1e-11);

odefun = @(t,r) EoMXYZ_dragJ2(t,r,kP,rP,J2,CD,AtoMratio,rotP);

[tcar,rcar] = ode113(odefun,tspan,condition0,options);

end