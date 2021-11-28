function [tgauss,kepgauss] = propagatorgauss(tspan,kep0,kP,rP,J2,CD,AtoMratio,rotP)

% propagatorgauss.m - computes the time evolution of the keplerian elements 
%                     for an orbit in which both the second zonal harmonic 
%                     and the drag due to atmosphere are important perturbations. 
%
% PROTOTYPE:
%   [tgauss,kepgauss] = propagatorgauss(tspan,kep0,kP,rP,J2,CD,AtoMratio,rotP)
%
% INPUT:
%   tspan       (nx1)   Vector of time at which the keplerian elements  [s]
%                       will be evaluated
%   kep0        (1x6)   Vector of initial keplerian elements  [km,-,rad,rad,rad,rad]
%   kP          (1x1)   Planetary gravitational constant                [km^3/s^2]
%   rP          (1x1)   Planetary radius                                [km]
%   J2          (1x1)   Second zonal harmonic                           [-]
%   CD          (1x1)   Drag coefficient                                [-]
%   AtoMratio   (1x1)   Area to mass ratio for the spacecraft           [km^2/kg]
%   rotP        (1x1)   Rotation rate of the planet                     [rad/s]
%
% OUTPUT:
%   tgauss      (nx1)   Vector of time at which the keplerian elements  [s]
%                       have been evaluated
%   kepgauss    (nx6)   Matrix of the keplerian elements evaluated at
%                       every given time instant              [km,-,rad,rad,rad,rad]
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version

% Definition of anonymous function containing both drag and J2
% perturbations
a_tot = @(t,kep) a_tot_RSW(t,kep,kP,rP,J2,CD,AtoMratio,rotP);

% Ordinary Differential Equation for the orbit in the RSW frame
odefun = @(t,kep) EoMRSW(t,kep,kP,a_tot);

% Solution of the differential equation
options = odeset('RelTol',1e-13,'AbsTol',1e-11);
[tgauss,kepgauss] = ode113(odefun,tspan,kep0,options);

end