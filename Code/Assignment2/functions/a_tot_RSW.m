function acc = a_tot_RSW(t,kep,kP,rP,J2,CD,AtoMratio,rotP)

% a_tot_RSW.m - computes the total acceleration due to perturbations 
%               (considering only drag and J2) in the
%               radial transversal out of plane RF.
%
% PROTOTYPE:
%   acc = a_tot_RSW(t,kep,kP,rP,J2,CD,AtoMratio,rotP) 
%
% INPUT:
%   t          (nx1)   Time instant at which the acceleration will be [s]
%                      computed
%   kep0       (1x6)   Vector of initial keplerian elements
%                      [km,-,rad,rad,rad,rad]
%   kP         (1x1)   Planetary gravitational constant              [km^3/s^2]
%   rP         (1x1)   Planetary radius                              [km]
%   J2         (1x1)   Second zonal harmonic                         [-]
%   CD         (1x1)   Drag coefficient                              [-]
%   AtoMratio  (1x1)   Area to mass ratio for the spacecraft         [km^2/kg]
%   rotP       (1x1)   Rotation rate of the planet                   [rad/s]
%
% OUTPUT:
%   acc        (3x1)   Vector of acceleration components in the RSW  [km/s^2]
%                       reference frame
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version

[a_dragRSW,~] = a_drag(t,kep,kP,rP,CD,AtoMratio,rotP);
% Definition of the acceleration due to the second zonal harmonic as a
% function of time and keplerian elements
a_J2RSW = a_J2_RSW(t,kep,kP,rP,J2);
% Total acceleration due to both perturbations
acc = a_dragRSW+a_J2RSW;

end