function ap_RSW = a_J2_RSW(t,kep,kP,rP,J2)


% a_J2_RSW.m - evaluates acceleration in a radial-transversal-out of plane 
%              reference frame. Said accelerations will contain the effects 
%              of J2 on the equations of motion and will have to be 
%              introduced in a KEPLERIAN element problem.
%
% PROTOTYPE:
%   ap_RSW = a_J2Gauss_RSW(t,kep,kP,rP,J2) 
%
% INPUT:
%   t     (1x1)   Time step                                       [s]
%   kep   (6x1)   Keplerian elements              [km,-,rad,rad,rad,rad,rad]
%   kP    (1x1)   Planetary gravitational constant                [km^3/s^2]
%   rP    (1x1)   Planetary radius                                [km]
%   J2    (1x1)   Second zonal harmonic                           [-]
%
% OUTPUT:
%   aJ2_RSW (3x1) Accelerations due to perturbations via J2 on the  [km/s^2]
%                 radial-transversal-oop RF 
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version

% Extraction of keplerian elements from state vector
a = kep(1);
e = kep(2);
i = kep(3);
om = kep(5);
th = kep(6);

% Definition of useful parameters
p = a*(1-e^2);
r = p/(1+e*cos(th));

% Acceleration, perturbed by the second zonal harmonic
ap_RSW = -1.5*J2*kP*rP^2/r^4.*[1 - 3*sin(i)^2*sin(th+om)^2;
                               sin(i)^2*sin(2*(om+th));
                               sin(2*i)*sin(om+th)];
end

