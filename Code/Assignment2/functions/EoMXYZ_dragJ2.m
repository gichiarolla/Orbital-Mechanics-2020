function dr = EoMXYZ_dragJ2(t,r,kP,rP,J2,CD,AtoMratio,rotP)

% EoMXYZ_dragJ2.m - set of ordinary differential equations in the time 
%                   dependent variables radius and its derivative.
%                   This function accounts for both the second zonal
%                   harmonic of the primary planet (which will be the Earth)  
%                   and the aerodynamic drag due to Earth's atmosphere
%
% PROTOTYPE:
%   dr = EoMXYZ_dragJ2(t,r,kP,rP,J2,CD,AtoMratio,rotP)
%
% INPUT:
%   t          (1x1)   Time used as an independent variable            [s]
%   r          (6x1)   Radius and its time derivative (velocity) as    [km]
%                      time-dependent variables                        [km/s]  
%   kP         (1x1)   Planetary gravitational constant                [km^3/s^2]
%   rP         (1x1)   Planetary radius                                [km]
%   J2         (1x1)   Second zonal harmonic                           [-]
%   CD         (1x1)   Drag coefficient                                [-]
%   AtoMratio  (1x1)   Area to mass ratio for the spacecraft           [km^2/kg]
%   rotP       (1x1)   Rotational speed of the planet                  [rad/s]
%
% OUTPUT:
%   dr         (6x1)   Time derivative of radius and velocity          [km/s]
%                      (velocity and acceleration)                     [km/s^2]
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


% Computation of radius magnitude
rmag = norm(r(1:3));

% Computation of the acceleration perturbation
[a,e,i,RAAN,omega,theta] = car2kep(r(1:3)',r(4:6)',kP);
kep = [a,e,i,RAAN,omega,theta];
[~,a_drag_XYZ] = a_drag(t,kep,kP,rP,CD,AtoMratio,rotP);
                                   
% Time derivative of radius and velocity
dr = [r(4);
      r(5);
      r(6);
      -kP*r(1)/(rmag)^3+(1.5*(J2*kP*rP^2)/rmag^4)*(r(1)/rmag)*(5*(r(3)/rmag)^2-1)+a_drag_XYZ(1);
      -kP*r(2)/(rmag)^3+(1.5*(J2*kP*rP^2)/rmag^4)*(r(2)/rmag)*(5*(r(3)/rmag)^2-1)+a_drag_XYZ(2);
      -kP*r(3)/(rmag)^3+(1.5*(J2*kP*rP^2)/rmag^4)*(r(3)/rmag)*(5*(r(3)/rmag)^2-3)+a_drag_XYZ(3)];
end