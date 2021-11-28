function [lon,lat]=groundTrack(a,e,i,OM,om,f,mu,om_E,t0,tspan,J2,RE)

% groundTrack.m - plot the unperturbed ground track.
%
% PROTOTYPE:
%   [lon,lat]=groundTrack(a,e,i,OM,om,f,mu,om_E,t0,tspan,J2,RE)
%
% INPUT:
%   a                   Semi-major axis                       [km]
%   e                   Eccentricity                          [-]
%   i                   Inclination                           [rad]
%   OM                  RAAN                                  [rad]
%   om                  Argument of periapsis                 [rad]
%   f                   True anomaly                          [rad]      
%   mu                  Planetary gravitational constant      [km^3/s^2]
%   om_E                Rotation velocity of Earth            [rad/s]
%   t0                  Initial time                          [s]
%   tspan       (1x2)   Vector of time                        [s]
%   J2                  Second zonal harmonic                      
%   RE                  Earth's radius                        [km]
%
% OUTPUT:
%   Lon         (1xn)   Longitude                             [rad]
%   Lat         (1xn)   Latitude                              [rad]
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


[r0,v0] = kep2car(a,e,i,OM,om,f, mu);
Y0 =[r0; v0];

options = odeset( 'Reltol', 1e-13, 'Abstol', 1e-14);
[ T, X ] = ode45( @(t,Y) odeOrbit(t,Y,mu,J2,RE), tspan, Y0,options);

alpha= zeros(1,length(T));
delta= zeros(1,length(T));
lon= zeros(1,length(T));
lat= zeros(1,length(T));

for k = 1:length(T)
    r = X(k,1:3);
    delta(k) = asin(r(3)/norm(r));
    
    if (r(2)/norm(r))>0
        alpha(k) = acos(r(1)/(norm(r)*cos(delta(k))));
    else
        alpha(k) = 2*pi - acos(r(1)/(norm(r)*cos(delta(k))));
    end
    
    lon(k) = alpha(k) - om_E*(T(k)-t0);
    lat(k) = delta(k);
end

return
