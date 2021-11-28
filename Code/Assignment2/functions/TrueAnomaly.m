function [f] = TrueAnomaly(t, e, a, mu, t0, f0)

% TrueAnomaly.m - computes the true anomaly for the actual time. 
%
% PROTOTYPE:
%   [f] = TrueAnomaly(t, e, a, mu, t0, f0)
%
% INPUT:
%   t                   Actual time                           [s]
%   e                   Eccentricity                          [-]
%   a                   Semi-major axis                       [km]
%   mu                  Planetary gravitational constant      [km^3/s^2]
%   t0                  Initial time                          [s]
%   f0                  Initial true anomaly                  [rad]      
%
% OUTPUT:
%   f                   True anomaly                          [deg]      
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


if nargin < 5
    % initial values of t0 and f0 assumed to be at perigee passage t0 = tp and f0 = 0
    f0 = 0;
end

n = sqrt(mu/a^3);               % mean motion       
M = n*(t-t0);                   % mean anomaly      

% initial Eccentric Anomaly at time t0
E0 = 2*atan2(tan(f0/2),sqrt((1+e)/(1-e)));

% M to m [0,2*pi]rad, k ∈ Z: M = M_prime + k*2*pi
m = wrapTo2Pi(M);
k = floor(M/(2*pi));

% solve Kepler's equation
E_kep = m + (e*sin(m))/(1- sin(m + e) + sin(m));
eq = @(E) E - e*sin(E) - (E0 - e*sin(E0)) - m;

options = optimoptions('fsolve','Display','none');
E_sol = fsolve(eq, E_kep, options);

% corresponding f [0,2*pi] rad
f = 2*atan2d(tan(E_sol/2),sqrt((1-e)/(1+e)));

if f < 0
    f = f + 360;
end

% add number of revolutions

f = f + k*360 ;


end
