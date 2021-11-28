function [] = GT_comparison (a,e,i,OM,om,f0,muE,om_E,t0,T,m,k,J2,RE)

% GT_comparison.m - plot a comparison of GT depending on orbit's 
%                   duration or type of GT.
%
% PROTOTYPE:
%   [] = GT_comparison (a,e,i,OM,om,f,muE,om_E,t0,T,m,k,J2,RE)
%
% INPUT:
%   a                   Semi-major axis                       [km]
%   e                   Eccentricity                          
%   i                   Inclination                           [rad]
%   OM                  RAAN                                  [rad]
%   om                  Argument of periapsis                 [rad]
%   f0                  True anomaly                          [rad]      
%   muE                 Planetary gravitational constant      [km^3/s^2]
%   om_E                Rotation velocity of Earth            [rad/s]
%   t0                  Initial time                          [s]
%   T                   Earth's orbit period                  [s]
%   m                   Rotations of the planet
%   k                   revolutions of the satellite
%   J2                  Second zonal harmonic                      
%   RE                  Earth's radius                        [km]
%
% OUTPUT:
%   Figure with ground tracks' rapresentation
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


% a) Ground Track
% not perturbed 
T_1Day_10Days_Unp(a,e,i,OM,om,f0,muE,om_E,t0,T,RE,'all');
    
% secular
T_1Day_10Days_Sec(a,e,i,OM,om,f0,muE,om_E,t0,T,m,k,J2,RE,'all');


% b) Repeating Ground Track
% not perturbed
T_1Day_10Days_Rep(a,e,i,OM,om,f0,muE,om_E,t0,T,m,k,RE,'all');

% secular
T_1Day_10Days_RepSec(a,e,i,OM,om,f0,muE,om_E,t0,T,m,k,J2,RE,'all');


% One figure for each period of orbit 

GT_Sec_Rep_RepSec(a,e,i,OM,om,f0,muE,om_E,t0,T,m,k,J2,RE);

return