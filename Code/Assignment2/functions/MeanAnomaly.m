function M = MeanAnomaly(e,f)

% MeanAnomaly.m - computes the value of the mean anomaly.
%
% PROTOTYPE:
%   M = MeanAnomaly(e,f)
%
% INPUT:
%   e       [1x1]   Eccentricity                          [-]
%   f       [1x1]   True Anomaly                          [rad]
%
% OUTPUT:
%   M               Mean anomaly                          [rad]      
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


E = 2*atan(sqrt((1-e)/(1+e))*tan(f/2));
M = E - e*sin(E);

end