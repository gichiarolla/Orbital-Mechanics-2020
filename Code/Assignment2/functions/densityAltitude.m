function rho = densityAltitude(t,h)

% densityAltitude.m - computes Earth atmosphere density given time 
%             and altitude of a satellite during its orbit.
%
% PROTOTYPE:
%   rho = densityAltitude(t,h) is a function that 
%
% INPUT:
%   t     [nx1]   Time as an independent variable                     [s]
%   h     [nx1]   Altitude from Earth's surface                       [km]
%
% OUTPUT:
%   rho   [nx1]   Atmosphere density                              [kg/m^3]
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


% Vectors containing given informations for a the interpolation to
% compute the density
href = [0;25;30;40;50;60;70;80;90;100;110;120;130;140;150;180;200;250;300;
        350;400;450;500;600;700;800;900;1000];          % [km]
rhoref = [1.225;3.899*1e-2;1.774*1e-2;3.972*1e-3;1.057*1e-3;3.206*1e-4;
          8.770*1e-5;3.396*1e-6;5.297*1e-7;9.661*1e-8;2.438*1e-8;
          8.484*1e-9;3.845*1e-9;2.070*1e-9;5.464*1e-10;7.248*1e-11;
          2.418*1e-11;9.158*1e-12;3.725*1e-12;1.585*1e-12;6.967*1e-13;
          1.454*1e-13;3.614*1e-14;1.170*1e-14;5.245*1e-15;
          3.019*1e-15];                                 % [kg/m^3]
H = [7.249;6.349;6.682;7.554;8.382;7.714;6.549;5.799;5.382;5.877;7.263;
    9.473;12.636;16.149;22.523;29.740;37.105;45.546;53.628;53.298;58.515;
    60.828;63.822;71.835;88.667;124.64;181.05;268];     % [km]

% Cycle relating the given altitude to the right height interval and
% computing its respective density. 1600km is the altitude at which the
% code considers the atmosphere completely negligible.
if h > href(end) && h < 1600
    rho = rhoref(end)*exp(-(h-href(end))/H(end));
elseif h >= 1600
    rho = 0;
else 
    k = 1;
    while h > href(k+1) 
        rho = rhoref(k)*exp(-(h-href(k))/H(k));
        k = k+1;
    end
end
end