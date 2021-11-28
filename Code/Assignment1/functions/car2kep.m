function [a,e,i,OMG,omg,f] = car2kep(r,v,mu)

% kep2car - Conversion of cartesian state vectosinto Keplerian elements
%
% PROTOTYPE:
%   [r,v] = kep2car(a,e,i,OMG,omg,f,mu)
% 
%  INPUT:
%    r[3x1]  position [rx;ry;rz] [km].
%    v[3x1]  velocity [vx;vy;vz] [km/s].
%    mu[1]   Gravitational constant of the planet [km^3/s^2].
%
%  OUTPUT:
%    a[1]    Semi-latus rectum [km].
%    e[1]    Eccemtricity.
%    i[1]    Inclination [rad].
%    OMG[1]  Right ascension of the ascending node [rad].
%    omg[1]  Argument of perigee [rad].
%    f[1]    True anomaly [rad].
%
%  FUNCTIONS CALLED:
%    (none)
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Menicaglia
%
% VERSIONS:
%   2021-01: Last version

% Preliminary computations
h=cross(r,v); % Orbital momentum vector
e_vect=1/mu*cross(v,h)-r/norm(r); % eccentricity vector
N=cross([0,0,1],h); % vector pointing in the ascending node direction
Nv=N/norm(N); % verson of N 

% True anomaly f
if (v*r'/norm(r))>=0
    f=acos(e_vect/norm(e_vect)*r'/norm(r));
else
    f=2*pi-acos(e_vect/norm(e_vect)*r'/norm(r));
end

% Orbit inclination i
i=acos(h(3)/norm(h));

% Eccentricity e
e=norm(e_vect);

% Longitude of the ascending node OMG
if any(isnan(Nv))
    if e_vect(1)<0
        OMG=pi;
    else 
        OMG=0;
    end
else
    if Nv(2)>=0
        OMG=acos(Nv(1));
    else
        OMG=2*pi-acos(Nv(1));  
    end
end
        
% Argument of the periapsis omg
if e_vect(3)==0
    omg=0;
elseif e_vect(3)>=0
    omg=acos((N*e_vect')/(norm(N)*norm(e_vect)));
else
    omg=2*pi-acos((N*e_vect')/(norm(N)*norm(e_vect)));
end

% Semi-major axis a
a=1/(2/norm(r)-norm(v)^2/mu);

end

