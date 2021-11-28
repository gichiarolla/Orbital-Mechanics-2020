function [r,v] = kep2car(a,e,i,OMG,omg,f,mu)

% kep2car - Conversion of keplerian orbit elements into cartesian state vectos.
%
% PROTOTYPE:
%   [r,v] = kep2car(a,e,i,OMG,omg,f,mu)
% 
%  INPUT:
%    a[1]    Semi-latus rectum [km].
%    e[1]    Eccemtricity.
%    i[1]    Inclination [rad].
%    OMG[1]  Right ascension of the ascending node [rad].
%    omg[1]  Argument of perigee [rad].
%    f[1]    True anomaly [rad].
%    mu[1]   Gravitational constant of the planet [km^3/s^2].
%
%  OUTPUT:
%    r[3]  position [rx;ry;rz] [km].
%    v[3]  velocity [vx;vy;vz] [km/s].
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

% Distance rr=norm(r)
rr=a*(1-e^2)/(1+e*cos(f));

% Specific angular momentum h
h=sqrt(mu*a*(1-e^2));

% Rotation matrixes
R3_OMG=[cos(-OMG),sin(-OMG),0;-sin(-OMG),cos(-OMG),0;0,0,1];
R1_i=[1,0,0;0,cos(-i),sin(-i);0,-sin(-i),cos(-i)];
R3_omg=[cos(-omg),sin(-omg),0;-sin(-omg),cos(-omg),0;0,0,1];

% Perifocal frame vectors
rp=rr*[cos(f);sin(f);0];
vp=mu/h*[-sin(f);e+cos(f);0];

% r & v
r=R3_OMG*R1_i*R3_omg*rp;
v=R3_OMG*R1_i*R3_omg*vp;

end