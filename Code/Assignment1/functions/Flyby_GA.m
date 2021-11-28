function [delta,rp,h_ga,DVp,DVtot,em,ep,am,ap,ym,yp] = Flyby_GA(vinf_m,vinf_p,id)

% Flyby_GA.m - Computes the parameters defining a powered gravity assist
%              fly-by given the incoming and outcoming excess velocity.
%
% PROTOTYPE:
%   [delta,rp,h_ga,DVp,DVtot,em,ep,am,ap,ym,yp] = Flyby_GA(vinf_m,vinf_p,planet_id)
%
%  INPUT:
%    vinf_m[3]        Incoming excess velocity [km/s].
%    vinf_p[3]        Outcoming excess velocity [km/s].
%    planet_id[1]     Integer number identifying the planet around which the fly-by.
%
%  OUTPUT:
%    delta[1]  Turning angle [rad]
%    rp[1]     Radius of periapsis [km].
%    h_ga[1]   Altitude of the closest approach for the fly-by [km].
%    DVp[1]    Cost of the powered gravity assist manouvre [km/s].
%    DVtot[1]  Total cost for the fly-by [km/s].
%    em[1]     Eccentricity of the incoming hyperbolic arc.
%    ep[1]     Eccentricity of the outcoming hyperbolic arc.
%    am[1]     Semi-latus rectum of the incoming hyperbolic arc [km].
%    ap[1]     Semi-latus rectum of the incoming hyperbolic arc [km].
%    ym[6]     Vector conteining the radius [km] and velocity [km/s] at
%              periapsis of the incoming hyperbolic arc [rx;ry;rz;vx;vy;vz].
%    yp[6]     Vector conteining the radius [km] and velocity [km/s] at
%              periapsis of the outcoming hyperbolic arc [rx;ry;rz;vx;vy;vz].
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Menicaglia
%
% VERSIONS:
%   2021-01: Last version

% Planet constants
R=astroConstants(20+id);
muP=astroConstants(10+id);

% Computing module of excess velocity
vinf_m_norm=norm(vinf_m);
vinf_p_norm=norm(vinf_p);

% Computing semi-latus rectum
am=-muP/(vinf_m_norm^2);
ap=-muP/(vinf_p_norm^2);

% Computing turning angle
delta=acos(dot(vinf_m,vinf_p)/(vinf_m_norm*vinf_p_norm));

% computing radius of perigee and heigth of the closest approach
rp_fun= @(rp) asin(1/(1+(rp*vinf_m_norm^2)/muP))+asin(1/(1+(rp*vinf_p_norm^2)/muP))-delta;

e=1/sin(delta/2);
rp_guess=am*(1-e);

options = optimoptions('fsolve','Display','none');
rp=fsolve(rp_fun,rp_guess,options);

h_ga=rp-R;

% Computing eccentricity
em=1+rp*vinf_m_norm^2/muP;
ep=1+rp*vinf_p_norm^2/muP;

% Computing total cost and graviy assist cost
DVtot=norm(vinf_p-vinf_m);

[rm,vm] = kep2car(am,em,0,0,0,0,muP);
[r_p,vp] = kep2car(ap,ep,0,0,0,0,muP);

DVp=norm(vp-vm);

% For plot
ym=[rm;vm];
yp=[r_p;vp];

end