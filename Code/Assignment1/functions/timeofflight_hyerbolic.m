function Dt_tot = timeofflight_hyerbolic(id,ym,yp,em,ep,am,ap)

% timeoffligth_hyperbolic.m - Computes the time duration of a fly-by 
%                             (considering a finite SOI) using the
%                             hyperbolic time law.
%
% PROTOTYPE:
%   [Dt_tot] = timeofflight_hyerbolic(id,ym,yp,em,ep,am,ap) 
%
%  INPUT:
%    id[1]    Integer number identifying the planet around which the fly-by
%             is performed: Mercury(1),Venus(2),Earth(3),Mars(4),Jupiter(5),Saturn(6),Uranus(7),Neptune(8)
%    ym[6]    Vector conteining the radius [km] and velocity [km/s] at
%             periapsis of the incoming hyperbolic arc [rx;ry;rz;vx;vy;vz].
%    yp[6]    Vector conteining the radius [km] and velocity [km/s] at
%             periapsis of the outcoming hyperbolic arc [rx;ry;rz;vx;vy;vz].
%    em[1]    Eccentricity of the incoming hyperbolic arc.
%    ep[1]    Eccentricity of the outcoming hyperbolic arc.
%    am[1]    Semi-latus rectum of the incoming hyperbolic arc [km].
%    ap[1]    Semi-latus rectum of the incoming hyperbolic arc [km].
% 
%  OUTPUT:
%    Dt_tot[1]  Time duration of the flyby [s].
%    
%  FUNCTIONS CALLED:
%    astroConstants.m
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Menicaglia
%
% VERSIONS:
%   2021-01: Last version

% Computation of the radius of the SOI
SOI_rR=[46.1,101.7,145.3,170,675.1,906.9,2024.6,3494.8]; % vector cointeining the ratios r_SOI/Rp fot the planets of the Solar system
R=astroConstants(21:28); % Vector conteining the radius of planets of the Solar system 
mu=astroConstants(11:18); % Vector conteining the gravity constants of planets of the Solar system 
r=R.*SOI_rR; % r=r_SOI

% Computation of angular momentum
hm=norm(cross(ym(1:3),ym(4:6)));
hp=norm(cross(yp(1:3),yp(4:6)));

% Computation of true anomaly
fm=2*pi-acos(1/em*(hm^2/(mu(id)*r(id))-1));
fp=acos(1/ep*(hp^2/(mu(id)*r(id))-1));

% Computation of eccentric anomaly
Fm=2*atan(sqrt((em-1)/(em+1))*tan(fm/2));
Fp=2*atan(sqrt((ep-1)/(ep+1))*tan(fp/2));

% Computation of time
Dtm=sqrt(-am^3/mu(id))*(-em*sinh(Fm)-Fm);
Dtp=sqrt(-ap^3/mu(id))*(em*sinh(Fp)+Fp);

Dt_tot=Dtm+Dtp;

end