function [a,e,i,OMG,omg,f] = heliocentric_trajectories(r_dep,r_fb,r_arr,vt1_arc1,vt2_arc1,vt1_arc2,vt2_arc2,muSun)

% timeoffligth_hyperbolic.m - computes the Keplerian elements of the
%                             heliocentric transfer arcs for the
%                             interplanetary explorer mission
%
% PROTOTYPE:
%   [a,e,i,OMG,omg,f] = heliocentric_trjectories(r_dep,r_fb,r_arr,vt1_arc1,vt2_arc1,vt1_arc2,vt2_arc2,muSun) 
%
%  INPUT:
%    r_dep[3x1]    Position at departure [km].
%    r_fb[3x1]     Position at fly-by [km].  
%    r_arr[3x1]    Position at arrival [km].  
%    vt1_arc1[3x1] velocity at departure [km/s].
%    vt2_arc1[3x1] velocity at fly-by for the first arc [km/s].
%    vt1_arc2[3x1] velocity at fly-by for the second arc [km/s].
%    vt2_arc2[3x1] velocity at arrival [km/s].
%    muSun[1]      Gravitational constant of the Sun [km^3/s^2].
%
%  OUTPUT:
%    a[1x4]    Semi-latus rectum [km].
%    e[1x4]    Eccemtricity.
%    i[1x4]    Inclination [rad].
%    OMG[1x4]  Right ascension of the ascending node [rad].
%    omg[1x4]  Argument of perigee [rad].
%    f[1x4]    True anomaly [rad].
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

[a(1),e(1),i(1),OMG(1),omg(1),f(1)]=car2kep(r_dep,vt1_arc1,muSun);
[a(2),e(2),i(2),OMG(2),omg(2),f(2)]=car2kep(r_fb,vt2_arc1,muSun);
[a(3),e(3),i(3),OMG(3),omg(3),f(3)]=car2kep(r_fb,vt1_arc2,muSun);
[a(4),e(4),i(4),OMG(4),omg(4),f(4)]=car2kep(r_arr,vt2_arc2,muSun);
end

