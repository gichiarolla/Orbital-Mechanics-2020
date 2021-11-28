function dy= eq_of_motion(~,y,mu)

% eq_of_motion - ODE system for the two body problem (Keplerian motion)
%
% PROTOTYPE:
%   dy= eq_of_motion(~,y,mu)
% 
% DESCRIPTION:
%   Implicit equation of motion to be numerically integrated in order to
%   propagate the orbit.
%   
%  INPUT:
%    t[1]    Time [s] (can be omitted, as the system is autonomous).
%    y[6]    Vector containing position [km] and velocity [km/s] at time t [rx,ry,rz,vx,vy,vz]'.
%    mu[1]   Gravitational constant of the planet [km^3/s^2].
%
%  OUTPUT:
%    dy[6]  Vector containing the prime and the second derivative of position.
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
%   2020-09-24: First version
%   2020-12-18: Second version

dy=[y(4),y(5),y(6),-mu/norm(y(1:3))^3*y(1),-mu/norm(y(1:3))^3*y(2),-mu/norm(y(1:3))^3*y(3)]';

end