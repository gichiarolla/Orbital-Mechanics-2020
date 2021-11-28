function [] = plot3D_initial_orbit (a,e,i,OM,om,f, muE, T, RE)

% plot3D_initial_orbit.m - plot initial S/C orbit.
%
% PROTOTYPE:
%   [] = plot3D_initial_orbit (a,e,i,OM,om,f, muE, T, RE)
%
% INPUT:
%   a                   Semi-major axis                       [km]
%   e                   Eccentricity                          [-]
%   i                   Inclination                           [rad]
%   OM                  RAAN                                  [rad]
%   om                  Argument of periapsis                 [rad]
%   f                   True anomaly                          [rad]      
%   muE                 Planetary gravitational constant      [km^3/s^2]
%   T                   Earth's orbit period                  [s]
%   RE                  Earth' radius			      [km]
%
% OUTPUT:
%   Figure with initial S/C orbit
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


[r_0,v_0] = kep2car(a,e,i,OM,om,f, muE);
Y_0 =[r_0; v_0];    %6x1
tsp = linspace( 0, T, 10000);

options = odeset( 'Reltol', 1e-13, 'Abstol', 1e-14);
[~, X_] = ode113( @(t,Y) ode_motion_equation(t,Y,muE), tsp, Y_0, options);

figure
Earth = imread('EarthTextureSwap.jpg', 'jpg');
props.FaceColor = 'texturemap';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Earth;
[X_E, Y_E, Z_E] = ellipsoid(0,0,0,RE,RE,RE);
surf(X_E,Y_E,Z_E, props);
set(gca,'color','k')
ax = gca;
ax.GridAlpha = 0.3;
ax.GridColor = [1,1,1];
ax.GridLineStyle = '--';
grid on
hold on
axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
hold on
P(1)=plot3(X_(:,1),X_(:,2),X_(:,3),'linewidth',3);
hold on 
[x, y, z]=ellipsoid(r_0(1), r_0(2), r_0(3), RE/15, RE/15, RE/15);
P(2)=surf(x,y,z,'EdgeColor','none','FaceColor','r');
legend(P,'Original orbit','Satellite initial position','Location','bestoutside','Color','w');
axis equal
axis xy
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Orbit representation');
grid on


return