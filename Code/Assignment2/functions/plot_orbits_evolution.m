function [] = plot_orbits_evolution (RE,kepgauss,ntotOrb,nOrb,tspan,muE,string)

% plot_orbits_evolution.m - plot all the orbits, for both Gauss end Cartesian. 
%
% PROTOTYPE:
%   [] = plot_orbits_evolution (RE,kepgauss,ntotOrb,nOrb,tspan,muE)
%
% INPUT:
%   RE                  Earth' radius                           [km]
%   kepgauss    [6x1]   Gauss keplerian elements          [km,-,rad,rad,rad,rad]
%   ntotOrb             Number of orbits                        [-]
%   nOrb                Step for orbits to plot                 [-]
%   tspan       [nx1]   Time Vector                             [s]
%   muE         [1x1]   Planetary gravitational constant     [km^3/s^2]
%   string              String for J2 or Drag legend
%
% OUTPUT:
%   Figure with all the orbits of both methods 
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


% Plot of orbit
figure
Earth = imread('EarthTextureSwap.jpg', 'jpg');
props.FaceColor = 'texturemap';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = Earth;
[X_E, Y_E, Z_E] = ellipsoid(0,0,0,RE,RE,RE);
Earth = surf(X_E,Y_E,Z_E, props);
set(gca,'color','k')
ax = gca;
ax.GridAlpha = 0.3;
ax.GridColor = [1,1,1];
ax.GridLineStyle = '--';
grid on
hold on
axis equal
hcb = colorbar;
caxis([0 tspan(end)]/(24*3600))
set(get(hcb,'Title'),'String',string, 'Interpreter', 'latex');

fullRevolution = [0:pi/50:2*pi];
X = [];
Y = [];
Z = [];
nstep = ntotOrb/nOrb;
Step = floor(length(tspan)/nstep);
cmap = parula(nstep);

for kPlot1 = 1:nstep
    kPlot2 = Step*(kPlot1-1)+1;
    for k1Orbit = 1:length(fullRevolution)
        [rPlot,~] = kep2car(kepgauss(kPlot2,1),kepgauss(kPlot2,2),kepgauss(kPlot2,3),kepgauss(kPlot2,4),kepgauss(kPlot2,5),fullRevolution(k1Orbit),muE);
        X = [X;rPlot(1)];
        Y = [Y;rPlot(2)];
        Z = [Z;rPlot(3)];
    end
    orbitcolor = cmap(kPlot1,:);
    plot3(X,Y,Z,'Color',orbitcolor);
    X = []; Y = []; Z = [];
end

xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');

return