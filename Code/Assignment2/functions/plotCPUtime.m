function [] = plotCPUtime(n,T,kep0,kP,rP,J2,CD,AtoMratio,rotP)

% plotCPUtime.m - evaluates the CPU time of both Gauss and Cartesian equations
% 		  for the propagation of an orbit, then plots the results, 
% 		  function of the number of time steps for the solution of the 
% 		  ODEs, using linear interpolation.
%
% PROTOTYPE:
%   [] = plotCPUtime(n,T,kep0,kP,rP,J2,CD,AtoMratio,rotP)
%
% INPUT:
%   n         [nx1]   Vector of number of steps of integration        [-]
%   T         [1x1]   Orbit Period                                    [s]
%   kep0      [1x6]   Vector of initial keplerian elements  [km,-,rad,rad,rad,rad]
%   kP        [1x1]   Planetary Constant                              [km^3/s^2]
%   rP        [1x1]   Planetary Radius                                [km]
%   J2        [1x1]   Second Zonal Harmonic                           [-]
%   CD        [1x1]   Drag Coefficient for S/C                        [-]
%   AtoMratio [1x1]   Area to Mass ratio for S/C                      [km^2/kg]
%   rotP      [1x1]   Planet Rotational Rate                          [rad/s]
%
% OUTPUT:
%   Figure with plot of CPU time as a function of number of integration steps
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


tG = zeros(1,length(n));
tC = zeros(1,length(n));
[r0,v0] = kep2car(kep0(1),kep0(2),kep0(3),kep0(4),kep0(5),kep0(6),kP);
condition0 = [r0;v0];

for kcpu = 1:length(n)
    
    t_kcpu = linspace(0,100*T,n(kcpu));
    
    tstartG = cputime;
    [tG_kcpu,kepG_kcpu] = propagatorgauss(t_kcpu,kep0,kP,rP,J2,CD,AtoMratio,rotP);
    tG(kcpu) = cputime - tstartG;
    
    tstartC = cputime;
    [tC_kcpu,rC_kcpu] = propagatorcar(t_kcpu,condition0,kP,rP,J2,CD,AtoMratio,rotP);
    tC(kcpu) = cputime - tstartC;
    
    tG_kcpu = 0;
    kepG_kcpu = 0;
    tC_kcpu = 0;
    rC_kcpu = 0;
end

cG = polyfit(n,tG,1);
cC = polyfit(n,tC,1);
x = n(20):20:n(end);
yG = polyval(cG,x);
yC = polyval(cC,x);

figure
dataG = plot(n(20:end),tG(20:end),'ro');
hold on
grid on
bestfitG = plot(x,yG,'r--');
dataC = plot(n(20:end),tC(20:end),'bo');
bestfitC = plot(x,yC,'b--');
xlabel('$n_{tSTEPS}$', 'interpreter', 'latex', 'FontSize', 12);
ylabel('$t_{CPU}$', 'interpreter', 'latex', 'FontSize',12);
title('CPU Time Comparison');
legend([dataG bestfitG dataC bestfitC], {'Gauss', 'Best-fit line (G)', 'Cartesian', 'Best-fit line (C)'});

return