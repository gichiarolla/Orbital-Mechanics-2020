function [] = plot_SC_el (dateMJD,KepMethod,KepSC,sizet,T,tgauss_SC,CartORGauss)

% plot_SC_el.m - plot Cartesian or Gauss elements and real S/C datas. 
%
% PROTOTYPE:
%   [] = plot_SC_el (dateMJD,KepMethod,KepSC,sizet,T,tgauss_SC,CartORGauss)
%
% INPUT:
%   dateMJD     [1xn]   Time vector in MJD                [days]
%   KepMethod   [6x1]   Method keplerian elements         [km,-,rad,rad,rad,rad]
%   KepSC       [6x1]   Real S/C keplerian elements       [km,-,rad,rad,rad,rad]
%   sizet       	Time vector dimension             [s]
%   T                   Period of the S/C orbit           [s]
%   tgauss_SC    [nx1]  Time vector from Gauss simulation [s]
%   CartORGauss         String to choose the method
%
% OUTPUT:
%   Figure with Cartesian or Gauss elements and real S/C datas 
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


%label for Legend (no needed to edit when changing type of method)
switch CartORGauss
    case 'Gauss'
        legend_string = {'Gauss Simulation', 'Gauss Filtered', 'Ephemerides'};
    case 'Cart'
        legend_string = {'Cart Simulation', 'Cart Filtered', 'Ephemerides'};
end


tplot = linspace(dateMJD(1),dateMJD(end),sizet);

Tfilter = 4*T;
nwindow = nearest(Tfilter/(sum(diff(tgauss_SC))/(numel(tgauss_SC)-1)));
kepGfiltered = movmean(KepMethod,nwindow,1);

fig = figure('WindowState','maximized');   % Open allscreen
pause(0.3);
fig.Position;
sgtitle('Keplerian Parameters Time Evolution for ARASE (ERG)');

% 1) Semi-major axis time evolution
subplot(2,3,1)
plot(tplot,KepMethod(:,1),'r','LineWidth',0.5);
hold on
plot(tplot,kepGfiltered(:,1),'g','LineWidth',0.5);
hold on
plot(dateMJD,KepSC(:,1),'b','LineWidth',1.1);
hold off
legend(legend_string,'location', 'best')
title('Semi-Major Axis');
xlabel('Time $[MJD]$','interpreter','latex');
ylabel('a $[km]$','interpreter','latex','FontSize',14);
grid on

% 2) Eccentricity time evolution
subplot(2,3,2)
plot(tplot,KepMethod(:,2),'r','LineWidth',0.5);
hold on
plot(tplot,kepGfiltered(:,2),'g','LineWidth',0.5);
hold on
plot(dateMJD,KepSC(:,2),'b','LineWidth',1.1);
hold off
legend(legend_string,'location', 'best')
title('Eccentricity');
xlabel('Time $[MJD]$','interpreter','latex');
ylabel('e $[-]$','interpreter','latex','FontSize',14);
grid on

% 3) Inclination time evolution
subplot(2,3,3)
plot(tplot,rad2deg(KepMethod(:,3)),'r','LineWidth',0.5);
hold on
plot(tplot,rad2deg(kepGfiltered(:,3)),'g','LineWidth',0.8);
hold on
plot(dateMJD,rad2deg(KepSC(:,3)),'b','LineWidth',1.1);
hold off
legend(legend_string,'location', 'best')
title('Inclination');
xlabel('Time $[MJD]$','interpreter','latex');
ylabel('i $[deg]$','interpreter','latex','FontSize',14);
grid on

% 4) Right Ascension of the Ascending Node time evolution
subplot(2,3,4)
plot(tplot,rad2deg(KepMethod(:,4)),'r','LineWidth',0.8);
hold on
plot(tplot,rad2deg(kepGfiltered(:,4)),'g','LineWidth',0.7);
hold on
plot(dateMJD,rad2deg(unwrap(KepSC(:,4))),'b','LineWidth',1.1);
hold off
legend(legend_string,'location', 'best')
title('Right Ascension of the Ascending Node');
xlabel('Time $[MJD]$','interpreter','latex');
ylabel('RAAN $[deg]$','interpreter','latex','FontSize',14);
grid on

% 5) Argument of periapsis time evolution
subplot(2,3,5)
plot(tplot,rad2deg(KepMethod(:,5)),'r','LineWidth',0.8);
hold on
plot(tplot,rad2deg(kepGfiltered(:,5)),'g','LineWidth',0.7);
hold on
plot(dateMJD,rad2deg(unwrap(KepSC(:,5))),'b','LineWidth',1.1);
hold off
legend(legend_string,'location', 'best')
title('Argument of Perigee');
xlabel('Time $[MJD]$','interpreter','latex');
ylabel('omega $[deg]$','interpreter','latex','FontSize',14);
grid on


return