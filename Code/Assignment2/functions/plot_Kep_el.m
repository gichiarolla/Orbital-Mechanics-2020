function [] = plot_Kep_el (tgauss,kepgauss,tcar,kepcar,T,filter)

% plot_Kep_el.m - plot Cartesian & Gauss elements and filtered Gauss results. 
%
% PROTOTYPE:
%   [] = plot_Kep_el (tgauss,kepgauss,tcar,kepcar,T,filter)
%
% INPUT:
%   tgauss        [nx1]   Gauss time vector                 [s]
%   kepgauss      [6x1]   Gauss keplerian elements          [km,-,rad,rad,rad,rad]
%   tcar          [nx1]   Cartesian time vector             [s]
%   kepcar        [6x1]   Cartesian keplerian elements      [km,-,rad,rad,rad,rad]
%   T                     Earth's orbit period              [s]
%   filter                Input that plot filter instead of cartesian
%
% OUTPUT:
%   Figure with kep els of both methods and real the Gauss filtered
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


fig = figure('WindowState','maximized');   % Open allscreen
pause(0.3);
fig.Position;

if nargin<6
    legend_strings = {'Gauss', 'Cartesian'};
    col_str = [1 0 0; 0 0 1];
    linewidth = [1.1, 0.8];
else
    legend_strings = {'Gauss', 'Filtered'};
    col_str = [0 0 1; 0 1 0];
    linewidth = [1.1, 1.3];
end
    
% 1) Semi-major axis time evolution
subplot(2,3,1)
plot(tgauss./T,kepgauss(:,1),'Color',col_str(1,:),'LineWidth',linewidth(1));
hold on
plot(tcar./T,kepcar(:,1),'Color',col_str(2,:),'LineWidth',linewidth(2));
hold off
title('Semi-Major Axis');
xlabel('Time $[T]$','interpreter','latex');
ylabel('a $[km]$','interpreter','latex','FontSize',14);
legend(legend_strings,'location', 'best')

% 2) Eccentricity time evolution
subplot(2,3,2)
plot(tgauss./T,kepgauss(:,2),'Color',col_str(1,:),'LineWidth',linewidth(1));
hold on
plot(tcar./T,kepcar(:,2),'Color',col_str(2,:),'LineWidth',linewidth(2));
hold off
title('Eccentricity');
xlabel('Time $[T]$','interpreter','latex');
ylabel('e $[-]$','interpreter','latex','FontSize',14);
legend(legend_strings,'location', 'best')

% 3) Inclination time evolution
subplot(2,3,3)
plot(tgauss./T,rad2deg(kepgauss(:,3)),'Color',col_str(1,:),'LineWidth',linewidth(1));
hold on
plot(tcar./T,rad2deg(kepcar(:,3)),'Color',col_str(2,:),'LineWidth',linewidth(2));
title('Inclination');
xlabel('Time $[T]$','interpreter','latex');
ylabel('i $[deg]$','interpreter','latex','FontSize',14);
legend(legend_strings,'location', 'best')

% 4) Right Ascension of the Ascending Node time evolution
subplot(2,3,4)
plot(tgauss./T,rad2deg(kepgauss(:,4)),'Color',col_str(1,:),'LineWidth',linewidth(1));
hold on
plot(tcar./T,rad2deg(kepcar(:,4)),'Color',col_str(2,:),'LineWidth',linewidth(2));
hold off
title('RAAN');
xlabel('Time $[T]$','interpreter','latex');
ylabel('RAAN $[deg]$','interpreter','latex','FontSize',14);
legend(legend_strings,'location', 'best')

% 5) Argument of periapsis time evolution
subplot(2,3,5)
plot(tgauss./T,rad2deg(kepgauss(:,5)),'Color',col_str(1,:),'LineWidth',linewidth(1));
hold on
plot(tcar./T,rad2deg(kepcar(:,5)),'Color',col_str(2,:),'LineWidth',linewidth(2));
hold off
title('Argument of Periapsis');
xlabel('Time $[T]$','interpreter','latex');
ylabel('omega $[deg]$','interpreter','latex','FontSize',14);
legend(legend_strings,'location', 'best')

% 6) True anomaly time evolution
subplot(2,3,6)
plot(tgauss./T,rad2deg(kepgauss(:,6)),'Color',col_str(1,:),'LineWidth',linewidth(1));
hold on
plot(tcar./T,rad2deg(unwrap(kepcar(:,6))),'Color',col_str(2,:),'LineWidth',linewidth(2));
hold off
title('True Anomaly');
xlabel('Time $[T]$','interpreter','latex');
ylabel('theta $[deg]$','interpreter','latex','FontSize',14);
legend(legend_strings,'location', 'best')


return