function [errR_a,errR_e,errR_i,errR_RAAN,errR_omega,errR_TA] = plot_comparison_Kep_el (a,T,kepgauss,kepcar,tspan)

% plot_comparison_Kep_el.m - compare Cartesian & Gauss methods 
%                            by plotting the difference between kep els.
%
% PROTOTYPE:
%   [] = plot_comparison_Kep_el (a,T,kepgauss,kepcar,tspan)
%
% INPUT:
%   a                 Semi-major axis                   [km]
%   T                 Period of the orbit               [s]
%   kepgauss  [6x1]   Gauss keplerian elements          [km,-,rad,rad,rad,rad]
%   kepcar    [6x1]   Cartesian keplerian elements      [km,-,rad,rad,rad,rad]
%   tspan     [nx1]   Time vector                       [s]
%
% OUTPUT:
%   errR_a = norm(err_a,Inf);		Norm infinte error of semi-major axis
%   errR_e = norm(err_e,Inf);		Norm infinte error of eccentricity
%   errR_i = norm(err_i,Inf);		Norm infinte error of inclination
%   errR_RAAN = norm(err_RAAN,Inf);	Norm infinte error of RAAN
%   errR_omega = norm(err_omega,Inf);	Norm infinte error of omega
%   errR_TA = norm(err_TA,Inf);		Norm infinte error of true anomaly
%   Figure with kep els of both methods	
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


err_a = (kepcar(:,1)-kepgauss(:,1))/a;
err_e = (kepcar(:,2)-kepgauss(:,2));
err_i = (kepcar(:,3)-kepgauss(:,3))/(2*pi);
err_RAAN = (kepcar(:,4)-kepgauss(:,4))/(2*pi);
err_omega = (kepcar(:,5)-kepgauss(:,5))/(2*pi);
err_TA = (unwrap(kepcar(:,6))-kepgauss(:,6))./kepgauss(:,6);

errR_a = norm(err_a,Inf);
errR_e = norm(err_e,Inf);
errR_i = norm(err_i,Inf);
errR_RAAN = norm(err_RAAN,Inf);
errR_omega = norm(err_omega,Inf);
errR_TA = norm(err_TA,Inf);
 
fig = figure('WindowState','maximized');   % Open allscreen
pause(0.3);
fig.Position;

subplot(2,3,1)
semilogy(tspan./T,abs(err_a),'k');
grid on
title('Semi-Major Axis');
xlabel('Periods $[T]$','interpreter','latex');
ylabel('${|a_{Cartesian}-a_{Gauss}| / a_{0}}$','interpreter','latex','FontSize',14);

subplot(2,3,2)
semilogy(tspan./T,abs(err_e),'k');
grid on
title('Eccentricity');
xlabel('Periods $[T]$','interpreter','latex');
ylabel('${|e_{Cartesian} - e_{Gauss}|}$','interpreter','latex','FontSize',14);

subplot(2,3,3)
semilogy(tspan./T,abs(err_i),'k');
grid on
title('Inclination');
xlabel('Periods $[T]$','interpreter','latex');
ylabel('${|i_{Cartesian} - i_{Gauss}| / 2\pi}$','interpreter','latex','FontSize',14);

subplot(2,3,4)
semilogy(tspan./T,abs(err_RAAN),'k');
grid on
title('RAAN');
xlabel('Periods $[T]$','interpreter','latex');
ylabel('${|RAAN_{Cartesian} - RAAN_{Gauss}| / 2\pi}$','interpreter','latex','FontSize',14);

subplot(2,3,5)
semilogy(tspan./T,abs(err_omega),'k');
grid on
title('Argument of Periapsis');
xlabel('Periods $[T]$','interpreter','latex');
ylabel('${|omega_{Cartesian} - omega_{Gauss}| / 2\pi}$','interpreter','latex','FontSize',14);

subplot(2,3,6)
semilogy(tspan./T,abs(err_TA),'k');
grid on
title('True Anomaly');
xlabel('Periods $[T]$','interpreter','latex');
ylabel('${|theta_{Cartesian}-theta_{Gauss}| / theta_{Gauss}}$','interpreter','latex','FontSize',14);


return