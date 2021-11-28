function [] = T_1Day_10Days_RepSec (~,e,i,OM,om,f,muE,om_E,t0,T,m,k,J2,RE,type)

% T_1Day_10Days_RepSec.m - plot repeated secular ground track for 
%                          each time duration . 
%
% PROTOTYPE:
%   [] = T_1Day_10Days_RepSec (~,e,i,OM,om,f,muE,om_E,t0,T,m,k,J2,RE,type)
%
% INPUT:
%   e                   Eccentricity                          [-]
%   i                   Inclination                           [rad]
%   OM                  RAAN                                  [rad]
%   om                  Argument of periapsis                 [rad]
%   f                   True anomaly                          [rad]      
%   muE                  Planetary gravitational constant      [km^3/s^2]
%   om_E                Rotation velocity of Earth            [rad/s]
%   t0                  Initial time                          [s]
%   T                   Earth's orbit period                  [s]
%   m                   Rotations of the planet
%   k                   revolutions of the satellite
%   J2                  Second zonal harmonic                      
%   RE                  Earth's radius                        [km]
%   type                String: choosing if plot 1 figure for each one or not
%
% OUTPUT:
%   Figure with the repeated secular ground track for each time duration. 
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version



N = [T, 3600*24, 3600*24*10];             %duration: 1 orbit, 1 day, 10 days
a_new = rep_orbit(om_E,m,k,muE);    %modified samimajor axis

switch type
    case 'all'
        figure

        for l = 1:3
            P=100000;        % number of steps

            tspan = linspace( 0, N(l), P);

            a_sec = secular_J2(a_new,e,i,om_E,m,k,muE,J2,RE);

            [lon_s,lat_s]=groundTrack(a_sec,e,i,OM,om,f,muE,om_E,t0,tspan,J2,RE);
            lon_s=rad2deg(wrapToPi(lon_s));
            lat_s=rad2deg(lat_s);

            for j=2:length(lon_s)
                if abs(lon_s(j) - lon_s(j-1)) > 200
                    lon_s = [lon_s(1:j-1) NaN lon_s(j:end)];
                    lat_s = [lat_s(1:j-1) NaN lat_s(j:end)];
                end
            end

            switch l
                case 1
                    plot(lon_s(2:end-1),lat_s(2:end-1),'r','linewidth',5)
                    hold on
                    plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',5)
                    plot(lon_s(end),lat_s(end),'sm','linewidth',2,'markersize',12)
                case 2
                    plot(lon_s(2:end-1),lat_s(2:end-1),'b','linewidth',3)
                    hold on
                    plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',5)
                    plot(lon_s(end),lat_s(end),'sc','linewidth',2,'markersize',9)
                otherwise
                    plot(lon_s(2:end-1),lat_s(2:end-1),'y','linewidth',1)
                    hold on
                    plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',5)
                    plot(lon_s(end),lat_s(end),'sy','linewidth',2,'markersize',9)
            end
        end

        legend('1 orbit','Start','End','1 day','Start','End','10 days','Start','End','location', 'bestoutside')
        axis([-180 180 -90 90])
        xticks(-180:30:180)
        yticks(-90:30:90)
        xlabel('Longitude [deg]')
        ylabel('Latitude [deg]')
        map=imread('EarthTexture.jpg');
        h = image(xlim,-ylim,map);
        uistack(h,'bottom')
        title('Repeating Secular Ground Track');
        grid on
        set(gca,'Layer','top')
        
    case 'one'
        for l = 1:3
            P=100000;        % number of steps

            tspan = linspace( 0, N(l), P);

            a_sec = secular_J2(a_new,e,i,om_E,m,k,muE,J2,RE);

            [lon_s,lat_s]=groundTrack(a_sec,e,i,OM,om,f,muE,om_E,t0,tspan,J2,RE);
            lon_s=rad2deg(wrapToPi(lon_s));
            lat_s=rad2deg(lat_s);

            for j=2:length(lon_s)
                if abs(lon_s(j) - lon_s(j-1)) > 200
                    lon_s = [lon_s(1:j-1) NaN lon_s(j:end)];
                    lat_s = [lat_s(1:j-1) NaN lat_s(j:end)];
                end
            end

            switch l
                case 1
                    figure
                    plot(lon_s(2:end-1),lat_s(2:end-1),'y','linewidth',1.5)
                    hold on
                    plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',15)
                    plot(lon_s(end),lat_s(end),'sr','linewidth',2,'markersize',20)
                    legend('1 orbit','Start','End','FontSize',12,'location', 'bestoutside')
                    axis([-180 180 -90 90])
                    xticks(-180:30:180)
                    yticks(-90:30:90)
                    xlabel('Longitude [deg]')
                    ylabel('Latitude [deg]')
                    map=imread('EarthTexture.jpg');
                    h = image(xlim,-ylim,map);
                    uistack(h,'bottom')
                    title('Repeating Secular Ground Track');
                    grid on
                    set(gca,'Layer','top')
                case 2
                    figure
                    plot(lon_s(2:end-1),lat_s(2:end-1),'y','linewidth',1.5)
                    hold on
                    plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',15)
                    plot(lon_s(end),lat_s(end),'sr','linewidth',2,'markersize',20)
                    legend('1 day','Start','End','FontSize',12,'location', 'bestoutside')
                    axis([-180 180 -90 90])
                    xticks(-180:30:180)
                    yticks(-90:30:90)
                    xlabel('Longitude [deg]')
                    ylabel('Latitude [deg]')
                    map=imread('EarthTexture.jpg');
                    h = image(xlim,-ylim,map);
                    uistack(h,'bottom')
                    title('Repeating Secular Ground Track');
                    grid on
                    set(gca,'Layer','top')
                otherwise
                    figure
                    plot(lon_s(2:end-1),lat_s(2:end-1),'y','linewidth',1)
                    hold on
                    plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',15)
                    plot(lon_s(end),lat_s(end),'sr','linewidth',2,'markersize',20)
                    legend('10 days','Start','End','FontSize',12,'location', 'bestoutside')
                    axis([-180 180 -90 90])
                    xticks(-180:30:180)
                    yticks(-90:30:90)
                    xlabel('Longitude [deg]')
                    ylabel('Latitude [deg]')
                    map=imread('EarthTexture.jpg');
                    h = image(xlim,-ylim,map);
                    uistack(h,'bottom')
                    title('Repeating Secular Ground Track');
                    grid on
                    set(gca,'Layer','top')
            end
        end
end


return
