function [] = T_1Day_10Days_Rep (~,e,i,OM,om,f,muE,om_E,t0,T,m,k,RE,type)

% T_1Day_10Days_Rep.m - plot the repeated ground track for each time duration. 
%
% PROTOTYPE:
%   [] = T_1Day_10Days_Rep (~,e,i,OM,om,f,muE,om_E,t0,T,m,k,RE,type)
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
%   RE                  Earth's radius                        [km]
%   type                String: choosing if plot 1 figure for each one or not
%
% OUTPUT:
%   Figure with the repeated ground track for each time duration. 
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


N = [T, 3600*24, 3600*24*10];       %duration: 1 orbit, 1 day, 10 days
a_new = rep_orbit(om_E,m,k,muE);    %modified samimajor axis

switch type
    case 'all'
        figure

        for l = 1:3
            P=100000;        % number of steps

            tspan = linspace( 0, N(l), P);

            [lon,lat]=groundTrack(a_new,e,i,OM,om,f,muE,om_E,t0,tspan,0,RE); 
            lon=rad2deg(wrapToPi(lon));
            lat=rad2deg(lat);

            for j=2:length(lon)
                if abs(lon(j) - lon(j-1)) > 200
                   lon = [lon(1:j-1) NaN lon(j:end)];
                   lat = [lat(1:j-1) NaN lat(j:end)];
                end
            end

            switch l
                case 1
                    plot(lon(2:end-1),lat(2:end-1),'r','linewidth',5)
                    hold on
                    plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',5)
                    plot(lon(end),lat(end),'sm','linewidth',2,'markersize',11)
                case 2
                    plot(lon(2:end-1),lat(2:end-1),'b','linewidth',3)
                    hold on
                    plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',5)
                    plot(lon(end),lat(end),'sc','linewidth',2,'markersize',9)
                otherwise
                    plot(lon(2:end-1),lat(2:end-1),'y','linewidth',1)
                    hold on
                    plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',5)
                    plot(lon(end),lat(end),'sy','linewidth',2,'markersize',9)
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
        title('Repeating Ground Track');
        grid on
        set(gca,'Layer','top')
        
    case 'one'
        for l = 1:3
            P=100000;        % number of steps

            tspan = linspace( 0, N(l), P);

            [lon,lat]=groundTrack(a_new,e,i,OM,om,f,muE,om_E,t0,tspan,0,RE); 
            lon=rad2deg(wrapToPi(lon));
            lat=rad2deg(lat);

            for j=2:length(lon)
                if abs(lon(j) - lon(j-1)) > 200
                   lon = [lon(1:j-1) NaN lon(j:end)];
                   lat = [lat(1:j-1) NaN lat(j:end)];
                end
            end

            switch l
                case 1
                    figure
                    plot(lon(2:end-1),lat(2:end-1),'y','linewidth',1.5)
                    hold on
                    plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',15)
                    plot(lon(end),lat(end),'sr','linewidth',2,'markersize',20)
                    legend('1 orbit','Start','End','FontSize',12,'location', 'bestoutside')
                    axis([-180 180 -90 90])
                    xticks(-180:30:180)
                    yticks(-90:30:90)
                    xlabel('Longitude [deg]')
                    ylabel('Latitude [deg]')
                    map=imread('EarthTexture.jpg');
                    h = image(xlim,-ylim,map);
                    uistack(h,'bottom')
                    title('Repeating Ground Track');
                    grid on
                    set(gca,'Layer','top')
                case 2
                    figure
                    plot(lon(2:end-1),lat(2:end-1),'y','linewidth',1.5)
                    hold on
                    plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',15)
                    plot(lon(end),lat(end),'sr','linewidth',2,'markersize',20)
                    legend('1 day','Start','End','FontSize',12,'location', 'bestoutside')
                    axis([-180 180 -90 90])
                    xticks(-180:30:180)
                    yticks(-90:30:90)
                    xlabel('Longitude [deg]')
                    ylabel('Latitude [deg]')
                    map=imread('EarthTexture.jpg');
                    h = image(xlim,-ylim,map);
                    uistack(h,'bottom')
                    title('Repeating Ground Track');
                    grid on
                    set(gca,'Layer','top')
                otherwise
                    figure
                    plot(lon(2:end-1),lat(2:end-1),'y','linewidth',1)
                    hold on
                    plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',15)
                    plot(lon(end),lat(end),'sr','linewidth',2,'markersize',20)
                    legend('10 days','Start','End','FontSize',12,'location', 'bestoutside')
                    axis([-180 180 -90 90])
                    xticks(-180:30:180)
                    yticks(-90:30:90)
                    xlabel('Longitude [deg]')
                    ylabel('Latitude [deg]')
                    map=imread('EarthTexture.jpg');
                    h = image(xlim,-ylim,map);
                    uistack(h,'bottom')
                    title('Repeating Ground Track');
                    grid on
                    set(gca,'Layer','top')
            end
        end
end

return