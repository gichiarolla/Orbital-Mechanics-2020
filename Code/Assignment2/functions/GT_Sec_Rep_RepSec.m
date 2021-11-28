function [] = GT_Sec_Rep_RepSec (a,e,i,OM,om,f,muE,om_E,t0,T,m,k,J2,RE)

% GT_Sec_Rep_RepSec.m - plot all types of ground track for each time duration.
%
% PROTOTYPE:
%   [] = GT_Sec_Rep_RepSec (a,e,i,OM,om,f,muE,om_E,t0,T,m,k,J2,RE)
%
% INPUT:
%   a                   Semi-major axis                       [km]
%   e                   Eccentricity                          [-]
%   i                   Inclination                           [rad]
%   OM                  RAAN                                  [rad]
%   om                  Argument of periapsis                 [rad]
%   f                   True anomaly                          [rad]      
%   mu                  Planetary gravitational constant      [km^3/s^2]
%   om_E                Rotation velocity of Earth            [rad/s]
%   t0                  Initial time                          [s]
%   T                   Earth's orbit period                  [s]
%   m                   Rotations of the planet
%   k                   revolutions of the satellite
%   J2                  Second zonal harmonic                      
%   RE                  Earth's radius                        [km]
%
% OUTPUT:
%   Figure with ground tracks' rapresentation
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


N = [T, 3600*24, 3600*24*10]; %duration: 1 orbit, 1 day, 10 days

for l = 1:3
    P=100000;        % number of steps
    
    tspan = linspace( 0, N(l), P);
    
    switch l
        case 1
            %GT
            [lon,lat]=groundTrack(a,e,i,OM,om,f,muE,om_E,t0,tspan,0,RE); 
            lon=rad2deg(wrapToPi(lon));
            lat=rad2deg(lat);

            for j=2:length(lon)
                if abs(lon(j) - lon(j-1)) > 200
                   lon = [lon(1:j-1) NaN lon(j:end)];
                   lat = [lat(1:j-1) NaN lat(j:end)];
                end
            end
            figure
            plot(lon(2:end-1),lat(2:end-1),'r','linewidth',2.5)
            hold on
            plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',7)
            plot(lon(end),lat(end),'sm','linewidth',2,'markersize',9)
            hold on
            
            %secular
            a_sec = secular_J2(a,e,i,om_E,m,k,muE,J2,RE);

            [lon_s,lat_s]=groundTrack(a_sec,e,i,OM,om,f,muE,om_E,t0,tspan,J2,RE);
            lon_s=rad2deg(wrapToPi(lon_s));
            lat_s=rad2deg(lat_s);

            for j=2:length(lon_s)
                if abs(lon_s(j) - lon_s(j-1)) > 200
                    lon_s = [lon_s(1:j-1) NaN lon_s(j:end)];
                    lat_s = [lat_s(1:j-1) NaN lat_s(j:end)];
                end
            end
            plot(lon_s(2:end-1),lat_s(2:end-1),'k','linewidth',2.5)
            hold on
            plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon_s(end),lat_s(end),'sk','linewidth',3,'markersize',14)
            hold on

            %repeated
            a_new = rep_orbit(om_E,m,k,muE);    %modified samimajor axis

            [lon,lat]=groundTrack(a_new,e,i,OM,om,f,muE,om_E,t0,tspan,0,RE); 
            lon=rad2deg(wrapToPi(lon));
            lat=rad2deg(lat);

            for j=2:length(lon)
                if abs(lon(j) - lon(j-1)) > 200
                   lon = [lon(1:j-1) NaN lon(j:end)];
                   lat = [lat(1:j-1) NaN lat(j:end)];
                end
            end

            plot(lon(2:end-1),lat(2:end-1),'b','linewidth',1.5)
            hold on
            plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon(end),lat(end),'sc','linewidth',2,'markersize',11)
            hold on
            
            %repeated secular
            a_new = rep_orbit(om_E,m,k,muE);
            
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

            plot(lon_s(2:end-1),lat_s(2:end-1),'y','linewidth',1)
            hold on
            plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon_s(end),lat_s(end),'sy','linewidth',2,'markersize',8)
            
                        
            legend('GroundTrack','Start','End','Secular','Start','End',...
                'Repeated','Start','End','Repeated Secular','Start','End','Location', 'bestoutside')
            axis([-180 180 -90 90])
            xticks(-180:30:180)
            yticks(-90:30:90)
            xlabel('Longitude [deg]')
            ylabel('Latitude [deg]')
            map=imread('EarthTexture.jpg');
            h = image(xlim,-ylim,map);
            uistack(h,'bottom')
            title('1 T GroundTrack');
            grid on
            set(gca,'Layer','top')
            
            
        case 2
            %GT
            [lon,lat]=groundTrack(a,e,i,OM,om,f,muE,om_E,t0,tspan,0,RE); 
            lon=rad2deg(wrapToPi(lon));
            lat=rad2deg(lat);

            for j=2:length(lon)
                if abs(lon(j) - lon(j-1)) > 200
                   lon = [lon(1:j-1) NaN lon(j:end)];
                   lat = [lat(1:j-1) NaN lat(j:end)];
                end
            end
            figure
            plot(lon(2:end-1),lat(2:end-1),'r','linewidth',2.5)
            hold on
            plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon(end),lat(end),'sm','linewidth',2,'markersize',11)
            hold on
            
            %secular
            a_sec = secular_J2(a,e,i,om_E,m,k,muE,J2,RE);

            [lon_s,lat_s]=groundTrack(a_sec,e,i,OM,om,f,muE,om_E,t0,tspan,J2,RE);
            lon_s=rad2deg(wrapToPi(lon_s));
            lat_s=rad2deg(lat_s);

            for j=2:length(lon_s)
                if abs(lon_s(j) - lon_s(j-1)) > 200
                    lon_s = [lon_s(1:j-1) NaN lon_s(j:end)];
                    lat_s = [lat_s(1:j-1) NaN lat_s(j:end)];
                end
            end
            plot(lon_s(2:end-1),lat_s(2:end-1),'k','linewidth',2.5)
            hold on
            plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon_s(end),lat_s(end),'sk','linewidth',3,'markersize',14)
            hold on

            %repeated
            a_new = rep_orbit(om_E,m,k,muE);    %modified samimajor axis

            [lon,lat]=groundTrack(a_new,e,i,OM,om,f,muE,om_E,t0,tspan,0,RE); 
            lon=rad2deg(wrapToPi(lon));
            lat=rad2deg(lat);

            for j=2:length(lon)
                if abs(lon(j) - lon(j-1)) > 200
                   lon = [lon(1:j-1) NaN lon(j:end)];
                   lat = [lat(1:j-1) NaN lat(j:end)];
                end
            end

            plot(lon(2:end-1),lat(2:end-1),'b','linewidth',1.5)
            hold on
            plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon(end),lat(end),'sc','linewidth',2,'markersize',11)
            hold on
            
            %repeated secular
            a_new = rep_orbit(om_E,m,k,muE);
            
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

            plot(lon_s(2:end-1),lat_s(2:end-1),'y','linewidth',1)
            hold on
            plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon_s(end),lat_s(end),'sy','linewidth',2,'markersize',8)
            
                        
            legend('GroundTrack','Start','End','Secular','Start','End',...
                'Repeated','Start','End','Repeated Secular','Start','End','location', 'bestoutside')
            axis([-180 180 -90 90])
            xticks(-180:30:180)
            yticks(-90:30:90)
            xlabel('Longitude [deg]')
            ylabel('Latitude [deg]')
            map=imread('EarthTexture.jpg');
            h = image(xlim,-ylim,map);
            uistack(h,'bottom')
            title('1 Day GroundTrack');
            grid on
            set(gca,'Layer','top')
            
            
        otherwise
           %GT
            [lon,lat]=groundTrack(a,e,i,OM,om,f,muE,om_E,t0,tspan,0,RE); 
            lon=rad2deg(wrapToPi(lon));
            lat=rad2deg(lat);

            for j=2:length(lon)
                if abs(lon(j) - lon(j-1)) > 200
                   lon = [lon(1:j-1) NaN lon(j:end)];
                   lat = [lat(1:j-1) NaN lat(j:end)];
                end
            end
            figure
            plot(lon(2:end-1),lat(2:end-1),'r','linewidth',2.5)
            hold on
            plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon(end),lat(end),'sm','linewidth',2,'markersize',14)
            hold on
            
            %secular
            a_sec = secular_J2(a,e,i,om_E,m,k,muE,J2,RE);

            [lon_s,lat_s]=groundTrack(a_sec,e,i,OM,om,f,muE,om_E,t0,tspan,J2,RE);
            lon_s=rad2deg(wrapToPi(lon_s));
            lat_s=rad2deg(lat_s);

            for j=2:length(lon_s)
                if abs(lon_s(j) - lon_s(j-1)) > 200
                    lon_s = [lon_s(1:j-1) NaN lon_s(j:end)];
                    lat_s = [lat_s(1:j-1) NaN lat_s(j:end)];
                end
            end
            plot(lon_s(2:end-1),lat_s(2:end-1),'k','linewidth',2.5)
            hold on
            plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon_s(end),lat_s(end),'sk','linewidth',3,'markersize',14)
            hold on

            %repeated
            a_new = rep_orbit(om_E,m,k,muE);    %modified samimajor axis

            [lon,lat]=groundTrack(a_new,e,i,OM,om,f,muE,om_E,t0,tspan,0,RE); 
            lon=rad2deg(wrapToPi(lon));
            lat=rad2deg(lat);

            for j=2:length(lon)
                if abs(lon(j) - lon(j-1)) > 200
                   lon = [lon(1:j-1) NaN lon(j:end)];
                   lat = [lat(1:j-1) NaN lat(j:end)];
                end
            end

            plot(lon(2:end-1),lat(2:end-1),'b','linewidth',1.5)
            hold on
            plot(lon(1),lat(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon(end),lat(end),'sc','linewidth',2,'markersize',11)
            hold on
            
            %repeated secular
            a_new = rep_orbit(om_E,m,k,muE);
            
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

            plot(lon_s(2:end-1),lat_s(2:end-1),'y','linewidth',1)
            hold on
            plot(lon_s(1),lat_s(1),'oc','linewidth',1.5,'markersize',5)
            plot(lon_s(end),lat_s(end),'sy','linewidth',2,'markersize',8)
            
                        
            legend('GroundTrack','Start','End','Secular','Start','End',...
                'Repeated','Start','End','Repeated Secular','Start','End','location', 'bestoutside')
            axis([-180 180 -90 90])
            xticks(-180:30:180)
            yticks(-90:30:90)
            xlabel('Longitude [deg]')
            ylabel('Latitude [deg]')
            map=imread('EarthTexture.jpg');
            h = image(xlim,-ylim,map);
            uistack(h,'bottom')
            title('10 Days GroundTrack');
            grid on
            set(gca,'Layer','top')
    end
end

return