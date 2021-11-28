function [] = orbitsMovie (R_Cart,R_CartF,R_Gauss,R_GaussF,tspan,RE,T,om_Edeg,M)

% orbitsMovie.m - create the movie of orbit rapresentation, 
%                 for both Gauss end Cartesian.
%
% PROTOTYPE:
%   [] = orbitsMovie (R_Cart,R_CartF,R_Gauss,R_GaussF,tspan,RE,T,om_Edeg,M)
%
% INPUT:
%   R_Cart    [3x1]   Radius of 1st Cartesian orbit           [km]
%   R_CartF   [3x1]   Radius of last Cartesian orbit          [km]
%   R_Gauss   [3x1]   Radius of 1st Gauss orbit               [km]
%   R_GaussF  [3x1]   Radius of last Gauss orbit              [km]
%   tspan     [nx1]   Time vector                             [s]
%   RE                Earth' radius                           [km]
%   T                 Period of the orbit                     [s]
%   om_Edeg           Rotation velocity of Earth              [deg/s]
%   M                 Number of points of 1 orbit             [-]
%
% OUTPUT:
%   Figures (movie) with orbits' evolution during time
%   File saved as AVI
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version



deg=360/T;                                              % deg for rotation to follow s/c

for i=1:3                                               % 3 draw (cart, kep, both)
    
    switch i
        
        case 1                                          % Cartesian orbit
            figC = figure('WindowState','maximized');   % Open allscreen
            pause(0.3);
            figC.Position;
            
            for k=1:M
                clf
                
                %extract data at current time
                tk = tspan (k);
                xk = R_Cart(k,1);
                yk = R_Cart(k,2);
                zk = R_Cart(k,3);
                xkF = R_CartF(k,1);
                ykF = R_CartF(k,2);
                zkF = R_CartF(k,3);

                %earth orbit
                Earth = imread('EarthTextureSwap.jpg', 'jpg');
                props.FaceColor = 'texturemap';
                props.EdgeColor = 'none';
                props.FaceLighting = 'phong';
                props.Cdata = Earth;
                [X_E, Y_E, Z_E] = ellipsoid(0,0,0,RE,RE,RE);
                S=surf(X_E,Y_E,Z_E, props);
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
                title(['Cartesian orbits at T = ', num2str(k/M),' period'])
                ax.TitleHorizontalAlignment = 'left';
                rotate(S,[0 0 1],om_Edeg*tk)
                 
                %plot current location of satellite
                P(1)=plot3(xk, yk, zk, 'bo', 'LineWidth', 2.5, 'MarkerSize', RE/900);
                hold on
                P(2)=plot3(R_Cart(1:k,1),R_Cart(1:k,2),R_Cart(1:k,3),'b', 'LineWidth', 2);
                hold on
                P(3)=plot3(xkF, ykF, zkF, 'o', 'Color', [0 0.89 1], 'LineWidth', 2.5, 'MarkerSize', RE/900);
                hold on
                P(4)=plot3(R_CartF(1:k,1),R_CartF(1:k,2),R_CartF(1:k,3),'LineStyle',':','Color','c', 'LineWidth', 2.5);

                view(230+deg*tk,15)         %view with rotation
                legend(P,'S/C position 1stT','First Orbit','S/C position lastT','Last Orbit','Location',[0.801099827737443,0.856983387235045,0.076562498509884,0.072170299161632],'Color','w');
                
%                 drawnow
                movieVec_C(k) = getframe(figC, [395 35 1320 900]);
            end
            
            %save movie
            myWriter_C = VideoWriter('Cart orbits');
            myWriter_C.FrameRate = 20;
            open(myWriter_C);
            writeVideo(myWriter_C, movieVec_C);
            close(myWriter_C);
            
        case 2
            figG = figure('WindowState','maximized');
            pause(0.3);
            figG.Position;

            for k=1:M
                clf
                
                %extract data at current time
                tk = tspan (k);
                xk = R_Gauss(k,1);
                yk = R_Gauss(k,2);
                zk = R_Gauss(k,3);
                xkF = R_GaussF(k,1);
                ykF = R_GaussF(k,2);
                zkF = R_GaussF(k,3);
                
                %earth orbit
                Earth = imread('EarthTextureSwap.jpg', 'jpg');
                props.FaceColor = 'texturemap';
                props.EdgeColor = 'none';
                props.FaceLighting = 'phong';
                props.Cdata = Earth;
                [X_E, Y_E, Z_E] = ellipsoid(0,0,0,RE,RE,RE);
                S=surf(X_E,Y_E,Z_E, props);
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
                title(['Gauss orbits at T = ', num2str(k/M),' period'])
                ax.TitleHorizontalAlignment = 'left';
                rotate(S,[0 0 1],om_Edeg*tk)
                 
                %plot current location of satellite
                P(1)=plot3(xk, yk, zk, 'ro', 'LineWidth', 2.5, 'MarkerSize', RE/900);
                hold on
                P(2)=plot3(R_Gauss(1:k,1),R_Gauss(1:k,2),R_Gauss(1:k,3),'r', 'LineWidth', 2);
                hold on
                P(3)=plot3(xkF, ykF, zkF, 'o', 'Color',[0.850980392156863 0.325490196078431 0.0980392156862745], 'LineWidth', 2.5, 'MarkerSize', RE/900);
                hold on
                P(4)=plot3(R_GaussF(1:k,1),R_GaussF(1:k,2),R_GaussF(1:k,3),'LineStyle',':','Color',[0.850980392156863 0.325490196078431 0.0980392156862745], 'LineWidth', 2.5);
                
                view(230+deg*tk,15)         %view with rotation
                legend(P,'S/C position 1stT','First Orbit','S/C position lastT','Last Orbit','Location',[0.801099827737443,0.856983387235045,0.076562498509884,0.072170299161632],'Color','w');

%                 drawnow
                movieVec_G(k) = getframe(figG, [395 35 1320 900]);
            end
            %save movie
            myWriter_G = VideoWriter('Gauss orbits');
            myWriter_G.FrameRate = 20;
            open(myWriter_G);
            writeVideo(myWriter_G, movieVec_G);
            close(myWriter_G);

        case 3
            figGC = figure('WindowState','maximized');
            pause(0.3);
            figGC.Position;

            for k=1:M
                clf
                
                %extract data at current time
                tk = tspan (k);
                xkG = R_Gauss(k,1);
                ykG = R_Gauss(k,2);
                zkG = R_Gauss(k,3);
                xkC = R_Cart(k,1);
                ykC = R_Cart(k,2);
                zkC = R_Cart(k,3);
                xkCF = R_CartF(k,1);
                ykCF = R_CartF(k,2);
                zkCF = R_CartF(k,3);
                xkGF = R_GaussF(k,1);
                ykGF = R_GaussF(k,2);
                zkGF = R_GaussF(k,3);
                %earth orbit
                Earth = imread('EarthTextureSwap.jpg', 'jpg');
                props.FaceColor = 'texturemap';
                props.EdgeColor = 'none';
                props.FaceLighting = 'phong';
                props.Cdata = Earth;
                [X_E, Y_E, Z_E] = ellipsoid(0,0,0,RE,RE,RE);
                S=surf(X_E,Y_E,Z_E, props);
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
                title(['Gauss & Cart orbits at T = ', num2str(k/M),' period'])
                ax.TitleHorizontalAlignment = 'left';
                rotate(S,[0 0 1],om_Edeg*tk)
                 
                %plot current location of satellite
                P(1)=plot3(xkC, ykC, zkC, 'bo', 'LineWidth', 2.5, 'MarkerSize', RE/900);
                hold on
                P(2)=plot3(R_Cart(1:k,1),R_Cart(1:k,2),R_Cart(1:k,3),'b', 'LineWidth', 2);
                hold on
                P(3)=plot3(xkCF, ykCF, zkCF, 'o', 'Color', [0 0.89 1], 'LineWidth', 2.5, 'MarkerSize', RE/900);
                hold on
                P(4)=plot3(R_CartF(1:k,1),R_CartF(1:k,2),R_CartF(1:k,3),'Color','c', 'LineWidth', 2);
                hold on
                P(5)=plot3(xkG, ykG, zkG, 'ro', 'LineWidth', 2.5, 'MarkerSize', RE/450);
                hold on
                P(6)=plot3(R_Gauss(1:k,1),R_Gauss(1:k,2),R_Gauss(1:k,3),'LineStyle',':','Color','r', 'LineWidth', 2.5);
                hold on
                P(7)=plot3(xkGF, ykGF, zkGF, 'o', 'Color',[0.850980392156863 0.325490196078431 0.0980392156862745], 'LineWidth', 2.5, 'MarkerSize', RE/450);
                hold on
                P(8)=plot3(R_GaussF(1:k,1),R_GaussF(1:k,2),R_GaussF(1:k,3),'LineStyle',':','Color',[0.850980392156863 0.325490196078431 0.0980392156862745], 'LineWidth', 2.5);
                
                view(230+deg*tk,15)         %view with rotation
                legend(P,'S/C_C position 1stT','Cart 1stT','S/C_C position lastT','Cart lastT',...
                    'S/C_G position 1stT','Gauss 1stT','S/C_G position lastT','Gauss lastT','Location',[0.801099827737443,0.856983387235045,0.076562498509884,0.072170299161632],'Color','w');
                                
%                 drawnow
                movieVec_GC(k) = getframe(figGC, [395 35 1350 925]);

            end
            %save movie
            myWriter_GC = VideoWriter('Gauss & Cart orbits');
            myWriter_GC.FrameRate = 20;
            open(myWriter_GC);
            writeVideo(myWriter_GC, movieVec_GC);
            close(myWriter_GC);
    end
            
end
  
    
 
    
    