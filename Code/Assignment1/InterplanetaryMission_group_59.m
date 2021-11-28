clear; clc;

path(pathdef);
addpath(genpath('functions\'));

% FIRST ASSIGNEMENT: Interplanetary Mission
%
% ASSIGNEMENT DESCRIPTION: 
% The PoliMi Space Agency is carrying out a feasibility study for a potential Interplanetary Explorer Mission visiting three planets in the Solar System.
% As part of the mission analysis team you are requested to perform the preliminary mission analysis.
% You have to study the transfer options from the departure planet to the arrival planet, with a powered gravity assist (fly-by) at the intermediate planet,
% and propose a solution based on the mission cost (measured through the total Δ𝑣).
% The departure, flyby, and arrival planets have been decided by the science team. 
% Constraints on earliest departure and latest arrival have also been set by the launch provider, the systems engineering team, and the Agency’s leadership.

%% Planets definition
% Departure planet: Mars(4)
% Fly-by planet:    Venus(2)               
% Target planet:    Mercury(1)

id_dep=4;
id_fb=2;
id_arr=1;

%% Time constrains definition
% Earliest departure date: 2031 Jun 1
% Latest departure date:   2037 Jun 1
%
% Earliest fly-by date:    2032 Feb 1
% Latest fly-by date:      2038 May 1
%
% Earliest arrival date:   2032 Apr 1
% Latest arrival date:     2038 Jun 1

tdep_1=[2031,6,1,0,0,0];
tdep_2=[2037,6,1,0,0,0];

tfb_1=[2032,2,1,0,0,0];
tfb_2=[2038,5,1,0,0,0];

tarr_1=[2032,4,1,0,0,0];
tarr_2=[2038,6,1,0,0,0];

tw_dep=[tdep_1;tdep_2];
tw_fb=[tfb_1;tfb_2];
tw_arr=[tarr_1;tarr_2];

%% Preliminary paramenters definition
muSun=astroConstants(4);         % Gravitational constant of Sun               
mu_fb=astroConstants(10+id_fb);  % Gravitational constant of the fly-by planet
AU=astroConstants(2);            % Astronomical Unit

Rdep=astroConstants(20+id_dep);  % Mean radius of departure planet [km]
Rfb=astroConstants(20+id_fb);    % Mean radius of fly-by planet [km]
Rarr=astroConstants(20+id_arr);  % Mean radius of arrival planet [km]

h_atm=350;  % Atmospheric height of the fly_by planet [km]

ndep=50;  % Number of elements for the departure grid
nfb=50;   % Number of elements for the fly-by grid
narr=50;  % Number of elements for the arrival grid
 
%% Time vectors definitions
tdep_1_mjd2000=date2mjd2000(tw_dep(1,:)); % Conversion from the date in the Gregorian calendar to the modified Julian day 2000 number of the given date
tdep_2_mjd2000=date2mjd2000(tw_dep(2,:));

tfb_1_mjd2000=date2mjd2000(tw_fb(1,:));
tfb_2_mjd2000=date2mjd2000(tw_fb(2,:));

tarr_1_mjd2000=date2mjd2000(tw_arr(1,:));
tarr_2_mjd2000=date2mjd2000(tw_arr(2,:));

tspan_dep=linspace(tdep_1_mjd2000,tdep_2_mjd2000,ndep); % Definition of the time vectors for departure, flyby and arrival
tspan_fb=linspace(tfb_1_mjd2000,tfb_2_mjd2000,nfb);
tspan_arr=linspace(tarr_1_mjd2000,tarr_2_mjd2000,narr);

%% Computation of all ephemeris
r_dep=zeros(3,ndep);
v_dep=zeros(3,ndep);

for k=1:ndep
    [kep,~]=uplanet(tspan_dep(k),id_dep);
    [r_dep(:,k),v_dep(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),muSun);
end

r_fb=zeros(3,nfb);
v_fb=zeros(3,nfb);

for k=1:nfb
    [kep,~]=uplanet(tspan_fb(k),id_fb);
    [r_fb(:,k),v_fb(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),muSun);
end

r_arr=zeros(3,narr);
v_arr=zeros(3,narr);

for k=1:narr
    [kep,~]=uplanet(tspan_dep(k),id_arr);
    [r_arr(:,k),v_arr(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),muSun);
end

%% Computation of all possible solutions
vt1_arc1=zeros(3,ndep,nfb);
vt2_arc1=zeros(3,ndep,nfb);

vt1_arc2=zeros(3,nfb,narr);
vt2_arc2=zeros(3,nfb,narr);

rp=zeros(ndep,nfb,narr);
em=zeros(ndep,nfb,narr);
ep=zeros(ndep,nfb,narr);
am=zeros(ndep,nfb,narr);
ap=zeros(ndep,nfb,narr);
ym=zeros(6,ndep,nfb,narr);
yp=zeros(6,ndep,nfb,narr);


Dv_dep=zeros(ndep,nfb,narr);
Dv_fb=zeros(ndep,nfb,narr);
Dvtot_fb=zeros(ndep,nfb,narr);
Dv_arr=zeros(ndep,nfb,narr);
Dv_tot=zeros(ndep,nfb,narr);

for i=1:ndep
    for j=1:nfb
        for k=1:narr
            Dt_arc1=(tspan_fb(j)-tspan_dep(i))*24*60*60; % Computation of Dt for the transfer arc before fly-by
            Dt_arc2=(tspan_arr(k)-tspan_fb(j))*24*60*60; % Computation of Dt for the transfer arc after fly-by
            
            if Dt_arc1 <= 0  % I cannot arrive before have been departed
                Dv_dep(i,j,k)=NaN;
                Dv_arr(i,j,k)=NaN;
                Dv_fb(i,j,k)=NaN;
                
            elseif Dt_arc2 <= 0 % I cannot arrive before have been departed
                    Dv_dep(i,j,k)=NaN;
                    Dv_arr(i,j,k)=NaN;
                    Dv_fb(i,j,k)=NaN;
            else
                [~,~,~,~,vt1_arc1(:,i,j),vt2_arc1(:,i,j),~,~]=lambertMR(r_dep(:,i),r_fb(:,j),Dt_arc1,muSun,0,0,0,1); % Computation of the transfer arc before fly-by
                Dv_dep(i,j,k)=norm(vt1_arc1(:,i,j)-v_dep(:,i)); 
                
                [~,~,~,~,vt1_arc2(:,j,k),vt2_arc2(:,j,k),~,~]=lambertMR(r_fb(:,j),r_arr(:,k),Dt_arc2,muSun,0,0,0,1); % Computation of the transfer arc after fly-by
                Dv_arr(i,j,k)=norm(v_arr(:,k)-vt2_arc2(:,j,k));
                
                vinf_m=vt2_arc1(:,i,j)-v_fb(:,j);
                vinf_p=vt1_arc2(:,j,k)-v_fb(:,j);
                [~,rp(i,j,k),~,Dv_fb(i,j,k),Dvtot_fb(i,j,k),em(i,j,k),ep(i,j,k),am(i,j,k),ap(i,j,k),ym(:,i,j,k),yp(:,i,j,k)] = Flyby_GA(vinf_m,vinf_p,2); % Computatio of fly-by
            end
            
            Dv_tot(i,j,k)=Dv_dep(i,j,k)+Dv_arr(i,j,k)+Dv_fb(i,j,k);
                   
        end
    end
end

%% Finding the best feasible solution
Dv0=Dv_tot(1,1,1);

for i=1:ndep
    for j=1:nfb
        for k=1:narr
            if Dv_tot(i,j,k) < Dv0 && rp(i,j,k) > Rfb+h_atm
                Dv_tot_min=Dv_tot(i,j,k);
                ii=i;
                jj=j;
                kk=k;
                Dv0=Dv_tot_min;
            else 
                Dv_tot_min=Dv0;
            end
        end
    end
end

x0_gradient=[tspan_dep(ii);tspan_fb(jj);tspan_arr(kk)]; % starting values for checking with fmincon

dep_date=mjd20002date(tspan_dep(ii));
fb_date=mjd20002date(tspan_fb(jj));
arr_date=mjd20002date(tspan_arr(kk));

rp_opt=rp(ii,jj,kk);
h_fb_opt=rp_opt-Rfb;
em_opt=em(ii,jj,kk);
ep_opt=ep(ii,jj,kk);
am_opt=am(ii,jj,kk);
ap_opt=ap(ii,jj,kk);

Dt_arc1=(tspan_fb(jj)-tspan_dep(ii))*24*60*60;
Dt_arc2=(tspan_arr(kk)-tspan_fb(jj))*24*60*60;
Dt_mission=(Dt_arc1+Dt_arc2)/(24*60*60);
Dt_fb = timeofflight_hyerbolic(id_fb,ym(:,ii,jj,kk),yp(:,ii,jj,kk),em_opt,ep_opt,am_opt,ap_opt);

Dv_dep_opt=Dv_dep(ii,jj,kk);
Dv_fb_opt=Dv_fb(ii,jj,kk);
Dvtot_fb_opt=Dvtot_fb(ii,jj,kk);
Dv_arr_opt=Dv_arr(ii,jj,kk);    

fprintf('\nDV_mission = %4.4f km/s\n',Dv_tot_min) % Display of significant results
fprintf('\nDt_mission = %4.2f days\n',Dt_mission)
fprintf('\nDepartue time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n',dep_date)
fprintf('\nFly-by time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n',fb_date)
fprintf('\nArrival time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n',arr_date)
fprintf('\nh_min_fb = %4.4f km\n',h_fb_opt)
fprintf('\nDt_fly-by = %4.2f s\n',Dt_fb)
fprintf('\nDV_dep = %4.4f km/s\n',Dv_dep_opt)
fprintf('\nDV_GA = %4.4f km/s\n',Dv_fb_opt)
fprintf('\nDVtot_fb = %4.4f km/s\n',Dvtot_fb_opt)
fprintf('\nDV_arr = %4.4f km/s\n',Dv_arr_opt)

%% optimization with the gradient base method
[Dv_tot_g,Dt_g,tdep_grad,tfb_grad,tarr_grad] = ThreePlanetsMission_gradient(4,2,1,tw_dep,tw_fb,tw_arr,x0_gradient);

tdep_g=date2mjd2000(tdep_grad);
tfb_g=date2mjd2000(tfb_grad);
tarr_g=date2mjd2000(tarr_grad);

[kep,~]=uplanet(tdep_g,id_dep);
[r_dep_g,v_dep_g]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),muSun);

[kep,~]=uplanet(tfb_g,id_fb);
[r_fb_g,v_fb_g]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),muSun);

[kep,~]=uplanet(tarr_g,id_arr);
[r_arr_g,v_arr_g]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),muSun);

Dt_arc1_g=(tfb_g-tdep_g)*24*60*60;
Dt_arc2_g=(tarr_g-tfb_g)*24*60*60;

[~,~,~,~,vt1_arc1_g,vt2_arc1_g,~,~]=lambertMR(r_dep_g,r_fb_g,Dt_arc1_g,muSun,0,0,0,1);
Dv_dep_g=norm(vt1_arc1_g'-v_dep_g); 
                
[~,~,~,~,vt1_arc2_g,vt2_arc2_g,~,~]=lambertMR(r_fb_g,r_arr_g,Dt_arc2_g,muSun,0,0,0,1);
Dv_arr_g=norm(v_arr_g-vt2_arc2_g');
                
vinf_m_g=vt2_arc1_g'-v_fb_g;
vinf_p_g=vt1_arc2_g'-v_fb_g;
[~,rp_g,~,Dv_fb_g,Dvtot_fb_g,em_g,ep_g,am_g,ap_g,ym_g,yp_g] = Flyby_GA(vinf_m_g,vinf_p_g,id_fb); 

hfb_g=rp_g-Rfb;
Dt_fb_g = timeofflight_hyerbolic(id_fb,ym_g,yp_g,em_g,ep_g,am_g,ap_g);

[ah,eh,ih,OMGh,omgh,fh] = heliocentric_trajectories(r_dep_g',r_fb_g',r_arr_g',vt1_arc1_g,vt2_arc1_g,vt1_arc2_g,vt2_arc2_g,muSun); % carachterization of the heliocentric 
                                                                                                                                 % trajectories
fprintf('\nDV_mission = %4.4f km/s\n',Dv_tot_g) 
fprintf('\nDt_mission = %4.2f days\n',Dt_g)
fprintf('\nDepartue time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n',tdep_grad)
fprintf('\nFly-by time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n',tfb_grad)
fprintf('\nArrival time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n',tarr_grad)
fprintf('\nh_min_fb = %4.4f km\n',hfb_g)
fprintf('\nDt_fly-by = %4.2f s\n',Dt_fb_g)
fprintf('\nDV_dep = %4.4f km/s\n',Dv_dep_g)
fprintf('\nDV_GA = %4.4f km/s\n',Dv_fb_g)
fprintf('\nDVtot_fb = %4.4f km/s\n',Dvtot_fb_g)
fprintf('\nDV_arr = %4.4f km/s\n',Dv_arr_g)

%% Plot trajectories and fly-by
options=odeset('RelTol',1e-13,'AbsTol',1e-14); % tolerance settings for integrations

% SOI of the fly-by planet
SOI=101.7; % Ratio r_SOI/Rp for the fly-by planet
p_lim=SOI+2; % limits values for the fly-by with SOI plot


% Departure planet orbit
[kep_dep,~]=uplanet(tdep_g,id_dep);
a_dep=kep_dep(1);

r0_dep=r_dep_g;
v0_dep=v_dep_g;
y0_dep=[r0_dep;v0_dep];

T_dep=2*pi*sqrt(a_dep^3/muSun);
t_dep=linspace(0,T_dep,1000);

[~,y_dep]=ode113(@(t,y)eq_of_motion(t,y,muSun),t_dep,y0_dep,options);
r_dep_plot=y_dep(:,1:3)/AU;

% Fly_by planet orbit
[kep_fb,~]=uplanet(tfb_g,id_fb);
a_fb=kep_fb(1);

r0_fb=r_fb_g;
v0_fb=v_fb_g;
y0_fb=[r0_fb;v0_fb];

T_fb=2*pi*sqrt(a_fb^3/muSun);
t_fb=linspace(0,T_fb,1000);

[~,y_fb]=ode113(@(t,y)eq_of_motion(t,y,muSun),t_fb,y0_fb,options);
r_fb_plot=y_fb(:,1:3)/AU;

% Arrival planet orbit
[kep_arr,~]=uplanet(tarr_g,id_arr);
a_arr=kep_arr(1);

r0_arr=r_arr_g;
v0_arr=v_arr_g;
y0_arr=[r0_arr;v0_arr];

T_arr=2*pi*sqrt(a_arr^3/muSun);
t_arr=linspace(0,T_arr,1000);

[~,y_arr]=ode113(@(t,y)eq_of_motion(t,y,muSun),t_arr,y0_arr,options);
r_arr_plot=y_arr(:,1:3)/AU;

% Transfer arc 1
vt0_arc1=vt1_arc1_g';
y0_arc1=[r0_dep;vt0_arc1];

tspan_arc1=linspace(0,Dt_arc1_g,1000);

[~,y_arc1]=ode113(@(t,y)eq_of_motion(t,y,muSun),tspan_arc1,y0_arc1,options);
r_arc1=y_arc1(:,1:3)/AU;

% Transfer arc 2
vt0_arc2=vt1_arc2_g';
y0_arc2=[r0_fb;vt0_arc2];

tspan_arc2=linspace(0,Dt_arc2_g,1000);

[~,y_arc2]=ode113(@(t,y)eq_of_motion(t,y,muSun),tspan_arc2,y0_arc2,options);
r_arc2=y_arc2(:,1:3)/AU;

% Fly-by incoming arc and asymptote
y0_in=ym_g;

rp_vers=y0_in(1:3)/rp_g;

tspan_in=linspace(0,-SOI*1000,1000);

[~,y_in]=ode45(@(t,y)eq_of_motion(t,y,mu_fb),tspan_in,y0_in,options);

r_in=y_in(:,1:3)/Rfb;

cm=(rp_g-am_g)*rp_vers;
bm=-am_g*sqrt(em_g^2-1);
fm=@(x) cm(2)+bm/(-am_g) *(x-cm(1));
xlim_m=cm(2)*am_g/bm+cm(1);
xm=linspace(-5*Rfb,xlim_m,1000);

% Fly-by outcoming arc and asymptote
y0_out=yp_g;

tspan_out=linspace(0,SOI*1000,1000);

[~,y_out]=ode45(@(t,y)eq_of_motion(t,y,mu_fb),tspan_out,y0_out,options);

r_out=y_out(:,1:3)/Rfb;

cp=(rp_g-ap_g)*rp_vers;
bp=-ap_g*sqrt(ep_g^2-1);
fp=@(x) cp(2)-bp/(-ap_g) *(x-cp(1));
xlim_p=-cp(2)*ap_g/bp+cp(1);
xp=linspace(-5*Rfb,xlim_p,1000);

% Departure, flyby and arrival positions for the 3 planets
[kep_dep2,~]=uplanet(tfb_g,id_dep);
[pos_dep2,~]=kep2car(kep_dep2(1),kep_dep2(2),kep_dep2(3),kep_dep2(4),kep_dep2(5),kep_dep2(6),muSun);
pos_dep2=pos_dep2/AU;

[kep_dep3,~]=uplanet(tarr_g,id_dep);
[pos_dep3,~]=kep2car(kep_dep3(1),kep_dep3(2),kep_dep3(3),kep_dep3(4),kep_dep3(5),kep_dep3(6),muSun);
pos_dep3=pos_dep3/AU;

[kep_fb1,~]=uplanet(tdep_g,id_fb);
[pos_fb1,~]=kep2car(kep_fb1(1),kep_fb1(2),kep_fb1(3),kep_fb1(4),kep_fb1(5),kep_fb1(6),muSun);
pos_fb1=pos_fb1/AU;

[kep_fb3,~]=uplanet(tarr_g,id_fb);
[pos_fb3,~]=kep2car(kep_fb3(1),kep_fb3(2),kep_fb3(3),kep_fb3(4),kep_fb3(5),kep_fb3(6),muSun);
pos_fb3=pos_fb3/AU;

[kep_arr1,~]=uplanet(tdep_g,id_arr);
[pos_arr1,~]=kep2car(kep_arr1(1),kep_arr1(2),kep_arr1(3),kep_arr1(4),kep_arr1(5),kep_arr1(6),muSun);
pos_arr1=pos_arr1/AU;

[kep_arr2,~]=uplanet(tfb_g,id_arr);
[pos_arr2,~]=kep2car(kep_arr2(1),kep_arr2(2),kep_arr2(3),kep_arr2(4),kep_arr2(5),kep_arr2(6),muSun);
pos_arr2=pos_arr2/AU;

% Celestial corps 
[xx,yy,zz]=sphere(100);

x_dep=xx*Rdep*1000/AU;
y_dep=yy*Rdep*1000/AU;
z_dep=zz*Rdep*1000/AU;

x_fb=xx*Rfb*1000/AU;
y_fb=yy*Rfb*1000/AU;
z_fb=zz*Rfb*1000/AU;

x_arr=xx*Rarr*1000/AU;
y_arr=yy*Rarr*1000/AU;
z_arr=zz*Rarr*1000/AU;

x_Sun=xx*6955000/AU;
y_Sun=yy*6955000/AU;
z_Sun=zz*6955000/AU;

% Total transfer plot
figure()
plot3(r_dep_plot(:,1),r_dep_plot(:,2),r_dep_plot(:,3),'Color','#A2142F','LineStyle','--','LineWidth',1.5)
hold on

plot3(r_fb_plot(:,1),r_fb_plot(:,2),r_fb_plot(:,3),'Color','#EDB120','LineStyle','--','LineWidth',1.5)
hold on

plot3(r_arr_plot(:,1),r_arr_plot(:,2),r_arr_plot(:,3),'Color','#77AC30','LineStyle','--','LineWidth',1.5)
hold on

plot3(r_arc1(:,1),r_arc1(:,2),r_arc1(:,3),'Color','#7E2F8E','LineStyle','-','LineWidth',1.5)
hold on

plot3(r_arc2(:,1),r_arc2(:,2),r_arc2(:,3),'Color','#D95319','LineStyle','-','LineWidth',1.5)
hold on

plot3(pos_arr1(1),pos_arr1(2),pos_arr1(3),'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k')
hold on

plot3(pos_dep2(1),pos_dep2(2),pos_dep2(3),'o','MarkerSize',6,'MarkerFaceColor','y','MarkerEdgeColor','k')
hold on

plot3(pos_fb3(1),pos_fb3(2),pos_fb3(3),'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on

plot3(pos_dep3(1),pos_dep3(2),pos_dep3(3),'o','MarkerSize',6,'MarkerFaceColor','r','MarkerEdgeColor','k')
hold on

plot3(pos_fb1(1),pos_fb1(2),pos_fb1(3),'o','MarkerSize',6,'MarkerFaceColor','g','MarkerEdgeColor','k')
hold on

plot3(pos_arr2(1),pos_arr2(2),pos_arr2(3),'o','MarkerSize',6,'MarkerFaceColor','y','MarkerEdgeColor','k')
hold on

map_dep=imread('MarsTexture.jpg');
map_dep=imrotate(map_dep,180);
warp(x_dep+r_dep_plot(1,1),y_dep+r_dep_plot(1,2),z_dep+r_dep_plot(1,3),map_dep);

map_fb=imread('VenusTexture.jpg');
map_fb=imrotate(map_fb,180);
warp(x_fb+r_fb_plot(1,1),y_fb+r_fb_plot(1,2),z_fb+r_fb_plot(1,3),map_fb);

map_arr=imread('MercuryTexture.jpg');
map_arr=imrotate(map_arr,180);
warp(x_arr+r_arr_plot(1,1),y_arr+r_arr_plot(1,2),z_arr+r_arr_plot(1,3),map_arr);

map_Sun=imread('SunTexture.jpg');
map_Sun=imrotate(map_Sun,180);
warp(x_Sun,y_Sun,z_Sun,map_Sun);


set(gca,'YDir','normal');
xlabel('x [AU]')
ylabel('y [AU]')
zlabel('z [AU]')
legend('Mars orbit','Venus orbit','Mercury orbit','Transfer arc before fly-by','Transfer arc after fly-by','Position at Departure','Position at fly-by',...
       'Position at arrival','FontSize',12,'LineWidth',1.2)
title('Total transfer trajectory','FontSize',18,'BackgroundColor','w','EdgeColor','k','Interpreter','latex','LineWidth',1.2)
grid on

% Fly-by plot
figure()
plot3(r_in(:,1),r_in(:,2),r_in(:,3),'LineWidth',1)
hold on

plot3(r_out(:,1),r_out(:,2),r_out(:,3),'LineWidth',1)
hold on

viscircles([0,0],SOI,'Color','#4DBEEE','LineWidth',0.5);

map_fb=imread('VenusTexture.jpg');
map_fb=imrotate(map_fb,180);
warp(xx,yy,zz,map_fb);

set(gca,'YDir','normal');
xlabel('x [R_{Venus}]')
ylabel('y [R_{Venus}]')
zlabel('z [R_{Venus}]')
legend('Incoming hyperbola','Outcoming hyperbola','FontSize',12,'LineWidth',1.2);
title('Fly-by within SOI','FontSize',18,'BackgroundColor','w','EdgeColor','k','Interpreter','latex','LineWidth',1.2)
grid on
xlim([-p_lim p_lim])
ylim([-p_lim p_lim])

figure()
plot3(r_in(:,1),r_in(:,2),r_in(:,3),'LineWidth',2)
hold on

plot3(r_out(:,1),r_out(:,2),r_out(:,3),'LineWidth',2)
hold on

plot(xm/Rfb,fm(xm)/Rfb,'LineWidth',1.1)
hold on

plot(xp/Rfb,fp(xp)/Rfb,'LineWidth',1.1)
hold on

map_fb=imread('VenusTexture.jpg');
map_fb=imrotate(map_fb,180);
warp(xx,yy,zz,map_fb);

set(gca,'YDir','normal');
xlabel('x [R_{Venus}]')
ylabel('y [R_{Venus}]')
zlabel('z [R_{Venus}]')
legend('Incoming hyperbola','Outcoming hyperbola','Incoming asymptote','Outcoming asymptote','FontSize',12,'LineWidth',1.2);
title('Fly-by zoom','FontSize',18,'BackgroundColor','w','EdgeColor','k','Interpreter','latex','LineWidth',1.2)
grid on
xlim([-2 2])
ylim([-4 4])

%% Video
video=0; % video = 1 to run the video section, video = 0 otherwise

if video == 1
    file_name='Assignement1_video_final';

% Step 1: generate data
dep=tdep_g;
arr=tarr_g;

n=1000; % number of steps
n1=round(Dt_arc1_g*n/(Dt_arc1_g+Dt_arc2_g));
n2=n-n1;
t=linspace(0,Dt_arc1_g+Dt_arc2_g,n);
t1=linspace(0,Dt_arc1_g,n1);
t2=linspace(0,Dt_arc2_g,n2);

date_mjd2000=linspace(dep,arr,n);
date=zeros(n,6);

options=odeset('RelTol',1e-13,'AbsTol',1e-14);

[~,y_dep_movie]=ode113(@(t,y)eq_of_motion(t,y,muSun),t,y0_dep,options); % Mars
r_dep_movie=y_dep_movie(:,1:3)/AU;

[kep_Ea,~]=uplanet(dep,3); % Earth
[r0_Ea_movie,v0_Ea_movie]=kep2car(kep_Ea(1),kep_Ea(2),kep_Ea(3),kep_Ea(4),kep_Ea(5),kep_Ea(6),muSun);
y0_Ea_movie=[r0_Ea_movie;v0_Ea_movie];
[~,y_Ea_movie]=ode113(@(t,y)eq_of_motion(t,y,muSun),t,y0_Ea_movie,options);
r_Ea_movie=y_Ea_movie(:,1:3)/AU;
a_Ea=kep_Ea(1);
T_Ea=2*pi*sqrt(a_Ea^3/muSun);
tspan_Ea=linspace(0,T_Ea,1000);
[~,y_Ea]=ode113(@(t,y)eq_of_motion(t,y,muSun),tspan_Ea,y0_Ea_movie,options);
r_Ea=y_Ea(:,1:3)/AU;
x_Ea=xx*6371000/AU;
y_Ea=yy*6371000/AU;
z_Ea=zz*6371000/AU;

[~,y_fb_movie1]=ode113(@(t,y)eq_of_motion(t,y,muSun),-t1,y0_fb,options); % Venus
r_fb_movie1=y_fb_movie1(:,1:3)/AU;
[~,y_fb_movie2]=ode113(@(t,y)eq_of_motion(t,y,muSun),t2,y0_fb,options); 
r_fb_movie2=y_fb_movie2(:,1:3)/AU;
r_fb_movie=[r_fb_movie1;r_fb_movie2];

[~,y_arr_movie]=ode113(@(t,y)eq_of_motion(t,y,muSun),-t,y0_arr,options); % Mercury
r_arr_movie=y_arr_movie(:,1:3)/AU;

[~,y_arc1_movie]=ode113(@(t,y)eq_of_motion(t,y,muSun),t1,y0_arc1,options); % Transfer arcs
r_arc1_movie=y_arc1_movie(:,1:3)/AU;
[~,y_arc2_movie]=ode113(@(t,y)eq_of_motion(t,y,muSun),t2,y0_arc2,options);
r_arc2_movie=y_arc2_movie(:,1:3)/AU;
r_arc_movie=[r_arc1_movie;r_arc2_movie];

% Step 2: draw the scenario
fig=figure('WindowState','maximized');
% fig.Position

for k=1:n
    % Time data for the title
    date(k,:)=mjd20002date(date_mjd2000(k));
    
    % Wipe the state clean in order to plot on a blank figure
    clf
    
    % Plot the entire orbits
    plot3(r_dep_plot(:,1),r_dep_plot(:,2),r_dep_plot(:,3),'Color','#A2142F','LineStyle','--','LineWidth',1.5)
    
    hold on
    plot3(r_Ea(:,1),r_Ea(:,2),r_Ea(:,3),'Color','#0072BD','LineStyle','--','LineWidth',1.5)
    
    hold on
    plot3(r_fb_plot(:,1),r_fb_plot(:,2),r_fb_plot(:,3),'Color','#EDB120','LineStyle','--','LineWidth',1.5)
    
    hold on
    plot3(r_arr_plot(:,1),r_arr_plot(:,2),r_arr_plot(:,3),'Color','#77AC30','LineStyle','--','LineWidth',1.5)
    
    % Plot the current location of the palnets & s/c
    hold on
    if k <= n1 
        plot3(r_arc_movie(1:k,1),r_arc_movie(1:k,2),r_arc_movie(1:k,3),'Color','#7E2F8E','LineStyle','-','LineWidth',1.5)
    else
        plot3(r_arc_movie(1:k,1),r_arc_movie(1:k,2),r_arc_movie(1:k,3),'Color','#D95319','LineStyle','-','LineWidth',1.5)
        plot3(r_arc_movie(1:n1,1),r_arc_movie(1:n1,2),r_arc_movie(1:n1,3),'Color','#7E2F8E','LineStyle','-','LineWidth',1.5)
    end
    
    hold on
    plot3(r_arc_movie(k,1),r_arc_movie(k,2),r_arc_movie(k,3),'Color','k','Marker','o','MarkerFaceColor','r','MarkerSize',5)
    
    map_dep=imread('MarsTexture.jpg');
    map_dep=imrotate(map_dep,180);
    warp(x_dep+r_dep_movie(k,1),y_dep+r_dep_movie(k,2),z_dep+r_dep_movie(k,3),map_dep);
    
    map_Ea=imread('EarthTexture.jpg');
    map_Ea=imrotate(map_Ea,180);
    warp(x_Ea+r_Ea_movie(k,1),y_Ea+r_Ea_movie(k,2),z_Ea+r_Ea_movie(k,3),map_Ea);
    
    if k <= n1 
        map_fb=imread('VenusTexture.jpg');
        map_fb=imrotate(map_fb,180);
        warp(x_fb+r_fb_movie(n1+1-k,1),y_fb+r_fb_movie(n1+1-k,2),z_fb+r_fb_movie(n1+1-k,3),map_fb);
    else
        map_fb=imread('VenusTexture.jpg');
        map_fb=imrotate(map_fb,180);
        warp(x_fb+r_fb_movie(k,1),y_fb+r_fb_movie(k,2),z_fb+r_fb_movie(k,3),map_fb);
    end
       
    map_arr=imread('MercuryTexture.jpg');
    map_arr=imrotate(map_arr,180);
    warp(x_arr+r_arr_movie(n+1-k,1),y_arr+r_arr_movie(n+1-k,2),z_arr+r_arr_movie(n+1-k,3),map_arr);
    
    map_Sun=imread('SunTexture.jpg');
    map_Sun=imrotate(map_Sun,180);
    warp(x_Sun,y_Sun,z_Sun,map_Sun);
 
    % Decorate the plot
    grid on
    set(gca,'YDir','normal');
    xlabel('x [AU]')
    ylabel('y [AU]')
    zlabel('z [AU]')
    
    if k <= n1 
        legend('Mars orbit','Venus orbit','Earth orbit','Mercury orbit','Transfer arc before fly-by','s/c','Location','northwest','FontSize',12,'LineWidth',1.2);
    else
         legend('Mars orbit','Venus orbit','Earth orbit','Mercury orbit','Transfer arc after fly-by','Transfer arc before fly-by','s/c',...
                                                                                                           'Location','northwest','FontSize',12,'LineWidth',1.2)
    end
    
    title(['Date:  ',num2str(date(k,1:3))],'FontSize',20,'Position',[2.2,0.3,0],'BackgroundColor','w','EdgeColor','k','Interpreter','latex','LineWidth',1.2)
    xlim([-2 2])
    ylim([-2 2])
    zlim([-0.1 0.1])
    view(-50,85)
    
    % Force matlab to draw the image at this point
%     movie_vector=zeros(1,n);
    movie_vector(k)=getframe(fig,[60 40 1400 710]);   
end   

% Step 5: save the movie
MOVIE=VideoWriter(file_name,'MPEG-4');

open(MOVIE);
writeVideo(MOVIE,movie_vector);
close(MOVIE);
end
