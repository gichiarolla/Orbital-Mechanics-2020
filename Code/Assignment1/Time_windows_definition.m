clear; clc;

% Sript to make the porkchop plots used to select the time windows

%% First arc: Mars-Venus

%  Departure planet:                Mars
%  Target planet:                   Venus
%  Earliest departure requirement:  2031 June 1
%  Later departure requirement:     2071 June 1
%  Earliest arrival requirement:    2031 June 1
%  Later arrival requirement:       2071 June 1

% Time limits definition
t1_start=[2031,6,1,0,0,0];
t1_end=[2071,6,1,0,0,0];
t2_start=[2031,6,1,0,0,0];
t2_end=[2071,6,1,0,0,0];

t1_start_mjd2000=date2mjd2000(t1_start);
t1_end_mjd2000=date2mjd2000(t1_end);
t2_start_mjd2000=date2mjd2000(t2_start);
t2_end_mjd2000=date2mjd2000(t2_end);

% Get departure ephemeris
ndep=8000;
t1span_mjd2000=linspace(t1_start_mjd2000,t1_end_mjd2000,ndep);

r1_vec=zeros(3,ndep);
v1_vec=zeros(3,ndep);

for k=1:ndep
    [kep,ksun]=uplanet(t1span_mjd2000(k),4);
    [r1_vec(:,k),v1_vec(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),ksun);
end

% Get arrival ephemeris
narr=8000;
t2span_mjd2000=linspace(t2_start_mjd2000,t2_end_mjd2000,narr);

r2_vec=zeros(3,narr);
v2_vec=zeros(3,narr);

for k=1:narr
    [kep,ksun]=uplanet(t2span_mjd2000(k),2);
    [r2_vec(:,k),v2_vec(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),ksun);
end

% Compute all possible transfers
vt1_vec=zeros(3,ndep,narr);
vt2_vec=zeros(3,ndep,narr);
Dvtot_vec=zeros(ndep,narr);

for k1=1:ndep
    for k2=1:narr
        Dt=(t2span_mjd2000(k2)-t1span_mjd2000(k1))*24*60*60;
        [~,~,~,~,vt1_vec(:,k1,k2),vt2_vec(:,k1,k2),~,~]=lambertMR(r1_vec(:,k1),r2_vec(:,k2),Dt,ksun,0,0,0,1);
        Dvtot_vec(k1,k2)=norm(vt1_vec(:,k1,k2)-v1_vec(:,k1))+norm(v2_vec(:,k2)-vt2_vec(:,k1,k2));
    end
end 

% Minimum cost trajectory
Dvmin=min(min(Dvtot_vec));
[k1,k2]=find(Dvtot_vec==Dvmin);
t1=mjd20002date(t1span_mjd2000(k1));
t2=mjd20002date(t2span_mjd2000(k2));
r1=r1_vec(:,k1);
v1=v1_vec(:,k1);
r2=r2_vec(:,k2);
v2=v2_vec(:,k2);
vt1=vt1_vec(:,k1,k2);
vt2=vt2_vec(:,k1,k2);

% Date in matlab format
 t1_plot=zeros(ndep,1);
for k=1:ndep
    t1_plot(k)=datenum(mjd20002date(t1span_mjd2000(k)));
end

t2_plot=zeros(narr,1);
for k=1:narr
    t2_plot(k)=datenum(mjd20002date(t2span_mjd2000(k)));
end

% Draw the porkchop plot
figure()
[C,h]=contour(t1_plot,t2_plot, Dvtot_vec',floor(Dvmin)+(0:1:10));
clabel(C,h,(6:1:10))
caxis(floor(Dvmin)+[0 5])
caxis('manual')
xtickangle(45)
ytickangle(45)
datetick('x','yyyy mmm dd','keeplimits')
datetick('y','yyyy mmm dd','keeplimits')
xlabel('Departure date')
ylabel('Arrival date')
title('Mars-Venus')
grid on

cb=colorbar;
cb.Title.String='$\Delta v$ [km/s]';
cb.Title.Interpreter='latex';
cb.Title.FontSize=15;

ax=gca;
ax.FontSize=12;

hold on
plot(t1_plot(k1),t2_plot(k2),'o','MarkerSize',3,'MarkerFaceColor','b')

fprintf('\nDeparture time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f \n\n',t1)
fprintf('\nArrival time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n',t2)
fprintf('\nDV: %4.4f \n\n',Dvmin)

%% First arc: Mars-Venus detail

%  Departure planet:                Mars
%  Target planet:                   Venus
%  Earliest departure requirement:  2044 December 1
%  Later departure requirement:     2046 Jenuary 1
%  Earliest arrival requirement:    2046 February 1
%  Later arrival requirement:       2046 July 1

% Time limits definition
t1_start=[2044,12,1,12,0,0];
t1_end=[2046,1,1,12,0,0];
t2_start=[2046,2,1,12,0,0];
t2_end=[2046,7,1,12,0,0];

t1_start_mjd2000=date2mjd2000(t1_start);
t1_end_mjd2000=date2mjd2000(t1_end);
t2_start_mjd2000=date2mjd2000(t2_start);
t2_end_mjd2000=date2mjd2000(t2_end);

% Get departure ephemeris
ndep=1000;
t1span_mjd2000=linspace(t1_start_mjd2000,t1_end_mjd2000,ndep);

r1_vec=zeros(3,ndep);
v1_vec=zeros(3,ndep);

for k=1:ndep
    [kep,ksun]=uplanet(t1span_mjd2000(k),4);
    [r1_vec(:,k),v1_vec(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),ksun);
end

% Get arrival ephemeris
narr=1000;
t2span_mjd2000=linspace(t2_start_mjd2000,t2_end_mjd2000,narr);

r2_vec=zeros(3,narr);
v2_vec=zeros(3,narr);

for k=1:narr
    [kep,ksun]=uplanet(t2span_mjd2000(k),2);
    [r2_vec(:,k),v2_vec(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),ksun);
end

% Compute all possible transfers
vt1_vec=zeros(3,ndep,narr);
vt2_vec=zeros(3,ndep,narr);
Dvtot_vec=zeros(ndep,narr);

for k1=1:ndep
    for k2=1:narr
        Dt=(t2span_mjd2000(k2)-t1span_mjd2000(k1))*24*60*60;
        [~,~,~,~,vt1_vec(:,k1,k2),vt2_vec(:,k1,k2),~,~]=lambertMR(r1_vec(:,k1),r2_vec(:,k2),Dt,ksun,0,0,0,1);
        Dvtot_vec(k1,k2)=norm(vt1_vec(:,k1,k2)-v1_vec(:,k1))+norm(v2_vec(:,k2)-vt2_vec(:,k1,k2));
    end
end 

% Minimum cost trajectory
Dvmin=min(min(Dvtot_vec));
[k1,k2]=find(Dvtot_vec==Dvmin);
t1=mjd20002date(t1span_mjd2000(k1));
t2=mjd20002date(t2span_mjd2000(k2));
r1=r1_vec(:,k1);
v1=v1_vec(:,k1);
r2=r2_vec(:,k2);
v2=v2_vec(:,k2);
vt1=vt1_vec(:,k1,k2);
vt2=vt2_vec(:,k1,k2);

% Date in matlab format
 t1_plot=zeros(ndep,1);
for k=1:ndep
    t1_plot(k)=datenum(mjd20002date(t1span_mjd2000(k)));
end

t2_plot=zeros(narr,1);
for k=1:narr
    t2_plot(k)=datenum(mjd20002date(t2span_mjd2000(k)));
end

% Draw the porkchop plot
figure()
[C,h]=contour(t1_plot,t2_plot, Dvtot_vec',floor(Dvmin)+(0:1:10));
clabel(C,h,(10:1:15))
caxis(floor(Dvmin)+[0 5])
caxis('manual')
xtickangle(45)
ytickangle(45)
datetick('x','yyyy mmm dd','keeplimits')
datetick('y','yyyy mmm dd','keeplimits')
xlabel('Departure date')
ylabel('Arrival date')
title('Mars-Venus')
grid on

cb=colorbar;
cb.Title.String='$\Delta v$ [km/s]';
cb.Title.Interpreter='latex';
cb.Title.FontSize=15;

ax=gca;
ax.FontSize=12;

hold on
plot(t1_plot(k1),t2_plot(k2),'o','MarkerSize',3,'MarkerFaceColor','b')

fprintf('\nDeparture time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f \n\n',t1)
fprintf('\nArrival time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n',t2)
fprintf('\nDV: %4.4f \n\n',Dvmin)

%% Second arc: Venus-Mercury

%  Departure planet:                Venus
%  Target planet:                   Mercury
%  Earliest departure requirement:  2031 June 1
%  Later departure requirement:     2071 June 1
%  Earliest arrival requirement:    2031 June 1
%  Later arrival requirement:       2071 June 1

% Time limits definition
t1_start=[2031,6,1,0,0,0];
t1_end=[2071,6,1,0,0,0];
t2_start=[2031,6,1,0,0,0];
t2_end=[2071,6,1,0,0,0];

t1_start_mjd2000=date2mjd2000(t1_start);
t1_end_mjd2000=date2mjd2000(t1_end);
t2_start_mjd2000=date2mjd2000(t2_start);
t2_end_mjd2000=date2mjd2000(t2_end);

% Get departure ephemeris
ndep=8000;
t1span_mjd2000=linspace(t1_start_mjd2000,t1_end_mjd2000,ndep);

r1_vec=zeros(3,ndep);
v1_vec=zeros(3,ndep);

for k=1:ndep
    [kep,ksun]=uplanet(t1span_mjd2000(k),2);
    [r1_vec(:,k),v1_vec(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),ksun);
end

% Get arrival ephemeris
narr=8000;
t2span_mjd2000=linspace(t2_start_mjd2000,t2_end_mjd2000,narr);

r2_vec=zeros(3,narr);
v2_vec=zeros(3,narr);

for k=1:narr
    [kep,ksun]=uplanet(t2span_mjd2000(k),1);
    [r2_vec(:,k),v2_vec(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),ksun);
end

% Compute all possible transfers
vt1_vec=zeros(3,ndep,narr);
vt2_vec=zeros(3,ndep,narr);
Dvtot_vec=zeros(ndep,narr);

for k1=1:ndep
    for k2=1:narr
        Dt=(t2span_mjd2000(k2)-t1span_mjd2000(k1))*24*60*60;
        [~,~,~,~,vt1_vec(:,k1,k2),vt2_vec(:,k1,k2),~,~]=lambertMR(r1_vec(:,k1),r2_vec(:,k2),Dt,ksun,0,0,0,1);
        Dvtot_vec(k1,k2)=norm(vt1_vec(:,k1,k2)-v1_vec(:,k1))+norm(v2_vec(:,k2)-vt2_vec(:,k1,k2));
    end
end 

% Minimum cost trajectory
Dvmin=min(min(Dvtot_vec));
[k1,k2]=find(Dvtot_vec==Dvmin);
t1=mjd20002date(t1span_mjd2000(k1));
t2=mjd20002date(t2span_mjd2000(k2));
r1=r1_vec(:,k1);
v1=v1_vec(:,k1);
r2=r2_vec(:,k2);
v2=v2_vec(:,k2);
vt1=vt1_vec(:,k1,k2);
vt2=vt2_vec(:,k1,k2);

% Date in matlab format
 t1_plot=zeros(ndep,1);
for k=1:ndep
    t1_plot(k)=datenum(mjd20002date(t1span_mjd2000(k)));
end

t2_plot=zeros(narr,1);
for k=1:narr
    t2_plot(k)=datenum(mjd20002date(t2span_mjd2000(k)));
end

% Draw the porkchop plot
figure()
[C,h]=contour(t1_plot,t2_plot, Dvtot_vec',floor(Dvmin)+(0:1:10));
clabel(C,h,(6:1:10))
caxis(floor(Dvmin)+[0 5])
caxis('manual')
xtickangle(45)
ytickangle(45)
datetick('x','yyyy mmm dd','keeplimits')
datetick('y','yyyy mmm dd','keeplimits')
xlabel('Departure date')
ylabel('Arrival date')
title('Venus-Mercury')
grid on

cb=colorbar;
cb.Title.String='$\Delta v$ [km/s]';
cb.Title.Interpreter='latex';
cb.Title.FontSize=15;

ax=gca;
ax.FontSize=12;

hold on
plot(t1_plot(k1),t2_plot(k2),'o','MarkerSize',3,'MarkerFaceColor','b')

fprintf('\nDeparture time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f \n\n',t1)
fprintf('\nArrival time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n',t2)
fprintf('\nDV: %4.4f \n\n',Dvmin)

%% Second arc: Venus-Mercury detail

%  Departure planet:                Venus
%  Target planet:                   Mercury
%  Earliest departure requirement:  2032 July 1
%  Later departure requirement:     2032 October 1
%  Earliest arrival requirement:    2032 October 1
%  Later arrival requirement:       2032 December 1

% Time limits definition
t1_start=[2032,7,1,12,0,0];
t1_end=[2032,10,1,12,0,0];
t2_start=[2032,10,1,12,0,0];
t2_end=[2032,12,1,12,0,0];

t1_start_mjd2000=date2mjd2000(t1_start);
t1_end_mjd2000=date2mjd2000(t1_end);
t2_start_mjd2000=date2mjd2000(t2_start);
t2_end_mjd2000=date2mjd2000(t2_end);

% Get departure ephemeris
ndep=1000;
t1span_mjd2000=linspace(t1_start_mjd2000,t1_end_mjd2000,ndep);

r1_vec=zeros(3,ndep);
v1_vec=zeros(3,ndep);

for k=1:ndep
    [kep,ksun]=uplanet(t1span_mjd2000(k),2);
    [r1_vec(:,k),v1_vec(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),ksun);
end

% Get arrival ephemeris
narr=1000;
t2span_mjd2000=linspace(t2_start_mjd2000,t2_end_mjd2000,narr);

r2_vec=zeros(3,narr);
v2_vec=zeros(3,narr);

for k=1:narr
    [kep,ksun]=uplanet(t2span_mjd2000(k),1);
    [r2_vec(:,k),v2_vec(:,k)]=kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),ksun);
end

% Compute all possible transfers
vt1_vec=zeros(3,ndep,narr);
vt2_vec=zeros(3,ndep,narr);
Dvtot_vec=zeros(ndep,narr);

for k1=1:ndep
    for k2=1:narr
        Dt=(t2span_mjd2000(k2)-t1span_mjd2000(k1))*24*60*60;
        [~,~,~,~,vt1_vec(:,k1,k2),vt2_vec(:,k1,k2),~,~]=lambertMR(r1_vec(:,k1),r2_vec(:,k2),Dt,ksun,0,0,0,1);
        Dvtot_vec(k1,k2)=norm(vt1_vec(:,k1,k2)-v1_vec(:,k1))+norm(v2_vec(:,k2)-vt2_vec(:,k1,k2));
    end
end 

% Minimum cost trajectory
Dvmin=min(min(Dvtot_vec));
[k1,k2]=find(Dvtot_vec==Dvmin);
t1=mjd20002date(t1span_mjd2000(k1));
t2=mjd20002date(t2span_mjd2000(k2));
r1=r1_vec(:,k1);
v1=v1_vec(:,k1);
r2=r2_vec(:,k2);
v2=v2_vec(:,k2);
vt1=vt1_vec(:,k1,k2);
vt2=vt2_vec(:,k1,k2);

% Date in matlab format
 t1_plot=zeros(ndep,1);
for k=1:ndep
    t1_plot(k)=datenum(mjd20002date(t1span_mjd2000(k)));
end

t2_plot=zeros(narr,1);
for k=1:narr
    t2_plot(k)=datenum(mjd20002date(t2span_mjd2000(k)));
end

% Draw the porkchop plot
figure()
[C,h]=contour(t1_plot,t2_plot, Dvtot_vec',floor(Dvmin)+(0:1:10));
clabel(C,h,(13:1:18))
caxis(floor(Dvmin)+[0 5])
caxis('manual')
xtickangle(45)
ytickangle(45)
datetick('x','yyyy mmm dd','keeplimits')
datetick('y','yyyy mmm dd','keeplimits')
xlabel('Departure date')
ylabel('Arrival date')
title('Venus-Mercury')
grid on

cb=colorbar;
cb.Title.String='$\Delta v$ [km/s]';
cb.Title.Interpreter='latex';
cb.Title.FontSize=15;

ax=gca;
ax.FontSize=12;

hold on
plot(t1_plot(k1),t2_plot(k2),'o','MarkerSize',3,'MarkerFaceColor','b')

fprintf('\nDeparture time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f \n\n',t1)
fprintf('\nArrival time: %4.0f %4.0f %4.0f %4.0f %4.0f %4.0f\n',t2)
fprintf('\nDV: %4.4f \n\n',Dvmin)



