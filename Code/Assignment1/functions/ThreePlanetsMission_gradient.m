function [Dv_tot,Dt,tdep,tfb,tarr] = ThreePlanetsMission_gradient(id_dep,id_fb,id_arr,tw_dep,tw_fb,tw_arr,x0)

% ThreePlanetsMission_gradient - Optimiozation of the Best departure, fly-by, 
%                                arrival time and total cost of the mission 
%                                obtained with the grid search method,
%                                throughout a gradient-based optimizer 
%
% PROTOTYPE:
%   [Dv_check,Dt_check,tdep_check,tfb_check,tarr_check] = ThreePlanetsMission_gradient(id_dep,id_fb,id_arr,tw_dep,tw_fb,tw_arr,x0_check)
%
%  INPUT:
%    id_dep[1]    Integer number identifying the departure planet. 
%    id_fb[1]     Integer number identifying the fly-by planet.
%    id_arr[1]    Integer number identifying the arrival planet. 
%                   Mercury(1),Venus(2),Earth(3),Mars(4),Jupiter(5),Saturn(6),Uranus(7),Neptune(8).
%    tw_dep[2x6]  Matrix containing the minimum and maximum departure dates in the Gregorian calendar,
%                 as 6-elements vectors [year, month, day, hour, minute, second].
%    tw_fb[2x6]   Matrix containing the minimum and maximum departure dates in the Gregorian calendar,
%                 as 6-elements vectors [year, month, day, hour, minute, second].
%    tw_arr[2x6]  Matrix containing the minimum and maximum departure dates in the Gregorian calendar,
%                 as 6-elements vectors [year, month, day, hour, minute, second].
%                   For dates before 1582, the resulting date components are valid only in the Gregorian proleptic calendar. 
%                   This is based on the Gregorian calendar but extended to cover dates before its introduction. 
%                   Date must be after 12:00 noon, 24 November -4713.
%    x0[3x1]      Vector containing the mjd2000 dates for the optimal solution obtained by the grid search method.
%
%  OUTPUT:
%    Dt[1]        Duration of the whole mission [days].
%    Dv_tot[1]    Total cost of the mission [km/s]
%    tdep[6]      Departure date [year, month, day, hour, minute, second].
%    tfb[6]       Fly-by date [year, month, day, hour, minute, second].
%    tarr[6]      Arrival date [year, month, day, hour, minute, second].
%
%  FUNCTIONS CALLED:
%    astroConstants.m
%    date2mjd2000.m
%    generaltransfer.m
%    constr.m
%    mjd20002date.m
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Menicaglia
%
% VERSIONS:
%   2021-01: Last version

%% Preliminary paramenters definition           
Rfb=astroConstants(20+id_fb);         % Mean radius of fly-by planet [km]
h_atm=[0,350,12,10.8,5000,59.5,19.7]; % Atmospheric height of the planets of Solar system [km]
rpmin=Rfb+h_atm(id_fb);               % Minimum feaseble rp at fly-by [km]

%% Time vectors definitions
tdep_1_mjd2000=date2mjd2000(tw_dep(1,:)); % Conversion from the date in the Gregorian calendar to the modified Julian day 2000 number of the given date
tdep_2_mjd2000=date2mjd2000(tw_dep(2,:));

tfb_1_mjd2000=date2mjd2000(tw_fb(1,:));
tfb_2_mjd2000=date2mjd2000(tw_fb(2,:));

tarr_1_mjd2000=date2mjd2000(tw_arr(1,:));
tarr_2_mjd2000=date2mjd2000(tw_arr(2,:));

%% Minimization
lb = [tdep_1_mjd2000;tfb_1_mjd2000;tarr_1_mjd2000]; 
ub = [tdep_2_mjd2000;tfb_2_mjd2000;tarr_2_mjd2000];

A = [];
b = [];
Aeq = [];
beq= [];

fun = @(x) generaltransfer(id_dep,id_fb,id_arr,x(1),x(2),x(3));
constr_nlin = @(x) constr(id_dep,id_fb,id_arr,x(1),x(2),x(3),rpmin);

options = optimoptions('fmincon','OptimalityTolerance',1e-14,'StepTolerance',1e-14);
x_opt = fmincon(@(x) fun(x),x0,A,b,Aeq,beq,lb,ub,constr_nlin,options);

%% Optimum solution
Dv_tot = generaltransfer(id_dep,id_fb,id_arr,x_opt(1),x_opt(2),x_opt(3));

tdep=mjd20002date(x_opt(1));
tfb=mjd20002date(x_opt(2));
tarr=mjd20002date(x_opt(3));

Dt=x_opt(3)-x_opt(1);

end

