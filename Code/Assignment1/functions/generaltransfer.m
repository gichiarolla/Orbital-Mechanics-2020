function [DV,DV_dep,DV_arr,flyby,transfer1,transfer2] = generaltransfer(departureID,flybyID,arrivalID,t_dep_vec_MJD2000,t_flyby_vec_MJD2000,t_arr_vec_MJD2000)

% generaltransfer.m - Computes the DV required for a trasfer from
%                     departureID planet to arrivalID planet with a flyby
%                     at the flybyID planet
%
% PROTOTYPE:
%   [DV,DV_dep,DV_arr,flyby,transfer1,transfer2] = generaltransfer(departureID,flybyID,arrivalID,t_dep_vec_MJD2000,t_flyby_vec_MJD2000,t_arr_vec_MJD2000) 
%
% DESCRIPTION:
%   This function Computes the overall DV, the departure DV the arrival DV
%   the powered flyby properties the transfer from departure planter to fly
%   by planet properties and the transfer from flyby plantet to arrival    
%   planet properties.
%   The inputs of this function are the deparutre planet, the arrival
%   planet, the flyby planet and the vector containing the  initial and final 
%   time of the transfer for both the first transfer arc and the second 
%   transfer arc afther the flyby.

%  INPUT:
%    deparutreID[1]           ID of the departure planet
%    flybyID[1]               ID of the flyby planet
%    arrivalID[1]             ID of the arrival planet
%    t_dep_vec_MJD2000[1xn]   Vector containing the departure time window in
%                             MJD2000 [s] (n = number of steps of the window)                  
%    t_flyby_vec_MJD2000[1xn] Vector containing the flyby time window in MJD2000 [s] (n = number of steps of the window)  
%    t_arr_vec_MJD2000[1xn]   Vector containing the arrival time window in
%                             MJD2000 [s] (n = number of steps of the window)  
%    rpmin[1]                 Minimum radius of perigee of the flyby [km]
%
%  OUTPUT:
%    DV[1]         Overall delta V  of the interplanetary trransfer [km/s]
%    DV_dep[1]     Delta V required at the departure in order to insert the
%                  spacecraft into an interplanetary orbit [km/s]
%    DV_arr[1]     Delta V required at the arrival in order to insert the
%                  spacecraft into an interplanetary orbit [km/s]
%    flyby         Handle with all flyby properties
%
%  FUNCTIONS CALLED:
%    singletransfer.m 
%    uplanet.m
%    flyby_GA.m 
%    kep2car.m
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Menicaglia
%
% VERSIONS:
%   2021-01: Last version

muSun = astroConstants(4);
muflyby = astroConstants(12);

for k1 = 1:length(t_dep_vec_MJD2000)
    kep = uplanet( t_dep_vec_MJD2000(k1), departureID);
    [r1,v1] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun);
    s = [r1',v1']';
    r_dep_vec(1:3,k1) = s(1:3);
    v_dep_vec(1:3,k1) = s(4:6);
end
% Now I get flyby ephemerides
% here I discretize the time window
for k1 = 1:length(t_flyby_vec_MJD2000)
    kep = uplanet( t_flyby_vec_MJD2000(k1), flybyID);
    [r3,v3] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6), muSun);
    s = [r3',v3']';
    r_flyby_vec(1:3,k1) = s(1:3);
    v_flyby_vec(1:3,k1) = s(4:6);
end
% Now I get arrival ephemerides

for k1 = 1:length(t_arr_vec_MJD2000)
    kep = uplanet( t_arr_vec_MJD2000(k1), arrivalID);
    [r2,v2] = kep2car(kep(1),kep(2),kep(3),kep(4),kep(5),kep(6),muSun);
    s = [r2',v2']';
    r_arr_vec(1:3,k1) = s(1:3);
    v_arr_vec(1:3,k1) = s(4:6);
end
%% Calculate first leg 
[ ~,DV1_1,~,transfer1.vdep,v_flyby_minus] = ...
    singletransfer(r_dep_vec, r_flyby_vec,v_dep_vec,v_flyby_vec,t_dep_vec_MJD2000,t_flyby_vec_MJD2000,muSun);
transfer1.vflyby = v_flyby_minus;
transfer1.rdep = r_dep_vec;
transfer1.rflyby = r_flyby_vec;
transfer1.tdep = t_dep_vec_MJD2000;
transfer1.tflyby = t_flyby_vec_MJD2000;


%% Calculate second leg 
[~,~,DV2_2,v_flyby_plus,transfer2.varr] = singletransfer(r_flyby_vec, r_arr_vec,v_flyby_vec,v_arr_vec,t_flyby_vec_MJD2000,t_arr_vec_MJD2000,muSun);
transfer2.vflyby = v_flyby_plus;
transfer2.rflyby = r_flyby_vec;
transfer2.rarr = r_arr_vec;
transfer2.tflyby = t_flyby_vec_MJD2000;
transfer2.tarr = t_arr_vec_MJD2000;

%% Caluclate powered flyby
[Delta,rp,~,Dv_p,~,em,ep,am,ap,ym,~] = Flyby_GA(v_flyby_minus'-v_flyby_vec,v_flyby_plus'-v_flyby_vec,2);
% [Dv_p, delta,a,e,Delta,rp,vp0] = flybypow(v_flyby_minus'-v_flyby_vec,v_flyby_plus'-v_flyby_vec,muflyby);
flyby.Dv_p = Dv_p;
flyby.a = [am,ap];
flyby.e = [em,ep];
flyby.Delta = Delta;
flyby.rp = rp;
flyby.vp0 = ym(3:6);

DV = DV1_1 +DV2_2+ abs(Dv_p);
DV_dep = 0;
DV_arr = 0;


