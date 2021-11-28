function [ DV,DV1,DV2,vti,vtf] = singletransfer(r1_vec, r2_vec,v1_vec,v2_vec,t1,t2,muSun) 
% singletransfer.m - Computes the DV required for a trasfer between two
%                    planets                          
%
% PROTOTYPE:
%   [ DV,DV1,DV2,vti,vtf] = singletransfer(r1_vec, r2_vec,v1_vec,v2_vec,t1,t2,muSun) 
%
% DESCRIPTION:
%   This function Computes the overall DV, for a transfer between two
%   different planets using the patched conics method 
%    
%  INPUT:
%    r1_vec[3x1]          position vector of the departure plantet at time t1 [km]                
%    r2_vec[3x1]          position vector of the arrival planet at time t2 [km]              
%    v1_vec[3x1]          velocity vector of the departure planet at time t1 [km/s]
%    v2_vec[3x1]          velocity vector of the arrival planet at time t2 [km/s]
%    t1[1]              initial time of the transfer orbit [s] 
%    t2[1]              final time of the transfer orbit [s]
%    muSun[1]           Sun gravitational parameter [km^3/s^2]

%  OUTPUT: 
%    DV[1]           Overall delta V  of the interplanetary transfer [km/s]
%    DV1[1]          Delta V required at the departure in order to insert the
%                    spacecraft into an interplanetary orbit  [km/s]
%    DV2[1]          Delta V required at the arrival in order to insert the
%                    spacecraft into an interplanetary orbit [km/s]
%    vti[3x1]        Velocity at the departure of transfer orbit [km/s]
%    vtf[3x1]        Velocity at the arrival of the transfer orbit [km/s]
%
%  FUNCTIONS CALLED:
%    lambertMR.m 
%
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Menicaglia
%
% VERSIONS:
%   2021-01: Last version

Day2Seconds =24*60*60;
Dt = (t2-t1)*Day2Seconds;
%% Lambert function
[~,~,~,~,VI,VF,~,~] = lambertMR(r1_vec,r2_vec,Dt,muSun);
vti = VI;
vtf = VF;
DV1 = norm(VI-v1_vec');
DV2 = norm(v2_vec'-VF );
DV = DV1 + DV2;
end
