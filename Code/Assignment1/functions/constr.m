function [c,ceq] = constr(departureID,flybyID,arrivalID, t_dep_vec_MJD2000, t_flyby_vec_MJD2000,t_arr_vec_MJD2000,rpmin)

% PROTOTYPE:
%   [c,ceq] = constr(departureID,flybyID,arrivalID, t_dep_vec_MJD2000, t_flyby_vec_MJD2000,t_arr_vec_MJD2000,rpmin)
%
% DESCRIPTION:
% This function generate constrains needed by the fmincon MATLAB function 
%
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
%    c[1]                Difference between the minimum radius of perigee 
%                        possible and the radius of perigee of the flyb [km]
%
%  FUNCTIONS CALLED:
%    generaltransfer.m 
%
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-01: Last version

[~,~,~,flyby] = generaltransfer(departureID,flybyID,arrivalID,t_dep_vec_MJD2000,t_flyby_vec_MJD2000,t_arr_vec_MJD2000);
c(1) = rpmin - flyby.rp;
ceq = [];
end
