function [a_vect,e_vect,i_vect,RAAN_vect,omega_vect,M_vect,dateSec,dateMJD] = realSC_Elem ()

% realSC_Elem.m - exctract ephemerides form txt file downloaded from SpaceTrack. 
%
% PROTOTYPE:
%   [a_vect,e_vect,i_vect,RAAN_vect,omega_vect,M_vect,dateSec,dateMJD] = realSC_Elem ()
%
% INPUT:
%   no inputs because it takes datas his own
%
% OUTPUT:
%   a_vect        [nx1]       Semi-major axis vector          [km]
%   e_vect        [nx1]       Eccentricity axis vector        [-]
%   i_vect        [nx1]       Inclination vector              [rad]
%   RAAN_vect     [nx1]       RAAN vector                     [rad]
%   omega_vect    [nx1]       Argument of periapsis vector    [rad]
%   M_vect        [nx1]       True anomaly vector             [rad]
%   dateSec       [nx1]       Time vector                     [s]
%   dateMJD       [1xn]       Time vector in MJD              [days]
%
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version


% Reading Ephemerides from the .txt file as obtained by 
% www.space-track.org
fileID = fopen('ephemerides41896.txt');
Ephcell = textscan(fileID, '%s %s');

% Definition of cells containing the useful datas as extracted from the
% Ephemerides and transposition of the cells
EpochCell = {Ephcell{2}{15:41:end}};
dateMJD = zeros(length(EpochCell),1);
a_cell = {Ephcell{2}{30:41:end}};
e_prev = {Ephcell{2}{17:41:end}};
i_cell = {Ephcell{2}{18:41:end}};
RAAN_cell = {Ephcell{2}{19:41:end}};
omega_cell = {Ephcell{2}{20:41:end}};
M_cell = {Ephcell{2}{21:41:end}};

a_cellT = cellfun(@transpose,a_cell,'UniformOutput',false);
i_cellT = cellfun(@transpose,i_cell,'UniformOutput',false);
RAAN_cellT = cellfun(@transpose,RAAN_cell,'UniformOutput',false);
omega_cellT = cellfun(@transpose,omega_cell,'UniformOutput',false);
M_cellT = cellfun(@transpose,M_cell,'UniformOutput',false);
Y0Cell = textscan(EpochCell{1}(2:5),'%f');
M0Cell = textscan(EpochCell{1}(7:8),'%f');
D0Cell = textscan(EpochCell{1}(10:11),'%f');
h0Cell = textscan(EpochCell{1}(13:14),'%f');
m0Cell = textscan(EpochCell{1}(16:17),'%f');
s0Cell = textscan(EpochCell{1}(19:end),'%f');
date0 = cell2mat([Y0Cell,M0Cell,D0Cell,h0Cell,m0Cell,s0Cell]);

% For cycle defining cells containg the time instants
dateSec = zeros(1,length(EpochCell));
for ktime = 1:length(EpochCell)
    % Definition of the time distance from the first date at which the
    % ephemerides are generated to the Epoch of the k-th registered set of
    % datas
    Epochk = EpochCell{ktime};
    Year = textscan(Epochk(2:5),'%f');
    Month = textscan(Epochk(7:8),'%f');
    Day = textscan(Epochk(10:11),'%f');
    Hour = textscan(Epochk(13:14),'%f');
    Minute = textscan(Epochk(16:17),'%f');
    Second = textscan(Epochk(19:end),'%f');
    dateCell = [Year,Month,Day,Hour,Minute,Second];
    date = cell2mat(dateCell);
    
    dateSec(ktime) = diffsecfromdate0(date,date0);

    dateMJD(ktime) = date2mjd2000(date);
end

% Definition of cells that can be converted in floating point numbers
for kdef = 1:length(dateSec)
    a_fpn(kdef) = textscan(a_cellT{kdef}(2:end),'%f');
    e_cell(kdef) = textscan(e_prev{kdef}(2:end),'%f');
    i_fpn(kdef) = textscan(i_cellT{kdef}(2:end),'%f');
    RAAN_fpn(kdef) = textscan(RAAN_cellT{kdef}(2:end),'%f');
    omega_fpn(kdef) = textscan(omega_cellT{kdef}(2:end),'%f');
    M_fpn(kdef) = textscan(M_cellT{kdef}(2:end),'%f');
end

% Definition of numerical vectors for the ephemerides
a_vectT = cell2mat(a_fpn);
e_vect = cell2mat(e_cell);
ideg_vect = cell2mat(i_fpn);
RAANdeg_vect = cell2mat(RAAN_fpn);
omegadeg_vect = cell2mat(omega_fpn);
Mdeg_vect = cell2mat(M_fpn);

a_vect = a_vectT';
i_vect = deg2rad(ideg_vect');
RAAN_vect = deg2rad(RAANdeg_vect');
omega_vect = deg2rad(omegadeg_vect');
M_vect = deg2rad(Mdeg_vect');

end

