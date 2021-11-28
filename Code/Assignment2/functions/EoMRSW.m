function dkep = EoMRSW(t,kep,kP,acc_per)

% EoMRSW.m - system of ordinary differential equations bonding the
%            keplerian elements, their time evolution, the properties of
%            the primary and the included perturbations to the time
%            derivative of the keplerian elements.
%
% PROTOTYPE:
%   dkep = EoMRSW(t,kep,kP,acc_per)
%
% INPUT:
%   t         (1x1)   Time instant at which the state is evaluated    [s] 
%   kep       (6x1)   Keplerian parameters as a variable for the 
%                     function [km,-,rad,rad,rad,rad]
%   kP                Planetary gravitational constant             [km^3/s^2]
%   acc_per           Function handle of the total perturbing acceleration
%
% OUTPUT:
%   dkep    (6x1)   Time derivative of the keplerian elements
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version

% Extraction of Keplerian Elements from the state vector
a = kep(1);
e = kep(2);
i = kep(3);
om = kep(5);
th = kep(6);

% Computation of useful parameters
p = a*(1-e^2);
r = p/(1+e*cos(th));
h = sqrt(p*kP);

% Extraction of perturbing acceleration from the function used as input
a_per = acc_per(t,kep);
ar = a_per(1);
as = a_per(2);
aw = a_per(3);

% Gauss formulae for the variation of the keplerian elements
da = 2*a^2*(e*sin(th)*ar + p*as/r)/h;

de = (p*sin(th)*ar + ((p+r)*cos(th) + r*e)*as)/h;

di = r*cos(th+om)*aw/h;

dOm = r*sin(th+om)*aw/(h*sin(i));

dom = (-p*cos(th)*ar+(p+r)*sin(th)*as)/(h*e)-r*sin(th+om)*cos(i)*aw/(h*sin(i));

dth = h/r^2+(p*cos(th)*ar-(p+r)*sin(th)*as)/(e*h);

dkep = [da; de; di; dOm; dom; dth];
end