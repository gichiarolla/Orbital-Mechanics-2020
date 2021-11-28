close all
clear 
clc
format long

path(pathdef);
addpath(genpath('functions\'));

% common parameteres and physical constants

muE = astroConstants(13);   %Earth's grav. param. [km^3/s^2]
RE = astroConstants(23);    %Earth's radius [km]
J2 = astroConstants(9);     

a = 11111;
e = 0.3939;
i = deg2rad(44.0523);
OM = deg2rad(60);
om = deg2rad(30);
f0 = deg2rad(60);
kep0 = [a, e, i, OM, om, f0];

t0=0;
om_Edeg = 15.04;                   %Earth's rotation vleocity [deg/h]
om_E = deg2rad(om_Edeg)/3600;      %Earth's rotation vleocity [rad/s]
T = 2*pi*sqrt((a^3)/muE);

k=15;
m=2;

CD = 2.1;                       % drag coeff
AtoMratio_km = 0.150*1e-6;      % A_cross: reference area [km^2/kg]

%orbit 3D 
plot3D_initial_orbit (a,e,i,OM,om,f0, muE, T, RE);

I = 1;          % param in order to loop
while I
    clc
    prompt_general = sprintf('\nOBTAIN DATA: \n"1" -> Ground Track\n"2" -> Orbit Propagation \n"3" -> Real Data Comparison - Arase (ERG) \n"4" -> CPU Time. \tN.B.: more than 30'' to run\nAny other input -> End the script \n\n Enter a data set: ');
    switchgeneral = input(prompt_general);

    switch(switchgeneral)

        case 1                  % general switch
        % Part II
        % a) Ground Track
        prompt_GT = sprintf('\nGROUND TRACK PLOTS: \n"1" -> One figure each GT\n"2" -> One figure each TYPE of GT and each TIME \n\n Enter GT number: ');
        switchGT = input(prompt_GT);

        switch switchGT 

            case 1
                % not perturbed
                T_1Day_10Days_Unp(a,e,i,OM,om,f0,muE,om_E,t0,T,RE,'one');

                % secular
                T_1Day_10Days_Sec(a,e,i,OM,om,f0,muE,om_E,t0,T,m,k,J2,RE,'one');


                % b) Repeating Ground Track
                % not perturbed
                T_1Day_10Days_Rep(a,e,i,OM,om,f0,muE,om_E,t0,T,m,k,RE,'one');

                % secular
                T_1Day_10Days_RepSec(a,e,i,OM,om,f0,muE,om_E,t0,T,m,k,J2,RE,'one');


            % One figure for each type of GT and for each time of orbit 
            case 2

                GT_comparison(a,e,i,OM,om,f0,muE,om_E,t0,T,m,k,J2,RE);

        end

        % Part IV
        case 2                  % general switch

            prompt_plot = sprintf('\nORBIT PROPAGATION\n"1" -> Plot Orbit evolution\n"2" -> Plot Keplerian Elements \n"3" -> Filtered Keplerian Elements\n"4" -> Error Gauss-Cartesian Equations\n"5" -> Orbit Evolution (Only Drag Effects). \tN.B.: more than 30'' to run\n"6" -> Movie\nAny other input -> Propagation OFF\n\n Enter an option: ');
            switchplot = input(prompt_plot);

            % set orbits' plots  
            switch(switchplot)
                case 5
                    ntotOrb = 1e5;
                    sizet = 1e5;
                    J2 = 0;
                otherwise
                    ntotOrb = 200;
                    sizet = 5e4;
            end

            % a) Propagate orbit GAUSS
            tspan = linspace(0,ntotOrb*T,sizet);
            [T_Gauss,Kep_Gauss] = propagatorgauss(tspan,kep0,muE,RE,J2,CD,AtoMratio_km,om_E);

            %compute r & v
            nvals = numel (T_Gauss);
            R_Gauss = zeros (nvals,3);
            V_Gauss = zeros (nvals,3);

            for j=1:nvals
                aa=Kep_Gauss(j,1);
                ee=Kep_Gauss(j,2);
                ii=Kep_Gauss(j,3);
                OO=Kep_Gauss(j,4);
                oo=Kep_Gauss(j,5);
                ff=Kep_Gauss(j,6);
                [R_Gauss(j,1:3),V_Gauss(j,1:3)] = kep2car (aa,ee,ii,OO,oo,ff, muE);
            end

            % b) Propagate orbit CARTESIAN
            switch(switchplot)              % compute Cart just when needed
                case 1
                case 3
                case 5
                otherwise
                    [r0,v0] = kep2car(a,e,i,OM,om,f0,muE);
                    condition0 = [r0;v0];
                    [T_Cart,Y_Cart] = propagatorcar(tspan,condition0,muE,RE,J2,CD,AtoMratio_km,om_E);

                    %compute Cartesian Elements
                    nvals = numel (T_Cart);
                    R_Cart = zeros (nvals,3);
                    V_Cart = zeros (nvals,3);
                    Kep_Cart = zeros (nvals,6);

                    for j=1:nvals
                            R_Cart(j,1:3) = Y_Cart(j,1:3);
                            V_Cart(j,1:3) = Y_Cart(j,4:6);
                            [aa,ee,ii,OO,oo,ff] = car2kep(R_Cart(j,1:3),V_Cart(j,1:3),muE);
                            Kep_Cart(j,:) = [aa,ee,ii,OO,oo,ff];
                    end
            end

            % Part V+VII 

            switch switchplot

                case 1 
                    nOrb = 10;
                    plot_orbits_evolution (RE,Kep_Gauss,ntotOrb,nOrb,tspan,muE,'$\Delta t$ [days]')

                case 2
                    % 6.a) 2 methods' elements comparison
                    plot_Kep_el (T_Gauss,Kep_Gauss,T_Cart,Kep_Cart,T);

                case 3
                    % 7.a) Filtering on Gauss equations
                    Tfilter = 4*T;
                    nwindow = nearest(Tfilter/(sum(diff(T_Gauss))/(numel(T_Gauss)-1)));
                    kepGfiltered = movmean(Kep_Gauss,nwindow,1);

                    % 6.a+7.b) plots for the Keplerian elements + Filtered
                    plot_Kep_el (T_Gauss,Kep_Gauss,T_Gauss,kepGfiltered,T,1);

                case 4
                    % 7) Comparisons and errors definition
                    [errR_a,errR_e,errR_i,errR_RAAN,errR_omega,errR_TA] = plot_comparison_Kep_el (a,T,Kep_Gauss,Kep_Cart,tspan);

                    fprintf('ERRORS in Norm Infinite: \na = %7.6e \ne = %7.6e \ni = %7.6e \nRAAN = %7.6e \nomega = %7.6e \nTA = %7.6e \n',errR_a,errR_e,errR_i,errR_RAAN,errR_omega,errR_TA);
                    pause(10);
                    
                case 5
                    nOrb = 10000;
                    plot_orbits_evolution (RE,Kep_Gauss,ntotOrb,nOrb,tspan,muE,'$\Delta t$ [years]');

                case 6
                    M = sizet/ntotOrb+1;   %1 orbit's points
                    orbitsMovie (R_Cart(1:M,:),R_Cart(end+1-M:end,:),R_Gauss(1:M,:),R_Gauss(end+1-M:end,:),tspan(1:M),RE,T,om_Edeg/3600,M);
            end

        % Part VIII

        case 3                  % general switch
            % a) stellite orbital elements
            [a_vect,e_vect,i_vect,RAAN_vect,omega_vect,M_vect,dateSec,dateMJD] = realSC_Elem;

            nvals = numel (dateSec);
            Kep_SC = zeros(nvals,6);
            for j=1:nvals
                Kep_SC(j,1) = a_vect(j);
                Kep_SC(j,2) = e_vect(j);
                Kep_SC(j,3) = i_vect(j);
                Kep_SC(j,4) = RAAN_vect(j);
                Kep_SC(j,5) = omega_vect(j);
                Kep_SC(j,6) = M_vect(j);
            end

            % b) propagate using the model
            a_SC = a_vect(1);
            e_SC = e_vect(1);
            irad_SC = i_vect(1);
            RAAN_SC = RAAN_vect(1);
            omega_SC = omega_vect(1);
            M_SC = M_vect(1);
            T_SC = 2*pi*sqrt(a_SC^3/muE);
            n_SC = 2*pi/T_SC;
            t0_SC = M_SC/n_SC;
            theta0_SC = TrueAnomaly(t0_SC,e_SC,a_SC,muE,0,0);
            kep0_SC = [a_SC,e_SC,irad_SC,RAAN_SC,omega_SC,wrapTo2Pi(deg2rad(theta0_SC))];
            CD = 2.1;                                       % [-]
            AtoMratioShort = (1.5*1.5+5*1)/350;             % [m^2/kg]
            AtoMratioLong = (1.5*2.7+5*1)/350;
            AtoMratio = (AtoMratioShort+AtoMratioLong)/2;
            AtoMratio_km = AtoMratio*1e-6;                  % [km^2/kg]
            
            sizet = 1e5;
            tspan = linspace(dateSec(1),dateSec(end),sizet);
            [tgauss_SC,kepgauss_SC] = propagatorgauss(tspan,kep0_SC,muE,RE,J2,CD,AtoMratio_km,om_E);

            % c) compare results with downloaded elements
            plot_SC_el (dateMJD,kepgauss_SC,Kep_SC,sizet,T,tgauss_SC,'Gauss');

       case 4                   % general switch

            % Cpu time comparison between the two methods
            n = linspace(1000,100000,100);
            plotCPUtime(n,T,kep0,muE,RE,J2,CD,AtoMratio_km,om_E);

        otherwise               % general switch
            I = 0;

    end
end