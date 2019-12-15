function [ETA XMASSFLOW DATEN DATEX DAT MASSFLOW COMBUSTION FIG] = ST(P_e,options,display)
% ST Steam power plants modelisation
% ST(P_e,options,display) compute the thermodynamics states for a Steam
% power plant (combustion, exchanger, cycle) turbine based on several
% inputs (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated). Refer to Fig 2.33 from reference book (in english)
% P_E = electrical power output target [kW]
% OPTIONS is a structure containing :
%   -options.nsout     [-] : Number of feed-heating
%   -options.reheat    [-] : Number of reheating
%   -options.T_max     [°C] : Maximum steam temperature
%   -options.T_cond_out[°C] : Condenseur cold outlet temperature
%   -options.p3_hp     [bar] : Maximum pressure
%   -options.drumFlag  [-] : if =1 then drum if =0 => no drum.
%   -options.eta_mec   [-] : mecanic efficiency of shafts bearings
%   -options.comb is a structure containing combustion data :
%       -comb.Tmax     [°C] : maximum combustion temperature
%       -comb.lambda   [-] : air excess
%       -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
%   -options.T_exhaust [°C] : Temperature of exhaust gas out of the chimney
%   -options.p4       [bar] : High pressure after last reheating
%   -options.x6        [-] : Vapor ratio [gaseous/liquid] (in french : titre)
%   -options.T_0       [°C] : Reference temperature
%   -options.TpinchSub [°C] : Temperature pinch at the subcooler
%   -options.TpinchEx  [°C] : Temperature pinch at a heat exchanger
%   -options.TpinchCond[°C] : Temperature pinch at condenser
%   -options.Tdrum     [°C] : minimal drum temperature
%   -options.eta_SiC    [-] : Internal pump efficiency
%   -options.eta_SiT    [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
%             	             eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others
% DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_gen, Steam generator energy efficiency
%   -eta(6) : eta_gex, Steam generator exergy efficiency
%   -eta(7) : eta_combex, Combustion exergy efficiency
%   -eta(8) : eta_chemex, Chimney exergy efficiency (losses)
%   -eta(9) : eta_condex, Condenser exergy efficiency
%   -eta(10): eta_transex, Steam bleeding heat exchanger overall exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% Xmassflow is a vector with each feedheating massflow [kg/s] (with respect to figure
%           2.33, page 91 "Thermal Power Plants" English version).
%           Xmassflow(1) = mass flow at 6_1 etc...
% DATEN is a vector with :
%   -daten(1) : perte_gen [kW]
%   -daten(2) : perte_mec [kW]
%   -daten(3) : perte_cond [kW]
% DATEX is a vector with :
%   -datex(1) : perte_mec    [kW]
%   -datex(2) : perte_totex  [kW]
%   -datex(3) : perte_rotex  [kW]
%   -datex(4) : perte_combex [kW]
%   -datex(5) : perte_condex [kW]
%   -datex(6) : perte_chemex [kW]
%   -datex(7) : perte_transex[kW]
% DAT is a matrix containing :
% dat = {T_1       , T_2       , ...       , T_6_I,     T_6_II, ... ;  [°C]
%        p_1       , p_2       , ...       , p_6_I,     p_6_II, ... ;  [bar]
%        h_1       , h_2       , ...       , h_6_I,     h_6_II, ... ;  [kJ/kg]
%        s_1       , s_2       , ...       , s_6_I,     s_6_II, ... ;  [kJ/kg/K]
%        e_1       , e_2       , ...       , e_6_I,     e_6_II, ... ;  [kJ/kg]
%        x_1       , x_2       , ...       , x_6_I,     x_6_II, ... ;   };[-]
% MASSFLOW is a vector containing :
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_v, water massflow at 2 [kg/s]
%   -massflow(3) = m_c, combustible massflow [kg/s]
%   -massflow(4) = m_f, exhaust gas massflow [kg/s]
% 
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combustible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas     [kJ/kg/K]
%   -combustion.fum  : is a vector of the exhaust gas composition :
%       -fum(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s]
%       -fum(2) = m_N2f  : massflow of N2 in exhaust gas [kg/s]
%       -fum(3) = m_CO2f : massflow of CO2 in exhaust gas [kg/s]
%       -fum(4) = m_H2Of : massflow of H2O in exhaust gas [kg/s]
%
% FIG is a vector of all the figure you plot. Before each figure, define a
% figure environment such as:
%  "FIG(1) = figure;
%  plot(x,y1);
%  [...]
%   FIG(2) = figure;
%  plot(x,y2);
%  [...]"
%  Your vector FIG will contain all the figure plot during the run of this
%  code (whatever the size of FIG).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS VERIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin<3
        display = 1;
        if nargin<2
            constr = "pressure";
            %constr = "vapor_ratio";
            options = struct();
            if nargin<1
                P_e = 288; % [MW] Puissance énergétique de l'installation
            end
        end
    end

    % options.reheat [-] : Number of reheating
    if isfield(options,'reheat') %OK
        reheat = options.reheat;
    else
        reheat = 1;
    end

    if reheat < 0
        fprintf('\n    WARNING : there cannot be a negative amount of reheating.');
        fprintf('\n              Its value is reinitialized to 0.');
        reheat = 0;
    end

    % options.nsout [-] : Number of feed-heating
    if isfield(options,'nsout') %OK
        nsout = options.nsout;
    else
        nsout = 8;
    end

    if nsout < reheat && isfield(options,'reheat')
        fprintf('\n    WARNING : It must be at least the same number of tapping as there are reheatings.');
        fprintf('\n              Its value is reinitialized to ''reheat''.');
        nsout = reheat;
    elseif nsout < reheat
        reheat = nsout;
    end

    % options.T_max [°C] : Maximum steam temperature
    if isfield(options,'T_max')
        T_max = options.T_max;
        if T_max > 1e3
            fprintf('\n    WARNING : The maximum temperature for the water is 560°C because');
            fprintf('\n              the container begin to smell at 565°C.');
            T_max = 560;
        end
    else
        T_max = 565;
    end

    % options.T_cond_out [°C] : Condenseur cold outlet temperature
    if isfield(options,'T_cond_out') %OK
        T_cond_out = options.T_cond_out;
        if T_cond_out <= 15
            fprintf('\n    WARNING : L''eau de sortie de refroidissement ne peut quitter le condenseur plus froide que ce qu''elle n''est arrivee');
            fprintf('\n              On retablira donc sa valeur à 20 degre Celsius');
            T_cond_out = 20
        end
    else
        T_cond_out = 23;
    end

    % options.p3_hp [bar] : Maximum pressure
    if isfield(options,'p3_hp') %OK
        p3_hp = options.p3_hp;
    else
        p3_hp = 350;
    end

    % options.drumFlag [-] : if =1 then drum if =0 => no drum.
    if isfield(options,'drumFlag')
        drumFlag = options.drumFlag;
        if drumFlag > .5
            drumFlag = 1;
        else
            drumFlag = 0;
        end
    else
        drumFlag = 1;
    end

    %   -options.eta_mec [-] : mecanic efficiency of shafts bearings
    if isfield(options,'eta_mec')
        eta_mec = options.eta_mec;
    else
        eta_mec = .98;
    end

    % options.x6 [-] : Vapor ratio [gaseous/liquid] (in french : titre)
    if isfield(options,'x_6')
        x_6 = options.x_6;
        constr = "vapor_ratio";
    else
        x_6 = .89;
    end

    % options.p4 [bar] : High pressure after last reheating
    if isfield(options,'p_4')
        p_4 = options.p_4;
        constr = "pressure";
    else
        p_4 = 70;
    end

    if ~exist('constr')
        constr = "pressure";
    end

    % options.T_0 [°C] : Reference temperature
    if isfield(options,'T_0')
        T_0 = options.T_0;
    else
        T_0 = 15;
    end

    % options.TpinchSub [°C] : Temperature pinch at the subcooler
    if isfield(options,'TpinchSub')
        TpinchSub = options.TpinchSub;
    else
        TpinchSub = 4;
    end


    % options.TpinchEx [°C] : Temperature pinch at a heat exchanger
    if isfield(options,'TpinchEx') %OK
        TpinchEx = options.TpinchEx;
    else
        TpinchEx = 1;
    end

    % options.TpinchCond [°C] : Temperature pinch at condenser
    if isfield(options,'TpinchCond') %OK
        TpinchCond = options.TpinchCond;
    else
        TpinchCond = 10;
    end

    % options.Tdrum [°C] : minimal drum temperature
    if isfield(options,'Tdrum')
        Tdrum = options.Tdrum;
    else
        Tdrum = 148.7;
    end

    % options.eta_SiC [-] : Internal pump efficiency
    if isfield(options,'eta_SiC')
        eta_SiC = options.eta_SiC;
    else
        eta_SiC = .85;  % [-]
    end

    %   -options.eta_SiT [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
    %                          eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others
    if isfield(options,'eta_SiT')
        eta_SiT = options.eta_SiT;
    else
        eta_SiT = [.927,.9];
    end

    if min(eta_SiT) < .88
        fprintf('\n    WARNING : if the isentropic efficiency of the turbines is lower than 88%,');
        fprintf('\n              because of the cavitation that can damage the turbine.');
        fprintf('\n          >>> Its value is reinitialized to 88%.');
    end
    while min(eta_SiT) < .88
        eta_SiT(find(eta_Sit == min(eta_Sit))) = .88;
    end

    if length(eta_SiT) == 2
        eta_SiT1 = eta_SiT(1);
        eta_SiT2 = eta_SiT(2);
    else
        eta_SiT1 = eta_SiT;
        eta_SiT2 = eta_SiT;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OTHER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M_O2  = 31.99800e-3;  % [kg/mol]
    M_N2  = 28.01400e-3;  % [kg/mol]
    M_CO2 = 44.00800e-3;  % [kg/mol]
    M_H2O = 18.01494e-3;  % [kg/mol]
    M_air = .21*M_O2 + .79*M_N2;  % [kg/mol]

    p_ext = 1.01325;  % [bar]
    T_sat_max = 373.9458;
    p_sat_max = XSteam('psat_T',T_sat_max);

    if isfield(options,'lam_chaud')
        lam_chaud = options.lam_chaud;
    else
        lam_chaud = 1/1.1290;
    end

    if isfield(options,'lam_exch')
        lam_exch = options.lam_exch;
    else
        lam_exch  = 1;%.98;
    end

    T_vect = [ones(1,300)*300,301:5000];
    Cp_O2  = janaf('O2',T_vect); % [J/(kg*K)]
    Cp_N2  = janaf('N2',T_vect);
    Cp_CO2 = janaf('CO2',T_vect);
    Cp_H2O = janaf('H2O',T_vect);

    function [X] = Cp_air(T)
        T = T +273.15;
        T(T<300) = 300; T(T>5000) = 5000;
        X = .21*janaf('O2',T) + .79*janaf('N2',T);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REFERENCE STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p_0 = XSteam('psat_T',T_0);
    h_0 = XSteam('hV_T',T_0);
    s_0 = XSteam('sV_T',T_0);
    e_0 = 0;

    Ta_0 = 15;
    pa_0 = 1e5;
    ha_0 = 1.006 * Ta_0;
    sa_0 = 0.054;
    ea_0 = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMBUSTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Verification of parameters
        % options.comb is a structure containing combustion data : 
        %   -comb.Tmax     [°C] : maximum combustion temperature
        %   -comb.lambda   [-] : air excess
        %   -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
        %   -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
        if isfield(options,'comb')
            comb = options.comb;
            if isfield(comb,'Tmax')
                Tmax = comb.Tmax;
            else
                Tmax = 600;
            end
            if isfield(comb,'lambda')
                lambda = comb.lambda;
            elseif ~isfield(comb,'Tmax')
                lambda = 1.3;
            end
            if isfield(comb,{'x','y'})
                x = options.comb.x;
                y = options.comb.y;
            else
                x = 0; y = 4;
            end
        else
            Tmax = 600;
            lambda = 1.05;
            x = 0;
            y = 4;
        end

        if exist('lambda')
            constr_comb = "lambda";
        else
            constr_comb = "Tmax";
        end
        %constr_comb = "";
        %constr_comb = "Tmax";

        % options.T_exhaust [°C] : Temperature of exhaust gas out of the chimney
        if isfield(options,'T_exhaust')
            T_exhaust = options.T_exhaust;
        else
            T_exhaust = 120;
        end

    M_c  = (12.01+1.01*y+16*x)*1e-3; % Molar mass of carburant [kg/mol]
    LHV  = (393.6 + 102.2*y - (110.6 + 204.4*y)*x/(1+y/2)) / M_c; % [kJ/kg]
    T_CH4 = 298.15;
    coeff1 = [-.703029 108.4773 -42.52157 5.862788 .678565];
    Cp_CH4 = integral(@(T) coeff1(1)+coeff1(2)*(T./1000)...
        +coeff1(3)*((T./1000).^2)+coeff1(4)*((T./1000).^3)...
        +(coeff1(5)./((T./1000).^2)),273.15,T_CH4)/(T_CH4-273.15) *1e-3;

    % Find the air excess
    if constr_comb == "Tmax"
        a = @(lam) (lam-1)*(1+y/4-x/2); b = y/2;
        w = @(lam) lam*(1+y/4-x/2); % Stoechiometric coefficients

        t = 273:ceil(Tmax+273);
        fun = @(lam) Tmax * (mean(Cp_CO2(t))*M_CO2 + b*mean(Cp_H2O(t))*M_H2O ...
            + a(lam)*mean(Cp_O2(t))*M_O2 + 3.76*w(lam)*mean(Cp_N2(t))*M_N2) ...
            - w(lam)*Cp_air(300)*M_air/.21 - LHV*M_c - Cp_CH4*25;
        % Bilan d'enthalpie sur la combustion
        lambda = fsolve(fun,1)
    end

    a = (lambda-1)*(1+y/4-x/2); b = y/2;
    w = lambda*(1+y/4-x/2); % Stoechiometric coefficients

    % Find the max temperature in the combustion chamber
    if constr_comb == "lambda"
        Tmax = 1400; iter = 0;
        while abs(Tmax - iter) > 1
            iter = Tmax; t = 273:int16(Tmax+273.15);
            Tmax = (w*1.006/15*M_air/.21 + LHV*M_c + Cp_CH4*25) /...
                (mean(Cp_CO2(t))*M_CO2 + b*mean(Cp_H2O(t))*M_H2O ...
                + a*mean(Cp_O2(t))*M_O2 + 3.76*w*mean(Cp_N2(t))*M_N2);
        end
    end

    comp_f_tot = M_CO2 + b*M_H2O + a*M_O2 + 3.76*w*M_N2; % [kg]
    comp_f_CO2 = M_CO2 / comp_f_tot; % [-]
    comp_f_H2O = b*M_H2O / comp_f_tot; % [-]
    comp_f_O2  = a*M_O2 / comp_f_tot; % [-]
    comp_f_N2  = 3.76*w*M_N2 / comp_f_tot; % [-]

    function T = Cp_f(T1,T2)
        if T1 < 300 
            T1 = 300;
        end
        if T2 < 300 
            T2 = 300;
        end
        if T2 < T1
            error('Error T2 < T1');
        elseif T1 == T2
            T = (comp_f_O2*janaf('O2',T1) + comp_f_N2*janaf('N2',T1)...
                + comp_f_CO2*janaf('CO2',T1) + comp_f_H2O*janaf('H2O',T1));
        else
            decT1 = T1 - floor(T1); decT2 = T2 - floor(T2);
            if decT1 > decT2
                t = floor(T2) - floor(T1);
            elseif decT1 < decT2
                t = ceil(T2) - floor(T1);
            else 
                t = T2 - T1;
            end 

            R = ones(t+1,1); S = ones(t,1);
            R(1) = (comp_f_O2*janaf('O2',T1) + comp_f_N2*janaf('N2',T1)...
                + comp_f_CO2*janaf('CO2',T1) + comp_f_H2O*janaf('H2O',T1));
            for i = 1:t
                if (T1 + i) <= T2
                    R(i +1) = (comp_f_O2*janaf('O2',T1 + i) + comp_f_N2*janaf('N2',T1 + i)...
                        + comp_f_CO2*janaf('CO2',T1 + i) + comp_f_H2O*janaf('H2O',T1 + i));
                    S(i) = ((R(i)+R(i +1))/2)*1;
                else
                    j = abs(T2-(T1+(i-1)));
                    R(i +1) = (comp_f_O2*janaf('O2',T2) + comp_f_N2*janaf('N2',T2)...
                        + comp_f_CO2*janaf('CO2',T2) + comp_f_H2O*janaf('H2O',T2));
                    S(i) = ((R(i)+R(i +1))/2)*j;
                end 
            end
            T = sum(S)/(T2 - T1);
        end
    end

    m_a1 = (1+y/4-x/2) * (M_air/.21)/M_c; % [mol]
    m_ac = lambda * m_a1;  % [-]
    m_ag = (1 + 1/m_ac)^(-1); % [-]
    Q_comb = LHV / m_ac;

    % First state of comburant
    f = 1.041; e_c = f * LHV;
    Tmax
    % State after combustion
    T_h = Tmax; % Temperature after combustion [K]
    h_h = Cp_f(273.15,T_h+273.15) * T_h; % Enthalpy after combustion [kJ/kg]
    s_h = sa_0 + Cp_f(Ta_0+273.15,T_h+273.15)*log((Ta_0+273.15)/(T_h+273.15));
    e_h = (h_h-ha_0) - Ta_0*(s_h-sa_0);

    % Exhaust gasses
    T_e = T_exhaust; % Temperature after the heat exchange [K]
    h_e = Cp_f(273.15,T_e+273.15) * T_e; % Enthalpy after the heat exchange [kJ/kg]
    s_e = sa_0 + Cp_f(Ta_0+273.15,T_e+273.15)*log((Ta_0+273.15)/(T_e+273.15));
    e_e = (h_e-ha_0) - Ta_0*(s_e-sa_0);

    eps_ech = (m_ac/m_ag*h_e - m_ac*ha_0)/LHV;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONDENSOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % State in
    T_i = T_0; p_i = 10;
    h_i = XSteam('h_pT',p_i,T_i);
    s_i = XSteam('s_pT',p_i,T_i);
    e_i = (h_i-h_0) + (T_0+273.15)*(s_i-s_0);

    % State out
    T_o = T_cond_out; p_o = 10;
    h_o = XSteam('h_pT',p_o,T_o);
    s_o = XSteam('s_pT',p_o,T_o);
    e_o = (h_o-h_0) + (T_0+273.15)*(s_o-s_0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CYCLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STATE 3 - GIVEN BY PRESSURE AND TEMPERATURE
    T_3 = T_max; % given
    p_3 = p3_hp *lam_chaud; % given
    h_3 = XSteam('h_pT',p_3,T_3);
    s_3 = XSteam('s_pT',p_3,T_3);
    e_3 = (h_3-h_0) - (T_0+273.15)*(s_3-s_0);
    x_3 = NaN;
    state_3 = [T_3;p_3;h_3;s_3;e_3;x_3];

% STATE 7 - GIVEN BY TPinch AT THE CONDENSOR AND THE VAPOR RATIO (=.0)
    T_7 = T_cond_out + TpinchCond;
    p_7 = XSteam('psat_T',T_7);
    h_7 = XSteam('hL_T',T_7);
    s_7 = XSteam('sL_T',T_7);
    e_7 = (h_7-h_0) - (T_0+273.15)*(s_7-s_0);
    x_7 = .0;
    state_7 = [T_7;p_7;h_7;s_7;e_7;x_7];

% STATES 4,5,6 - WE KNOW THE TEMPERATURE AT 5 / WE KNOW OR VR AT 6 OR PRESSION AT 4
% WARNING : There would be an overstress if x_6 & p_4 are known
% Transfo 3 - 4 : isentropic relaxation with isentropic efficiency
% Transfo 4 - 5 : isobaric reheating
% Transfo 5 - 6 : isentropic relaxation with isentropic efficiency
    if reheat == 0
        T_6  = T_7;
        p_6  = XSteam('psat_T',T_6);
        h_6s = XSteam('h_ps',p_6,s_3);
        h_6  = h_3 - eta_SiT2 * (h_3 - h_6s); % eq2.9
        s_6  = XSteam('s_ph',p_6,h_6);
        e_6  = (h_6-h_0) - (T_0+273.15)*(s_6-s_0);
        x_6  = XSteam('x_ph',p_6,h_6);

        p_4 = []; p_5 = [];
        T_4 = []; T_5 = [];
        h_4 = []; h_5 = [];
        s_4 = []; s_5 = [];
        e_4 = []; e_5 = [];
        x_4 = []; x_5 = [];
        X_4 = []; X_5 = [];
        st_4 = {}; st_5 = {};
    elseif constr == "pressure"
        r = (p_4/p_3).^(1/reheat) * lam_chaud^(1-1/reheat);

        p_4  = p_3*r;
        h_4s = XSteam('h_ps',p_4,s_3);
        h_4  = h_3 - eta_SiT1 * (h_3 - h_4s); % eq2.9
        T_4  = XSteam('T_ph',p_4,h_4);
        s_4  = XSteam('s_ph',p_4,h_4);
        e_4  = (h_4-h_0) - (T_0+273.15)*(s_4-s_0);
        x_4  = NaN; st_4 = {'4'};

        p_5 = p_4 *lam_chaud; % isobare
        T_5 = T_max; % réchauffeur
        h_5 = XSteam('h_pT',p_5,T_5);
        s_5 = XSteam('s_pT',p_5,T_5);
        e_5 = (h_5-h_0) - (T_0+273.15)*(s_5-s_0);
        x_5 = NaN; st_5 = {'5'};

        n_rh = 1;
        while n_rh < reheat
            n_rh = n_rh +1;

            p_4  = [p_4,p_5(end)*r];
            h_4s = XSteam('h_ps',p_4(end),s_5(end));
            h_4  = [h_4,h_5(end)-eta_SiT1*(h_5(end)-h_4s)];
            T_4  = [T_4,XSteam('T_ph',p_4(end),h_4(end))];
            s_4  = [s_4,XSteam('s_ph',p_4(end),h_4(end))];
            e_4  = [e_4,(h_4(end)-h_0)-(T_0+273.15)*(s_4(end)-s_0)];

            p_5 = [p_5,p_4(end) *lam_chaud]; % isobare
            T_5 = [T_5,T_max]; % réchauffeur
            h_5 = [h_5,XSteam('h_pT',p_5(end),T_5(end))];
            s_5 = [s_5,XSteam('s_pT',p_5(end),T_5(end))];
            e_5 = [e_5,(h_5(end)-h_0) - (T_0+273.15)*(s_5(end)-s_0)];

            x_4  = [x_4,NaN]; x_5 = [x_5,NaN];
            if n_rh == 2
                st_4 = {'4_I','4_II'};
                st_5 = {'5_I','5_II'};
            else
                c = convertCharsToStrings(num2roman(n_rh));
                st_4 = cat(2,st_4,cellstr(strcat("4_",c)));
                st_5 = cat(2,st_5,cellstr(strcat("5_",c)));
            end
        end

        T_6  = T_7;
        p_6  = XSteam('psat_T',T_6);
        h_6s = XSteam('h_ps',p_6,s_5(end));
        h_6  = h_5(end) - eta_SiT2 * (h_5(end) - h_6s); % eq2.9
        s_6  = XSteam('s_ph',p_6,h_6);
        e_6  = (h_6-h_0) - (T_0+273.15)*(s_6-s_0);
        x_6  = XSteam('x_ph',p_6,h_6);
        v_6  = XSteam('v_ph',p_6,h_6);
    elseif constr == "vapor_ratio"
        T_6 = T_7;
        p_6 = XSteam('psat_T',T_6);
        h_6 = XSteam('h_Tx',T_6,x_6);
        s_6 = XSteam('sL_T',T_6)*(1-x_6) + XSteam('sV_T',T_6)*x_6;
        e_6 = (h_6-h_0) - (T_0+273.15)*(s_6-s_0);

        T_5e = T_max;
        h_5s = adjust(T_5e,0,s_6);
        h_5e = h_6 + eta_SiT2 * (h_5s-h_6);
        s_5e = adjust(T_5e,h_5e,0);
        p_5e = XSteam('p_hs',h_5e,s_5e);
        e_5e = (h_5e-h_0) - (T_0+273.15)*(s_5e-s_0);
        x_5e = NaN;

        r = (p_5e*lam_chaud/p_3).^(1/reheat) * lam_chaud^(1-1/reheat);

        p_4 = p_3 * r
        h_4s = XSteam('h_ps',p_4,s_3);
        h_4  = h_3 - eta_SiT1 * (h_3-h_4s); % eq2.9
        T_4  = XSteam('T_ph',p_4,h_4);
        s_4  = XSteam('s_ph',p_4,h_4);
        e_4  = (h_4-h_0) - (T_0+273.15)*(s_4-s_0);
        x_4  = NaN; st_4 = {'4'};

        T_5 = []; p_5 = []; h_5 = []; s_5 = []; e_5 = []; x_5 = []; st_5 = {};
        n_rh = 1;
        while n_rh < reheat
            n_rh = n_rh +1;

            p_5 = [p_5,p_4(end)*lam_chaud]; % isobare
            T_5 = [T_5,T_max]; % réchauffeur
            h_5 = [h_5,XSteam('h_pT',p_5(end),T_5(end))];
            s_5 = [s_5,XSteam('s_pT',p_5,T_5)];
            e_5 = [e_5,(h_5(end)-h_0) - (T_0+273.15)*(s_5(end)-s_0)];

            p_4  = [p_4,p_5(end)*r];
            h_4s = XSteam('h_ps',p_4(end),s_5(end));
            h_4  = [h_4,h_5(end)-eta_SiT1*(h_5(end)-h_4s)];
            T_4  = [T_4,XSteam('T_ph',p_4(end),h_4(end))];
            s_4  = [s_4,XSteam('s_ph',p_4(end),h_4(end))];
            e_4  = [e_4,(h_4(end)-h_0)-(T_0+273.15)*(s_4(end)-s_0)];

            x_4  = [x_4,NaN]; x_5 = [x_5,NaN];
            if n_rh == 2
                st_4 = {'4_I','4_II'};
                st_5 = {'5_I'};
            else
                c = convertCharsToStrings(num2roman(n_rh));
                st_4 = cat(2,st_4,cellstr(strcat("4_",c)));
                st_5 = cat(2,st_5,cellstr(strcat("5_",c)));
            end
        end

        T_5 = [T_5,T_5e]; p_5 = [p_5,p_5e];
        h_5 = [h_5,h_5e]; s_5 = [s_5,s_5e];
        e_5 = [e_5,e_5e]; x_5 = [x_5,x_5e];
        if isempty(st_5);
            st_5 = {'5'};
        else
            st_5 = cat(2,s_5,cellstr(strcat("5",...
                convertCharsToStrings(num2roman(n_rh)))));
        end
    else
        fprintf('\n  The parameter you choose is not resolvable.\n  Please set no more than 2 reheating, and define :\n');
        fprintf('    - if reheat == 0 : //\n    - if reheat == 1 : whether p_4 or x_6\n    - if reheat == 2 : p_4 and x_6.\n\n');
        return
    end
    state_4 = [T_4;p_4;h_4;s_4;e_4;x_4];
    state_5 = [T_5;p_5;h_5;s_5;e_5;x_5];
    state_6 = [T_6;p_6;h_6;s_6;e_6;x_6];

if nsout >= 1
% STATES 6_i - ISOENTHALPIC TAPPING BETWEEN STATES 5 AND 6
% WARNING : Before each reheating, we perform a tapping to decrease the quantity of
%           vapor we have to reheat.
    if reheat == 0
        h_6i = linspace(h_3,h_6,nsout-reheat+2)'; h_6i = h_6i(2:end-1);
        T_6i  = zeros(nsout-reheat,1); p_6i  = zeros(nsout-reheat,1);
        s_6i  = zeros(nsout-reheat,1); e_6i  = zeros(nsout-reheat,1);
        x_6i  = zeros(nsout-reheat,1); st_6i = cell(1,nsout-reheat);

        h_6is = h_3 + (h_6i - h_3(end))/eta_SiT2; % eq2.9
        for i = 1:nsout-reheat
            p_6i(i) = soutirage(s_3,h_3,h_6i(i),p_3,p_6,eta_SiT1);
            T_6i(i) = XSteam('T_ph',p_6i(i),h_6i(i));
            s_6i(i) = XSteam('s_ph',p_6i(i),h_6i(i));
            e_6i(i) = (h_6i(i)-h_0) - (T_0+273.15)*(s_6i(i)-s_0);
            if abs(T_6i(i) - XSteam('Tsat_p',p_6i(i))) < 1e-3
                x_6i(i) = XSteam('x_ph',p_6i(i),h_6i(i));
            else
                x_6i(i) = NaN;
            end
            c = convertCharsToStrings(num2roman(nsout-i-reheat+1));
            st_6i(i) = cellstr(strcat("6_",c));
        end
    else
        if nsout - reheat == 0
            p_6i = []; T_6i = []; h_6i = []; s_6i = []; e_6i = []; x_6i = []; st_6i = {};
        elseif nsout - reheat >= 1
            h_6i  = linspace(h_5(end),h_6,nsout-reheat+2)'; h_6i = h_6i(2:end-1);
            T_6i  = zeros(nsout-reheat,1); p_6i  = zeros(nsout-reheat,1);
            s_6i  = zeros(nsout-reheat,1); e_6i  = zeros(nsout-reheat,1);
            x_6i  = zeros(nsout-reheat,1); st_6i = cell(1,nsout-reheat);

            h_6is = h_5(end) + (h_6i - h_5(end))/eta_SiT2; % eq2.9
            for i = 1:nsout-reheat
                p_6i(i) = soutirage(s_5(end),h_5(end),h_6i(i),p_5(end),p_6,eta_SiT2);
                T_6i(i) = XSteam('T_ph',p_6i(i),h_6i(i));
                s_6i(i) = XSteam('s_ph',p_6i(i),h_6i(i));
                e_6i(i) = (h_6i(i)-h_0) - (T_0+273.15)*(s_6i(i)-s_0);
                if abs(T_6i(i) - XSteam('Tsat_p',p_6i(i))) < 1e-3
                    x_6i(i) = XSteam('x_ph',p_6i(i),h_6i(i));
                else
                    x_6i(i) = NaN;
                end
                c = convertCharsToStrings(num2roman(nsout-i-reheat+1));
                st_6i(i) = cellstr(strcat("6_",c));
            end
        end

        T_6i = [T_4';T_6i]; p_6i = [p_4';p_6i];
        h_6i = [h_4';h_6i]; s_6i = [s_4';s_6i];
        e_6i = [e_4';e_6i]; x_6i = [x_4';x_6i];
        for i = nsout-reheat+1:nsout
            c = convertCharsToStrings(num2roman(i));
            st_6i = cat(2,cellstr(strcat("6_",c)),st_6i);
        end
    end

% FRACTIONS OF TAPPING X_6i

    % STATES 7_i - DEFINE BY STATES 6_i AS AN ISOBARIC CONDENSATION. STATE 7_i ARE
    %               SATURATED LIQUID
    p_7i = zeros(nsout,1); h_7i = zeros(nsout,1);
    for i = 1:nsout
        p_7i(i) = p_6i(i) * lam_exch;
        h_7i(i) = XSteam('hL_p',p_7i(i));
    end

    % Temperature at states 9_i are determined by the temperature pinch at a heat
    % exchanger and the temparute of states 7_i
    T_9i = zeros(nsout,1);
    for i = 1:nsout
        T_9i(i) = XSteam('Tsat_p',p_7i(i)) - TpinchEx;
    end

    p_2 = p3_hp;

    if drumFlag
        drum = p_7i - XSteam('psat_T',Tdrum);
        drum(x_6i<=1) = NaN; drum(drum<=0) = NaN;
        drum = find(drum == min(drum));

        p_8 = XSteam('psat_T',Tdrum) * lam_exch^(nsout-drum);
        h_8 = h_7 + (p_8-p_7)*eta_SiC/10;
        T_8 = XSteam('T_ph',p_8,h_8);
        s_8 = XSteam('s_ph',p_8,h_8);
        e_8 = (h_8-h_0) - (T_0+273.15)*(s_8-s_0);
        x_8 = NaN;

        % State 7_0
        T_70 = T_8 + TpinchSub;
        h_70 = XSteam('h_pT',p_7i(end),T_70);

        p_7i(drum) = XSteam('psat_T',Tdrum);
        h_7i(drum) = XSteam('hL_T',Tdrum);
        T_7i(drum) = Tdrum;

        T_1 = T_9i(1);
        p_1 = XSteam('psat_T',T_1)*1.4251;
        h_1 = XSteam('h_pT',p_1,T_1);
        s_1 = XSteam('s_pT',p_1,T_1);
        e_1 = (h_1-h_0) - (T_0+273.15)*(s_1-s_0);
        x_1 = NaN;

        p_9i = zeros(nsout,1); h_9i = zeros(nsout,1);
        for i = 1:nsout
            if i <= drum
                p_9i(i) = p_1 / lam_exch^(i-1);
            else
                p_9i(i) = p_8 * lam_exch^(nsout-i+1);
            end
            h_9i(i) = XSteam('h_pT',p_9i(i),T_9i(i));
        end
        h_9i(drum) = h_7i(drum) + (p_9i(drum)-p_7i(drum))*eta_SiC/10;
        T_9i(drum) = XSteam('T_ph',p_9i(drum),h_9i(drum));

        T_90 = (T_8+T_9i(end))/2; p_90 = p_8;
        h_90 = XSteam('h_pT',p_90,T_90);
        iter = 0;
        while abs(h_90-iter) > 1
            [A,b] = sys_drumFlag(nsout,drum,h_6i,h_7i,h_9i,h_90);
            X_6i  = A\b; iter = h_90;
            h_90  = sum(X_6i(drum+2:end))*(h_7i(end)-h_70)/(1+sum(X_6i(drum+2:end))) + h_8;
        end
        T_90 = XSteam('T_ph',p_90,h_90);
        
        h_2 = h_1 + (p_2-p_1)*eta_SiC/10;
        T_2 = XSteam('T_ph',p_2,h_2);
        s_2 = XSteam('s_ph',p_2,h_2);
        e_2 = (h_2-h_0) - (T_0+273.15)*(s_2-s_0);
        x_2 = NaN;
    else
        T_1 = T_9i(1);
        p_1 = XSteam('psat_T',T_1)*1.4251;
        h_1 = XSteam('h_pT',p_1,T_1);
        s_1 = XSteam('s_pT',p_1,T_1);
        e_1 = (h_1-h_0) - (T_0+273.15)*(s_1-s_0);
        x_1 = NaN;

        p_8 = p_1 / lam_exch^(nsout);
        h_8 = h_7 + (p_8-p_7)*eta_SiC/10;
        T_8 = XSteam('T_ph',p_8,h_8);
        s_8 = XSteam('s_ph',p_8,h_8);
        e_8 = (h_8-h_0) - (T_0+273.15)*(s_8-s_0);
        x_8 = NaN;

        T_70 = T_8 + TpinchSub;
        h_70 = XSteam('h_pT',p_7i(end),T_70);

        p_9i = zeros(nsout,1); h_9i = zeros(nsout,1);
        for i = 1:nsout
            p_9i(i) = p_1 / lam_exch^(i-1);
            h_9i(i) = XSteam('h_pT',p_9i(i),T_9i(i));
        end

        T_90 = (T_8+T_9i(end))/2; p_90 = p_8;
        h_90 = XSteam('h_pT',p_90,T_90); iter = 0;
        while abs(h_90-iter) > 1
            [A,b] = sys_drumOut(nsout,h_6i,h_7i,h_9i,h_90);
            X_6i = A\b; iter = h_90;
            h_90 = sum(X_6i)*(h_7i(end)-h_70)/(1+sum(X_6i)) + h_8;
        end
        T_90 = XSteam('T_ph',p_1,h_90);
        
        h_2 = h_1 + (p_2-p_1)*eta_SiC/10;
        T_2 = XSteam('T_ph',p_2,h_2);
        s_2 = XSteam('s_ph',p_2,h_2);
        e_2 = (h_2-h_0) - (T_0+273.15)*(s_2-s_0);
        x_2 = NaN;
    end

% FRACTIONS OF ALL STATES
    X_6i = X_6i / (1 + sum(X_6i));
    X_1  = 1; X_2  = X_1; X_3 = X_1;
    X_4 = []; X_5 = [];
    if reheat > 0
        X_4 = X_1;
        X_5 = X_4 - X_6i(1);
        for i = 2:reheat
            X_4 = [X_4,X_5(i-1)];
            X_5 = [X_5,X_4(i) - X_6i(i)];
        end
    end
    X_6 = 1 - sum(X_6i);
    if drumFlag
        X_7 = X_6 + sum(X_6i(drum:end));
    else
        X_7 = X_1;
    end
    X_8 = X_7;

% STATES 7_i & 9_i KNOWING ALL ABOUT THE TAPPING
    p_7i = [p_7i;p_7i(end)]; h_7i = [h_7i;h_70];
    T_7i = zeros(nsout +1,1); s_7i = zeros(nsout +1,1);
    x_7i = zeros(nsout +1,1); X_7i = zeros(nsout +1,1);

    T_9i = [T_9i;T_90]; p_9i = [p_9i;p_90];
    h_9i = [h_9i;h_90]; s_9i = zeros(nsout +1,1);
    x_9i = zeros(nsout +1,1); X_9i = zeros(nsout +1,1);

    st_7i = cell(1,nsout +1); st_9i = cell(1,nsout +1);

    for i = 1:nsout +1
        T_7i(i) = XSteam('T_ph',p_7i(i),h_7i(i));
        s_7i(i) = XSteam('s_ph',p_7i(i),h_7i(i));
        h_9i(i) = XSteam('h_pT',p_9i(i),T_9i(i));
        s_9i(i) = XSteam('s_pT',p_9i(i),T_9i(i));

        if drumFlag
            if i < drum
                X_9i(i) = X_7;
                X_7i(i) = sum(X_6i(1:i));
            elseif i == drum
                X_9i(i) = X_1;
                X_7i(i) = X_1;
            else
                X_9i(i) = X_1;
                X_7i(i) = sum(X_6i(drum+1:min(i,nsout)));
            end
        else
            X_9i(i) = X_1;
            if i == nsout +1
                X_7i(i) = sum(X_6i);
            else
                X_7i(i) = sum(X_6i(1:i));
            end
        end

        if i == nsout +1
            x_7i(i) = NaN;
            st_7i(i) = {'7_0'};
            st_9i(i) = {'9_0'};
        else
            c = convertCharsToStrings(num2roman(nsout-i+1)); 
            st_7i(i) = cellstr(strcat("7_",c));
            st_9i(i) = cellstr(strcat("9_",c));
        end
        x_9i(i) = NaN;
    end
    e_7i = (h_7i-h_0) - (T_0+273.15)*(s_7i-s_0);
    e_9i = (h_9i-h_0) - (T_0+273.15)*(s_9i-s_0);


    RowNames = cat(2,{'1','2','3'},st_4,st_5,st_6i,{'6','7'},st_7i,{'8'},st_9i);
else
    s_2 = s_7; p_2 = p_3;% *1.129;
    h_2 = h_7 + (p_2 - p_7)/10;
    T_2 = XSteam('T_hs',h_2,s_2);
    e_2 = (h_2-h_0) - (T_0+273.15)*(s_2-s_0);
    x_2 = NaN;

    T_1 = [];T_6i = []; T_7i = []; T_8 = []; T_9i = [];
    p_1 = [];p_6i = []; p_7i = []; p_8 = []; p_9i = [];
    h_1 = [];h_6i = []; h_7i = []; h_8 = []; h_9i = [];
    s_1 = [];s_6i = []; s_7i = []; s_8 = []; s_9i = [];
    e_1 = [];e_6i = []; e_7i = []; e_8 = []; e_9i = [];
    x_1 = [];x_6i = []; x_7i = []; x_8 = []; x_9i = [];
    X_1 = [];X_6i = []; X_7i = []; X_8 = []; X_9i = [];
    if reheat == 0
        X_2 = 1; X_3 = 1; X_6 = 1; X_7 = 1;
        RowNames = {'2';'3';'4';'1'};
    else
        X_2 = 1; X_3 = 1; X_4 = 1; X_5 = 1; X_6 = 1; X_7 = 1;
        RowNames = {'2';'3';'4';'5';'6';'1'};
    end
end

state_1  = [T_1;p_1;h_1;s_1;e_1;x_1];
state_2  = [T_2;p_2;h_2;s_2;e_2;x_2];
state_8  = [T_8;p_8;h_8;s_8;e_8;x_8];
state_6i = [T_6i';p_6i';h_6i';s_6i';e_6i';x_6i'];
state_7i = [T_7i';p_7i';h_7i';s_7i';e_7i';x_7i'];
state_9i = [T_9i';p_9i';h_9i';s_9i';e_9i';x_9i'];

DAT = [state_1,state_2,state_3,state_4,state_5,state_6,state_7,state_8,...
    state_6i,state_7i,state_9i];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EFFICIENCY, MASSFLOWS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Works, Heat Exchanges
        P_e = P_e *1e3;
        P_m = P_e / eta_mec;


        W_mT_en = (h_3 + sum(h_5 - h_4) - h_6)*(1-sum(X_6i));
        W_mT_ex = (e_3 + sum(e_5 - e_4) - e_6)*(1-sum(X_6i));
        for i = 1:nsout
            if reheat == 0
                W_mT_en = W_mT_en + (h_3 - h_6i(i))*X_6i(i);
                W_mT_ex = W_mT_ex + (e_3 - e_6i(i))*X_6i(i);
            else
                if i == 1
                    W_mT_en = W_mT_en + (h_3 - h_4(1))*X_6i(1);
                    W_mT_ex = W_mT_ex + (e_3 - e_4(1))*X_6i(1);
                elseif i <= reheat
                    W_mT_en = W_mT_en + (h_3 - sum(h_4(1:i)) + sum(h_5(1:i-1)))*X_6i(i);
                    W_mT_ex = W_mT_ex + (e_3 - sum(e_4(1:i)) + sum(e_5(1:i-1)))*X_6i(i);
                else
                    W_mT_en = W_mT_en + (h_3 + sum(h_5-h_4) - h_6i(i))*X_6i(i);
                    W_mT_ex = W_mT_ex + (e_3 + sum(e_5-e_4) - e_6i(i))*X_6i(i);
                end
            end
        end

        W_mC_en = (h_2-h_1) + (h_8-h_7)*X_7;
        W_mC_ex = (e_2-e_1) + (e_8-e_7)*X_7;
        if drumFlag
            W_mC_en = W_mC_en + (h_9i(drum +1) - h_7i(drum +1))*X_7i(drum +1);
            W_mC_ex = W_mC_ex + (e_9i(drum +1) - e_7i(drum +1))*X_7i(drum +1);
        end

        W_m_en  = W_mT_en - W_mC_en;
        W_m_ex  = W_mT_ex - W_mC_ex;

        Q_I = (h_3-h_2);
        n_rh = 1; E_I = (e_3-e_2);
        while n_rh < reheat
            Q_I  = Q_I  + (h_5(n_rh) - h_4(n_rh))*(1-sum(X_6i(1:n_rh)));
            E_I  = E_I  + (e_5(n_rh) - e_4(n_rh))*(1-sum(X_6i(1:n_rh)));
            n_rh = n_rh +1;
        end
        Q_II = h_6*X_6 + h_7i(end)*X_7i(end) - h_7*X_7;

    if drumFlag
        eta_cyclen  = W_m_en / (Q_I + e_2-e_9i(drum))
    else
        eta_cyclen  = W_m_en / (Q_I + e_2-e_8);
    end
    eta_gen     = 1 - 1e-2 - eps_ech
    eta_toten   = eta_mec * eta_cyclen * eta_gen

    m_v = P_e / eta_mec / W_m_en
    m_c = m_v * Q_I / (eta_gen * LHV)
    m_a = m_ac * m_c
    m_f = m_a  / m_ag
    m_w = m_v * Q_II / (h_o-h_i)

    eta_rotex   = W_m_en / W_m_ex

    eta_condex  = (e_6-e_7) / e_6

    eta_gex     = 1/f * eta_gen * E_I/Q_I
    eta_chemex  = (e_h-e_e) / e_h
    eta_combex  = (m_f*(e_h-.04)) / (m_c*e_c);
    eta_transex = m_v*E_I / (m_f*(e_h-e_e));

    if drumFlag
        eta_cyclex  = W_m_en / (e_3-e_9i(drum+1))
    else
        eta_cyclex  = W_m_en / (e_3-e_1);
    end
    eta_totex   = eta_mec * eta_gex * eta_cyclex
    
    e_c = P_e / m_c / eta_totex;
    eta_combex  = (m_f*e_h) / (m_c*e_c)
    eta_transex = eta_gex / eta_chemex / eta_combex

    ETA = [eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_gen, ...
            eta_gex,eta_combex,eta_chemex,eta_condex,eta_transex];

    Xmassflow = flip(X_6i') * m_v;

    MASSFLOW = [m_a,m_v,m_c,m_f];


    COMBUSTION = struct();
        COMBUSTION.LHV = LHV;
        COMBUSTION.e_c = e_c;
        COMBUSTION.lambda = lambda;
        COMBUSTION.Cp_g = Cp_f(0,T_e);
        COMBUSTION.fum = m_f * [comp_f_O2,comp_f_N2,comp_f_CO2,comp_f_H2O];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    perte_mec = (P_m - P_e);

    perte_gen    = (h_h - h_e)*m_f * (1-eta_gen); %(h_h - h_e)*m_f;
    perte_conden = (h_o - h_i)*m_w;

    DATEN = [perte_gen,perte_mec,perte_conden];

    perte_totex  = m_c*e_c * (1-eta_totex);
    perte_rotex  = W_m_en*m_v * (1-eta_rotex);
    perte_combex = m_c*e_c * (1-eta_combex);
    perte_condex = (e_o - e_i)*m_w *(1-eta_condex);
    perte_chemex = (e_e/m_ag - ea_0)*m_a;
    perte_transex = m_f*(e_h-e_e) * (1-eta_transex);

    DATEX = [perte_mec,perte_totex,perte_rotex, ... 
        perte_combex,perte_condex,perte_chemex,perte_transex];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FIG = [];

if display
    % Table of values
        T = table([T_1;T_2;T_3;T_4';T_5';T_6i;T_6;T_7;T_7i;T_8;T_9i],...
            [p_1;p_2;p_3;p_4';p_5';p_6i;p_6;p_7;p_7i;p_8;p_9i],...
            [h_1;h_2;h_3;h_4';h_5';h_6i;h_6;h_7;h_7i;h_8;h_9i],...
            [s_1;s_2;s_3;s_4';s_5';s_6i;s_6;s_7;s_7i;s_8;s_9i],...
            [e_1;e_2;e_3;e_4';e_5';e_6i;e_6;e_7;e_7i;e_8;e_9i],...
            [x_1;x_2;x_3;x_4';x_5';x_6i;x_6;x_7 ;x_7i;x_8;x_9i],...
            [X_1;X_2;X_3;X_4';X_5';X_6i;X_6;X_7;X_7i;X_8;X_9i]*m_v,...
            'RowNames',RowNames);
        T.Properties.VariableNames = {'Temperature','Pressure',...
            'Enthalpy','Entropy','Exergy','Ratio','Debit'};
        disp(T)

    samp = 20; lw = .75;
    %TS - graph
        FIG(1) = figure('Color','w'); hold on; grid;
        xlabel('entropie, $s$ [kJ/(kg K)]','Interpreter','latex');
        ylabel('temperature, $T$ [$^\circ$C]','Interpreter','latex');

        T = linspace(0,373.9458,1e2);
        sL = []; sV = [];
        hL = []; hV = [];
        for t = T
            sL = [sL,XSteam('sL_T',t)]; sV = [sV,XSteam('sV_T',t)];
            hL = [hL,XSteam('hL_T',t)]; hV = [hV,XSteam('hV_T',t)];
        end
        plot([sL,flip(sV)],[T,flip(T)],'k-.');

        T_23 = linspace(T_2,T_3,samp)'; p_23 = linspace(p_2,p_3,samp)';
        s_23 = zeros(samp,1); h_23 = zeros(samp,1);
        T_34 = zeros(samp,1); s_34 = zeros(samp,1);
        h_34 = linspace(h_3,h_4(1),samp)';
        T_56 = zeros(samp,1); s_56 = zeros(samp,1);
        h_56 = linspace(h_5(end),h_6,samp)';
        for i = 1:samp
            s_23(i) = XSteam('s_pT',p_23(i),T_23(i));
            h_23(i) = XSteam('h_pT',p_23(i),T_23(i));

            p_34 = soutirage(s_3,h_3,h_34(i),p_3,p_4(1),eta_SiT1);
            T_34(i) = XSteam('T_ph',p_34,h_34(i));
            s_34(i) = XSteam('s_ph',p_34,h_34(i));

            p_56 = soutirage(s_5(end),h_5(end),h_56(i),p_5(end),p_6,eta_SiT2);
            T_56(i) = XSteam('T_ph',p_56,h_56(i));
            s_56(i) = XSteam('s_ph',p_56,h_56(i));
        end
        my_fig = plot([s_56;s_6;s_7;s_8;flip(s_9i);s_1;s_23;s_34],...
            [T_56;T_6;T_7;T_8;flip(T_9i);T_1;T_23;T_34]);
        my_fig.Color = [0.8500 0.3250 0.0980]; my_fig.LineWidth = lw;

        for j = 1:reheat-1
            T_45 = linspace(T_4(j),T_5(j),samp); p_45 = linspace(p_4(j),p_5(j),samp); s_45 = zeros(1,samp);
            T_54 = zeros(1,samp); h_54 = linspace(h_5(j),h_4(j+1),samp); s_54 = zeros(1,samp);
            for i = 1:samp
                s_45(i) = XSteam('s_pT',p_45(i),T_45(i));

                p_54 = soutirage(s_5(j),h_5(j),h_54(i),p_5(j),p_4(j+1),eta_SiT2);
                T_54(i) = XSteam('T_ph',p_54,h_54(i));
                s_54(i) = XSteam('s_ph',p_54,h_54(i));
            end
            my_fig = plot([s_45,s_54],[T_45,T_54]);
            my_fig.Color = [0.8500 0.3250 0.0980]; my_fig.LineWidth = lw;
        end

        T_45 = linspace(T_4(end),T_5(end),samp)'; p_45 = linspace(p_4(end),p_5(end),samp)'; s_45 = zeros(samp,1);
        for i = 1:samp
            s_45(i) = XSteam('s_pT',p_45(i),T_45(i));
        end
        my_fig = plot(s_45,T_45); my_fig.Color = [0.8500 0.3250 0.0980]; my_fig.LineWidth = lw;

        my_fig = plot([s_1;s_2;s_3;s_4';s_5';s_6;s_7;s_8],...
            [T_1;T_2;T_3;T_4';T_5';T_6;T_7;T_8]);
        my_fig.Marker = 'o'; my_fig.MarkerSize = 3; my_fig.LineStyle = 'none';
        my_fig.Color = [0.8500 0.3250 0.0980]; my_fig.LineWidth = lw;

        for i = 1:nsout
            if i == nsout
                my_fig = plot([s_6i(i),s_7i(i),s_7i(i+1),s_8],[T_6i(i),T_7i(i),T_7i(i+1),T_8]);
            elseif drumFlag && i == drum
                p_67 = ones(1,samp)*p_7i(i); T_67 = zeros(1,samp);
                s_67 = linspace(XSteam('s_ph',p_7i(i),h_6i(i)),XSteam('sV_T',T_7i(i)),samp);
                for j = 1:samp
                    T_67(j) = XSteam('T_ps',p_67(j),s_67(j));
                end
                my_fig = plot([s_6i(i),s_67,s_7i(i),adjust(T_7i(i+1),h_7i(i),0)],...
                    [T_6i(i),T_67,T_7i(i),T_7i(i+1)]);
            elseif x_6i(i) <= 1
                my_fig = plot([s_6i(i),s_7i(i),adjust(T_7i(i+1),h_7i(i),0)],...
                    [T_6i(i),T_7i(i),T_7i(i+1)]);
            else
                T_67 = zeros(1,samp); p_67 = linspace(p_6i(i),p_7i(i),samp);
                s_67 = linspace(s_6i(i),XSteam('sV_T',T_7i(i)),samp);
                for j = 1:samp
                    T_67(j) = XSteam('T_ps',p_67(j),s_67(j));
                end
                my_fig = plot([s_67,s_7i(i),adjust(T_7i(i+1),h_7i(i),0)],...
                    [T_67,T_7i(i),T_7i(i+1)]);
            end
            my_fig.Color = [0.8500 0.3250 0.0980];
            my_fig.LineStyle = '--';
            if drumFlag
                if i == drum
                    my_fig.LineWidth = lw;
                end
            end
        end

    %sh - graph
        FIG(2) = figure('Color','w'); hold on; grid;
        xlabel('entropie, $s$ [kJ/(kg K)]','Interpreter','latex');
        ylabel('enthalpie, $h$ [kJ/kg]','Interpreter','latex');

        plot([sL,flip(sV)],[hL,flip(hV)],'k-.');

        my_fig = plot([s_56;s_7;s_8;flip(s_9i);s_1;s_23;s_34],...
            [h_56;h_7;h_8;flip(h_9i);h_1;h_23;h_34]);
        my_fig.Color = [0.8500 0.3250 0.0980]; my_fig.LineWidth = lw;

        for j = 1:reheat-1
            T_45 = linspace(T_4(j),T_5(j),samp); p_45 = linspace(p_4(j),p_5(j),samp);
            h_45 = zeros(1,samp); s_45 = zeros(1,samp);
            h_54 = linspace(h_5(j),h_4(j+1),samp); s_54 = zeros(1,samp);
            for i = 1:samp
                h_45(i) = XSteam('h_pT',p_45(i),T_45(i));
                s_45(i) = XSteam('s_pT',p_45(i),T_45(i));

                p_54 = soutirage(s_5(j),h_5(j),h_54(i),p_5(j),p_4(j+1),eta_SiT2);
                s_54(i) = XSteam('s_ph',p_54,h_54(i));
            end
            my_fig = plot([s_45,s_54],[h_45,h_54]);
            my_fig.Color = [0.8500 0.3250 0.0980]; my_fig.LineWidth = lw;
        end

        T_45 = linspace(T_4(end),T_5(end),samp)'; p_45 = linspace(p_4(end),p_5(end),samp)';
        h_45 = zeros(1,samp); s_45 = zeros(1,samp);
        for i = 1:samp
            h_45(i) = XSteam('h_pT',p_45(i),T_45(i));
            s_45(i) = XSteam('s_pT',p_45(i),T_45(i));
        end
        my_fig = plot(s_45,h_45); my_fig.Color = [0.8500 0.3250 0.0980]; my_fig.LineWidth = lw;

        my_fig = plot([s_1;s_2;s_3;s_4';s_5';s_6;s_7;s_8],...
            [h_1;h_2;h_3;h_4';h_5';h_6;h_7;h_8]);
        my_fig.Marker = 'o'; my_fig.MarkerSize = 3; my_fig.LineStyle = 'none';
        my_fig.Color = [0.8500 0.3250 0.0980]; my_fig.LineWidth = lw;

        for i = 1:nsout
            if i == nsout
                my_fig = plot([s_6i(i),s_7i(i),s_7i(i+1),s_8],[h_6i(i),h_7i(i),h_7i(i+1),h_8]);
            elseif x_6i(i) <= 1
                my_fig = plot([s_6i(i),s_7i(i),adjust(T_7i(i+1),h_7i(i),0)],...
                    [h_6i(i),h_7i(i),h_7i(i)]);
            else
                h_67 = zeros(1,samp); p_67 = linspace(p_6i(i),p_7i(i),samp);
                s_67 = linspace(s_6i(i),XSteam('sV_T',T_7i(i)),samp);
                for j = 1:samp
                    h_67(j) = XSteam('h_ps',p_67(j),s_67(j));
                end
                my_fig = plot([s_67,s_7i(i),adjust(T_7i(i+1),h_7i(i),0)],...
                    [h_67,h_7i(i),h_7i(i)]);
            end
            my_fig.Color = [0.8500 0.3250 0.0980];
            my_fig.LineStyle = '--';
            if drumFlag
                if i == drum
                    my_fig.LineWidth = lw;
                end
            end
        end

    % pie chart - energy
        FIG(3) = figure('Color','w'); grid;
        X = [P_e DATEN];
        pie(X,{'','','',''});
        legend(sprintf('Puissance effective : %g [MW]',P_e*1e-3),...
            sprintf('Pertes de combustion : %.4g [MW]',DATEN(1)*1e-3),...
            sprintf('Pertes mecaniques : %.3g [MW]',DATEN(2)*1e-3),...
            sprintf('Pertes au condenseur : %.5g [MW]',DATEN(3)*1e-3),...
            'Location','southoutside');

    % pie chart - energy
        FIG(4) = figure('Color','w'); grid;
        Irr_van = abs(X_7i(end)*(h_70-h_0 + T_0*(XSteam('s_ph',p_7,h_70)-s_0)));
        for i = 1:nsout
            e = (h_7i(i)-h_0) + T_0*(XSteam('s_ph',p_7i(i+1),h_7i(i))-s_0);
            Irr_van = Irr_van + abs(X_7i(i)*(e - e_7i(i)));
        end
        if drumFlag
            Irr_drum = m_v*(X_9i(drum+1)*e_9i(drum+1) + X_7i(drum-1)*e_7i(drum-1) ...
                + X_6i(drum)*e_6i(drum) - e_7i(drum));

            X = [P_e DATEX(1) DATEX(3:end),Irr_drum,Irr_van];
            pie(X,{'','','','','','','','',''});
            legend(sprintf('Puissance effective : %g [MW]',P_e*1e-3),...
                sprintf('Pertes mecaniques : %.3g [MW]',DATEX(1)*1e-3),...
                sprintf('Irr. aux turbines : %.5g [MW]',DATEX(3)*1e-3),...
                sprintf('Irr. a la combustion : %.5g [MW]',DATEX(4)*1e-3),...
                sprintf('Irr. au condenseur : %.5g [MW]',DATEX(5)*1e-3),...
                sprintf('Pertes a la cheminee : %.5g [MW]',DATEX(6)*1e-3),...
                sprintf('Irr. aux echangeurs : %.5g [MW]',DATEX(7)*1e-3),...
                sprintf('Irr. a la bache : %5g [MW]',Irr_drum*1e-3),...
                sprintf('Irr. aux vannes : %5g [MW]',Irr_van*1e-3),...
                'Location','southoutside');
        else
            X = [P_e DATEX(1) DATEX(3:end),Irr_van];
            pie(X,{'','','','','','','',''});
            legend(sprintf('Puissance effective : %g [MW]',P_e*1e-3),...
                sprintf('Pertes mecaniques : %.3g [MW]',DATEX(1)*1e-3),...
                sprintf('Irr. aux turbines : %.5g [MW]',DATEX(3)*1e-3),...
                sprintf('Irr. a la combustion : %.5g [MW]',DATEX(4)*1e-3),...
                sprintf('Irr. au condenseur : %.5g [MW]',DATEX(5)*1e-3),...
                sprintf('Pertes a la cheminee : %.5g [MW]',DATEX(6)*1e-3),...
                sprintf('Irr. aux echangeurs : %.5g [MW]',DATEX(7)*1e-3),...
                sprintf('Irr. aux vannes : %5g [MW]',Irr_van*1e-3),...
                'Location','southoutside');
        end
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p] = soutirage(s0,h0,hf,p0,pf,eta)
    tol = 1e-3; nmax = 5e2;
    err = 1; n = 0;
    while (abs(err) > tol && n < nmax)
        p = (p0 + pf)/2;
        hs  = XSteam('h_ps',p,s0);
        err = hf - h0 + eta * (h0 - hs); % eq2.9
        if err > 0
            pf = p;
        else
            p0 = p;
        end
        n = n +1;
    end
end

function [A,b] = sys_drumFlag(nsout,drum,h_6i,h_7i,h_9i,h_90)
    A = zeros(nsout); b = zeros(nsout,1);

    % Premiers soutirages [Condenseur ; Drum]
    if drum < nsout
        for i = drum+1:nsout
            A(i,drum+1:i) = h_6i(drum+1:i) - h_7i(i);
            if i == nsout
                dh = h_9i(drum +1) - h_90;
            else
                dh = h_9i(drum +1) - h_9i(i +1);
            end
            A(i,drum+1:end) = A(i,drum+1:end) - dh;
            b(i) = dh;
        end
        A(drum,drum+1:end) = h_7i(drum) - h_9i(drum +1);
        b(drum) = -h_7i(drum) + h_9i(drum +1);
    else
        b(drum) = -h_7i(drum) + h_90;
    end

    % Soutirage [Drum]
    A(drum,drum) = h_7i(drum) - h_6i(drum);

    % Soutirages [Drum ; Combustion]
    if drum > 1
        for i = 1:drum-1
            A(i,1:i) = h_7i(i) - h_6i(1:i);
            dh = h_9i(1) - h_9i(i +1);
            A(i,:) = A(i,:) + dh;
            b(i) = -dh;
        end
        A(drum,1:drum-1) = h_7i(drum) - h_7i(drum -1);
    end
end

function [A,b] = sys_drumOut(nsout,h_6i,h_7i,h_9i,h_90)
    A = zeros(nsout); b = zeros(nsout,1);

    for i = 1:nsout
        A(i,1:i) = h_6i(1:i) - h_7i(1);
        if i == nsout
            dh = h_9i(1) - h_90;
        else
            dh = h_9i(1) - h_9i(i+1);
        end
        A(i,:) = A(i,:) - dh; b(i) = dh;
    end
end

