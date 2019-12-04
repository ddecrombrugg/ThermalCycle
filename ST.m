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


%% PARAMETERS VERIFICATION 

    if nargin<3
        display = 1;
        if nargin<2
            constr = 'pressure';
            %constr = 'vap_ratio';
            options = struct();
            if nargin<1
                P_e = 35e3; % [kW] Puissance énergétique de l'installation
            end
        end
    end

    % options.nsout [-] : Number of feed-heating
    if isfield(options,'nsout') %OK
        nsout = options.nsout;
    else
        nsout = 8;
    end

    % options.reheat [-] : Number of reheating
    if isfield(options,'reheat') %OK
        reheat = options.reheat;
    else
        reheat = 1;
    end

    % options.T_max [°C] : Maximum steam temperature
    if isfield(options,'T_max') %OK
        T_max = options.T_max;
    else
        T_max = 565;
    end

    % options.T_cond_out [°C] : Condenseur cold outlet temperature
    if isfield(options,'T_cond_out') %OK
        T_cond_out = options.T_cond_out;
    else
        T_cond_out = 26;
    end

    % options.p3_hp [bar] : Maximum pressure
    if isfield(options,'p3_hp') %OK
        p3_hp = options.p3_hp;
    else
        p3_hp = 310;
    end

    % options.drumFlag [-] : if =1 then drum if =0 => no drum.
    if isfield(options,'drumFlag')
        drumFlag = options.drumFlag;
    else
        drumFlag = 1;
    end

    %   -options.eta_mec [-] : mecanic efficiency of shafts bearings
    if isfield(options,'eta_mec')
        eta_mec = options.eta_mec;
    else
        eta_mec = .98;
    end

    % options.comb is a structure containing combustion data : 
    %   -comb.Tmax     [°C] : maximum combustion temperature
    %   -comb.lambda   [-] : air excess
    %   -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
    %   -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
    if isfield(options,'comb')
        Tmax = options.comb.Tmax;
        lambda = options.comb.lambda;
        x = options.comb.x;
        y = options.comb.y;
    else
        Tmax = 1000;
        lambda = 1.05;
        x = 0;
        y = 4;
    end

    % options.T_exhaust [°C] : Temperature of exhaust gas out of the chimney
    if isfield(options,'T_exhaust')
        T_exhaust = options.T_exhaust;
    else
        T_exhaust = 120;
    end

    % options.p4 [bar] : High pressure after last reheating
    % options.x6 [-] : Vapor ratio [gaseous/liquid] (in french : titre)
    if isfield(options,'p_4') && ~isfield(options,'x_6')
        p_4 = options.p_4;
        constr = 'pressure';
    elseif isfield(options,'x_6')
        x_6 = options.x_6;
        constr = 'vap_ratio';
    elseif constr == 'pressure'
        p_4 = 70;
    else
        x_6 = .8633;
    end

    % options.T_0 [°C] : Reference temperature
    if isfield(options,'T_0')
        T_0 = options.T_0;
    else
        T_0 = 15;
    end

    % options.TpinchSub [°C] : Temperature pinch at the subcooler
    if isfield(options,'TpinchSub')
        TpinchSub = options.TpinchSub
    else
        TpinchSub = 30;
    end


    % options.TpinchEx [°C] : Temperature pinch at a heat exchanger
    if isfield(options,'TpinchEx') %OK
        TpinchEx = options.TpinchEx;
    else
        TpinchEx = 20;
    end

    % options.TpinchCond [°C] : Temperature pinch at condenser
    if isfield(options,'TpinchCond') %OK
        TpinchCond = options.TpinchCond;
    else
        TpinchCond = 7;
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
        eta_SiT = .927;  % [-]
    end

    if length(eta_SiT) == 2
        eta_SiT1 = eta_SiT(1);
        eta_SiT2 = eta_SiT(2);
    else
        eta_SiT1 = eta_SiT;
        eta_SiT2 = eta_SiT;
    end


%% OTHER PARAMETERS 

    M_O2  = 31.99800e-3;  % [kg/mol]
    M_N2  = 28.01400e-3;  % [kg/mol]
    M_CO2 = 44.00800e-3;  % [kg/mol]
    M_H2O = 18.01494e-3;  % [kg/mol]
    M_air = .21*M_O2 + .79*M_N2;  % [kg/mol]

    p_ext = 1.01325;  % [bar]

    R = 8.314472e-3;  % The ideal gas's constant [kJ/mol/K]
    R_O2  = R / M_O2;   % [kJ/(kg*K)]
    R_N2  = R / M_N2;   % [kJ/(kg*K)]
    R_CO2 = R / M_CO2;  % [kJ/(kg*K)]
    R_H2O = R / M_H2O;  % [kJ/(kg*K)]
    R_air = 287.058e-3; % [kJ/(kg*K)]

    function [X] = Cp_CO2(T)
        T = T +273.15;
        T(T<300) = 300; T(T>5000) = 5000;
        X = janaf('CO2',T); % [kJ/(kg*K)]
    end

    function [X] = Cp_H2O(T)
        T = T +273.15;
        T(T<300) = 300; T(T>5000) = 5000;
        X = janaf('H2O',T);
    end

    function [X] = Cp_O2(T)
        T = T +273.15;
        T(T<300) = 300; T(T>5000) = 5000;
        X = janaf('O2',T);
    end

    function [X] = Cp_N2(T)
        T = T +273.15;
        T(T<300) = 300; T(T>5000) = 5000;
        X = janaf('N2',T);
    end

    function [X] = Cp_air(T)
        T = T +273.15;
        T(T<300) = 300; T(T>5000) = 5000;
        X = .21*janaf('O2',T) + .79*janaf('N2',T);
    end

    function [X] = Cp_vap(p,T)
        X = zeros(size(p));
        for i = 1:length(p(:,1))
            for j = 1:length(p(i,:))
                X(i,j) = XSteam('Cp_pT',p(i,j),T(i,j));
            end
        end
    end


%% REFERENCE STATE 

    p_0 = XSteam('psat_T',T_0);
    h_0 = XSteam('hV_T',T_0);
    s_0 = XSteam('sV_T',T_0);
    e_0 = 0;


%% COMBUSTION 

    a = (lambda-1)*(1+y/4-x/2); b = y/2;
    w = lambda*(1+y/4-x/2); % Stoechiometric coefficients

    M_c = (12.01+1.01*y+16*x)*1e-3; % Molar mass of carburant [kg/mol]
    LHV = (393.6 + 102.2*y - (110.6 + 204.4*y)*x/(1+y/2)) / M_c; % [kJ/kg]

    comp_f_tot = M_CO2 + b*M_H2O + a*M_O2 + 3.76*w*M_N2; % [kg]
    comp_f_CO2 = M_CO2 / comp_f_tot; % [-]
    comp_f_H2O = b*M_H2O / comp_f_tot; % [-]
    comp_f_O2  = a*M_O2 / comp_f_tot; % [-]
    comp_f_N2  = 3.76*w*M_N2 / comp_f_tot; % [-]
    R_f = R / comp_f_tot * (1+b+a+3.76*w); % [kJ/(kg*K)]

    function [X] = Cp_f(T)
        T = T +273.15;
        T(T<300) = 300; T(T>5000) = 5000;
        X = comp_f_CO2*janaf('CO2',T) + comp_f_H2O*janaf('H2O',T) + ...
            comp_f_O2*janaf('O2',T) + comp_f_N2*janaf('N2',T); % [kJ/(kg*K)]
    end

    T_c1 = Tmax; % Temperature after combustion [K]
    h_c1 = h_0 + integral(@Cp_f,T_0,T_c1)*(T_c1 - T_0); % Enthalpy after combustion [kJ/kg]
    T_ce = T_exhaust; % Temperature after the heat exchange [K]
    h_ce = h_0 + integral(@Cp_f,T_0,T_ce)*(T_ce - T_0); % Enthalpy after the heat exchange [kJ/kg]

    m_a1 = (1+y/4-x/2) * (M_air/.21)/M_c; % [mol]
    m_ac = lambda * m_a1;  % [-]
    m_ag = (1 + 1/m_ac)^(-1); % [-]
    Q_comb = LHV / m_ac;


%% CYCLE

T_3 = T_max; % given
p_3 = p3_hp; % given
h_3 = XSteam('h_pT',p_3,T_3);
s_3 = XSteam('s_pT',p_3,T_3);
e_3 = (h_3-h_0) - (T_0+273.15)*(s_3-s_0);
v_3 = XSteam('v_pT',p_3,T_3);
state_3 = [T_3;p_3;h_3;s_3;e_3;NaN];

T_7 = T_cond_out + TpinchCond;
p_7 = XSteam('psat_T',T_7);
h_7 = XSteam('hL_T',T_7);
s_7 = XSteam('sL_T',T_7);
e_7 = (h_7-h_0) - (T_0+273.15)*(s_7-s_0);
state_7 = [T_7;p_7;h_7;s_7;e_7;.0];


if constr == 'pressure'
    h_4s = XSteam('h_ps',p_4,s_3);
    h_4  = h_3 - eta_SiT1 * (h_3 - h_4s); % eq2.9
    T_4  = XSteam('T_ph',p_4,h_4);
    s_4  = XSteam('s_ph',p_4,h_4);
    e_4 = (h_4-h_0) - (T_0+273.15)*(s_4-s_0);
    state_4 = [T_4;p_4;h_4;s_4;e_4;NaN];

    p_5 = p_4 /1.1; % isobare
    T_5 = T_max; % réchauffeur
    h_5 = XSteam('h_pT',p_5,T_5);
    s_5 = XSteam('s_pT',p_5,T_5);
    e_5 = (h_5-h_0) - (T_0+273.15)*(s_5-s_0);
    state_5 = [T_5;p_5;h_5;s_5;e_5;NaN];

    T_6  = T_7;
    p_6  = XSteam('psat_T',T_6);
    h_6s = XSteam('h_ps',p_6,s_5);
    h_6  = h_5 - eta_SiT2 * (h_5 - h_6s); % eq2.9
    s_6  = XSteam('s_ph',p_6,h_6);
    e_6  = (h_6-h_0) - (T_0+273.15)*(s_6-s_0);
    x_6  = XSteam('x_ph',p_6,h_6);
    v_6  = XSteam('v_ph',p_6,h_6);
    state_6 = [T_6;p_6;h_6;s_6;e_6;x_6];
else
    T_6 = T_7;
    p_6 = XSteam('psat_T',T_6);
    h_6 = XSteam('h_Tx',T_6,x_6);
    s_6 = XSteam('sL_T',T_6)*(1-x_6) + XSteam('sV_T',T_6)*x_6;
    e_6 = (h_6-h_0) - (T_0+273.15)*(s_6-s_0);
    v_6  = XSteam('v_ph',p_6,h_6);
    state_6 = [T_6;p_6;h_6;s_6;e_6;x_6];

    T_5 = T_3;
    s_5s = s_6;
    h_5s = adjust(T_5,[1e3,5e3],s_5s);
    h_5  = h_6 + eta_SiT*(h_5s-h_6);
    s_5  = adjust(T_5,h_5,[7,10]);
    p_5  = XSteam('p_hs',h_5,s_5);
    e_5 = (h_5-h_0) - (T_0+273.15)*(s_5-s_0);
    state_5 = [T_5;p_5;h_5;s_5;e_5;NaN];

    p_4  = p_5 *1.1;
    h_4s = XSteam('h_ps',p_4,s_3);
    h_4  = h_3 - eta_SiT1 * (h_3 - h_4s); % eq2.9
    T_4  = XSteam('T_ph',p_4,h_4);
    s_4  = XSteam('s_ph',p_4,h_4);
    e_4  = (h_4-h_0) - (T_0+273.15)*(s_4-s_0);
    state_4 = [T_4;p_4;h_4;s_4;e_4;NaN];
end


% Soutirages isoenthalpiques
if nsout == 1
    h_6i  = linspace(h_5,h_6,3)'; h_6i  = h_6i(2);
    h_6is = h_5 + (h_6i - h_5)/eta_SiT; % eq2.9
    p_6i  = soutirage(s_5,h_5,h_6i,p_5,p_6,eta_SiT2);
    T_6i  = XSteam('T_ph',p_6i,h_6is);
    s_6i  = XSteam('s_ph',p_6i,h_6i);
    e_6i  = (h_6i-h_0) - (T_0+273.15)*(s_6i-s_0);
    if abs(p_6i - XSteam('psat_T',T_6i)) < 1e-3
        x_6i = XSteam('x_ph',p_6i,h_6i);
    else
        x_6i = NaN;
    end
    st_6i = {'6_I'};
elseif nsout > 1
    h_6i  = linspace(h_5,h_6,nsout+1)'; h_6i = h_6i(2:end-1);
    h_6is = h_5 + (h_6i - h_5)/eta_SiT; % eq2.9
    T_6i  = zeros(nsout-1,1); p_6i  = zeros(nsout-1,1);
    s_6i  = zeros(nsout-1,1); e_6i  = zeros(nsout-1,1);
    x_6i  = zeros(nsout-1,1); st_6i = cell(1,nsout-1);
    for i = 1:nsout-1
        p_6i(i) = soutirage(s_5,h_5,h_6i(i),p_5,p_6,eta_SiT2);
        T_6i(i) = XSteam('T_ph',p_6i(i),h_6is(i));
        s_6i(i) = XSteam('s_ph',p_6i(i),h_6i(i));
        e_6i(i) = (h_6i(i)-h_0) - (T_0+273.15)*(s_6i(i)-s_0);
        if abs(T_6i(i) - XSteam('T_ph',p_6i(i),h_6i(i))) < 1e-3
            x_6i(i) = XSteam('x_ph',p_6i(i),h_6i(i));
        else
            x_6i(i) = NaN;
        end
        c = convertCharsToStrings(num2roman(nsout-i));
        st_6i(i) = cellstr(strcat("6_",c));
    end
end

p_1 = 100;
% Fractions des soutirages
if drumFlag
    T_8 = T_7;
    p_8 = XSteam('psat_T',Tdrum);
    h_8 = XSteam('h_pT',p_8,T_8);
    s_8 = XSteam('s_pT',p_8,T_8);
    e_8 = (h_8-h_0) - (T_0+273.15)*(s_8-s_0);
    state_8 = [T_8;p_8;h_8;s_8;e_8;NaN];

    eta_Ex = .98; drum = round(nsout/2); drum_r = nsout-drum;
    p_7i = [eta_Ex*p_4;zeros(nsout-1,1)];
    h_7i = [XSteam('hL_p',eta_Ex*p_4);zeros(nsout-1,1)];
    for i = 1:nsout-1
        p_7i(i+1) = p_6i(i) * eta_Ex;
        h_7i(i+1) = XSteam('hL_p',p_7i(i+1));
    end
    h_7i(drum +1) = p_8;
    h_7i(drum +1) = XSteam('hL_T',Tdrum);

    T_9i = zeros(nsout,1);
    for i = 1:nsout
        T_9i(i) = XSteam('Tsat_p',p_7i(i)) - TpinchEx;
    end
    T_9i(drum +1) = Tdrum;

    T_70 = T_7;
    h_70 = XSteam('h_pT',p_7i(end),T_70);
    T_90 = T_70 - TpinchSub;
    h_90 = XSteam('h_pT',p_8,T_90);

    A = zeros(nsout);
    b = zeros(nsout,1);

    % Premiers soutirages [Condenseur ; Drum]
    for i = drum+2:nsout
        A(i,drum+2:i) = h_6i(drum+1:i-1) - h_7i(i);
        if i == nsout
            dh = XSteam('h_pT',p_8,T_9i(drum+2)) - h_90;
        else
            dh = XSteam('h_pT',p_8,T_9i(drum+2)) - XSteam('h_pT',p_8,T_9i(i+1));
        end
        A(i,drum+2:end) = A(i,drum+2:end) - dh*ones(1,drum-1);
        b(i) = dh;
    end

    A(drum+1,drum+2:end) = h_7i(drum +1) - XSteam('h_pT',p_8,T_9i(drum+2));
    A(drum+1,drum+1) = h_7i(drum +1) - h_6i(drum);
    A(drum+1,1:drum) = h_7i(drum +1) - h_7i(drum);
    b(drum+1) = -h_7i(drum +1) + XSteam('h_pT',p_8,T_9i(drum+2));

    % Soutirages [Drum ; Combustion]
    for i = 1:drum
        if i == 1
            A(1,1) = h_4 - h_7i(1);
        else
            A(i,1:i) = [h_4;h_6i(1:i-1)] - h_7i(1);
        end
        dh = XSteam('h_pT',p_1,T_9i(1)) - XSteam('h_pT',p_1,T_9i(i+1));
        A(i,:) = A(i,:) - dh;
        b(i) = dh;
    end

    X_6i = A\b;

    T_1 = T_9i(1);
end
m_vap = 1+sum(X_6i);
m_vap_7 = 1+sum(X_6i(drum+2:end));

p_2 = p_3 *1.2;
h_2 = h_3 + (1 + sum(X_6i(1:end-1))) * (h_5 - h_4 - Q_comb)/(1 + sum(X_6i));
T_2 = XSteam('T_ph',p_2,h_2);
s_2 = XSteam('s_ph',p_2,h_2);
e_2 = (h_2-h_0) - (T_0+273.15)*(s_2-s_0);
state_2 = [T_2;p_2;h_2;s_2;e_2;NaN];

h_1 = XSteam('h_pT',p_1,T_1);
s_1 = XSteam('s_pT',p_1,T_1);
e_1 = (h_1-h_0) - (T_0+273.15)*(s_1-s_0);
state_1 = [T_1;p_1;h_1;s_1;e_1;NaN];

T = table([T_1;T_2;T_3;T_4;T_5;T_6i;T_6;T_7;T_8],...
        [p_1;p_2;p_3;p_4;p_5;p_6i;p_6;p_7;p_8],...
        [h_1;h_2;h_3;h_4;h_5;h_6i;h_6;h_7;h_8],...
        [s_1;s_2;s_3;s_4;s_5;s_6i;s_6;s_7;s_8],...
        [e_1;e_2;e_3;e_4;e_5;e_6i;e_6;e_7;e_8],...
        [NaN;NaN;NaN;NaN;NaN;x_6i;x_6;.0 ;NaN],...
        [m_vap;m_vap;m_vap;m_vap;m_vap-X_6i(end);...
        X_6i(2:end);1;m_vap_7;m_vap_7],...
        'RowNames',cat(2,{'1','2','3','4','5'},st_6i,{'6','7','8'}));
T.Properties.VariableNames = {'Temperature','Pressure',...
    'Enthalpy','Entropy','Exergy','Ratio','Debit'};
disp(T)

end

function [Y] = adjust(T,h,s)
    tol = 1e-3; nmax = 5e2;

    if length(h) == 2
        T1 = XSteam('T_hs',h(1),s);
        T2 = XSteam('T_hs',h(2),s);
        if (T1 < T && T2 < T)
            Y2 = max(h)*2; Y1 = max(h);
        elseif(T1 > T && T2 > T)
            Y1 = min(h)/2; Y2 = min(h);
        else
            Y1 = min(h); Y2 = max(h);
        end
    else
        T1 = XSteam('T_hs',h,s(1));
        T2 = XSteam('T_hs',h,s(2));
        if (T1 < T && T2 < T)
            Y1 = min(s)/2; Y2 = min(s);
        elseif(T1 > T && T2 > T)
            Y2 = max(s)*2; Y1 = max(s);
        else
            Y1 = min(s); Y2 = max(s);
        end
    end

    err = 1; n = 0;
    while (abs(err) > tol && n < nmax)
        Y = (Y1 + Y2)/2;
        if length(h) == 2
            err = T - XSteam('T_hs',Y,s);
        else
            err = T - XSteam('T_hs',h,Y);
        end
        if err > 0
            Y1 = Y;
        else
            Y2 = Y;
        end
        n = n +1;
    end
end

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


