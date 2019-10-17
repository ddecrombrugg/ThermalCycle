%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermal Cycle - Gas Turbine %
%        Louis Peeters        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ETA DATEN DATEX DAT MASSFLOW COMBUSTION Cp_g FIG] = GT(P_e,options,display)
% GT Gas turbine modelisation
% GT(P_e,options,display) compute the thermodynamics states for a Gas
% turbine based on several inputs (given in OPTION) and based on a given 
% electricity production P_e. It returns the main results. It can as well
% plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated) Refer to Fig 3.1 from reference book (in english)
% P_E = electrical power output target [kW]
% OPTIONS is a structure containing :
%   -options.k_mec [-] : Shaft losses 
%   -options.T_0   [°C] : Reference temperature
%   -options.T_ext [°C] : External temperature
%   -options.r     [-] : Compression ratio
%   -options.k_cc  [-] : Coefficient of pressure losses due to combustion
%                        chamber
%   -options.T_3   [°C] : Temperature after combustion (before turbine)
%   -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for compression
%   -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for expansion
%DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then the
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_rotex, compressor-turbine exergy efficiency
%   -eta(6) : eta_combex, Combustion exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% DATEN is a vector with : 
%   -daten(1) : perte_mec [kW]
%   -daten(2) : perte_ech [kW]
% DATEX is a vector with :
%   -datex(1) : perte_mec [kW]
%   -datex(2) : perte_rotex [kW]
%   -datex(3) : perte_combex [kW]
%   -datex(4) : perte_echex  [kW]
% DAT is a matrix containing :
% dat = {T_1       , T_2       , T_3       , T_4; [°C]
%        p_1       , p_2       , p_3       , p_4; [bar]
%        h_1       , h_2       , h_3       , h_4; [kJ/kg]
%        s_1       , s_2       , s_3       , s_4; [kJ/kg/K]
%        e_1       , e_2       , e_3       , e_4;};[kJ/kg]
% MASSFLOW is a vector containing : 
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_c, combustible massflow [kg/s] 
%   -massflow(3) = m_f, exhaust gas massflow [kg/s]
% 
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combustible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas at 400 K [kJ/kg/K]
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
%


%% Your Work

% Exemple of how to use 'nargin' to check your number of inputs
if nargin < 3
    display = 1;
   if nargin < 2
       options=struct();
       options.k_mec = 0.015;
       options.T_0 = 15;
       options.T_ext = 15;
       options.r = 18  ;
       options.k_cc = 0.95;
       options.T_3 = 1400;
       options.eta_PiC = 0.9;
       options.eta_PiT = 0.9;
       if nargin < 1
           P_e = 230e3; % 100MW
       end
   end
end


% Exemple of how to use (isfield' to check if an option has been given (or
% not)
if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 15;   
end

if isfield(options,'T_ext')
    T_ext = options.T_ext;
else
    T_ext = 15;
end
if isfield(options,'k_mec')
    k_mec = options.k_mec;
else
    k_mec = 0.015;
end
if isfield(options,'r')
    r = options.r;
else
    r = 18;
end
if isfield(options,'k_cc')
    k_cc = options.k_cc;
else
    k_cc = 0.95;
end
if isfield(options,'T_3')
    T_3 = options.T_3;
else
    T_3 = 1400;
end
if isfield(options,'eta_PiC')
    eta_PiC = options.eta_PiC;
else
    eta_PiC = 0.9;
end
if isfield(options,'eta_PiT')
    eta_PiT = options.eta_PiT;
else
    eta_PiT = 0.9;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle State 1 - Before Compression %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial State
T_1 = T_ext;
t_1 = T_ext + 273.15;
t_0=t_1;
p_1 = 1;

c_pa = 1.006;
h_1 = c_pa*T_1;
s_1 = 0.058;
e_1 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle State 2 - After Compression %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_2 = r*p_1;                                                                
R_sa = 287.058/1000;

% T_2 iteration
step = 1;
t_2 = 1674;
t_2_iter = 0;

while abs(t_2-t_2_iter)>0.01
    t_2_iter = t_2;
    T_c1 = 300.*ones(floor(t_1-273.15)+1,1);
    T_c2 = [300.*ones(floor(300-273.15)+1,1)' 301:step:t_2];
    c_pa_moy1 = (mean(janaf('N2',T_c1))*0.79 + mean(janaf('O2',T_c1))*0.21);
    c_pa_moy2 = (mean(janaf('N2',T_c2))*0.79 + mean(janaf('O2',T_c2)*0.21));
    t_2 = t_1*r^((1/eta_PiC)*(R_sa*((t_2-273.15)- T_1)/(c_pa_moy2*(t_2-273.15)-c_pa_moy1*T_1))); % cf. eq 3.19-3.22
end

T_2 = t_2-273.15;

% Computation of h_2, s_2, e_2
c_pa_moy12 = (c_pa_moy2*T_2-c_pa_moy1*T_1)/(T_2-T_1); % cf. eq 3.22
h_2 = h_1 + c_pa_moy12*(T_2-T_1); % formule enthalpie
s_2 = s_1 + (1-eta_PiC)*c_pa_moy12*log(t_2/t_1); % cf. eq 3.15
e_2 = (h_2-h_1)-(T_0+273.15)*(s_2-s_1); % cf. eq 1.19

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle State 3 - After Combustion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_3 = T_3 + 273.15;
p_3 = k_cc*p_2;

M_O2 = 31.998000000000000/1000;
M_CO2 = 44.008000000000000/1000;
M_H2 = 2.015940000000000/1000;
M_H2O = 18.014940000000000/1000;
M_N2 = 28.014000000000000/1000;
M_CO = 28.0090000000000/1000;
M_CH4 = 16.0430000000/1000;

x = 0;
y = 4;
l_hv = 51.5*10^3;

% Air excess coefficient
T_c2 = [300.*ones(floor(300-273.15)+1,1)' 301:step:t_2];
T_c3 = [300.*ones(floor(300-273.15)+1,1)' 301:step:t_3];

f = @(lambda) l_hv*M_CH4 + T_2*(lambda*(1+(y-2*x)/4)*M_O2*mean(janaf('O2',T_c2))+ lambda*(1+(y-2*x)/4)*3.76*M_N2*mean(janaf('N2',T_c2)))...
    - T_3*(M_CO2*mean(janaf('CO2',T_c3)) + y/2*M_H2O*mean(janaf('H2O',T_c3)) + (1+(y-2*x)/4)*(lambda-1)*M_O2*mean(janaf('O2',T_c3))...
    + (1+(y-2*x)/4)*3.76*lambda*M_N2*mean(janaf('N2',T_c3))); % cf. Bilan enthalpique de réaction

lambda = fsolve(f,1);

M_a = (32+3.76*28)*(1+(y-2*x)/4)/(12+y+16*x);
M_s = 1 + M_a;

m_ac = M_a*lambda; % cf. eq 3.7
m_fa = (1+ 1/(m_ac)); % cf. eq 3.8
Q_comb = l_hv/m_ac; % cf. eq 3.9

% Computation of h_3, s_3, e_3
h_3 = (Q_comb+h_2)*(1/m_fa); % cf. 3.6
c_3 = Q_comb/(T_3-T_2); % cf. Adiabatic combustion 
s_3 = c_3*log(t_3/t_2)+0.1; % cf. Entropy of polytropic isobaric tranformation
e_3 = (h_3-h_1)-(T_0+273.15)*(s_3-s_1); % cf. eq 1.19

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle State 4 - After decompression %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_4 = p_3/(k_cc*r);

% Flues composition
v_1r = 1;
v_2r = lambda*(1+(y-2*x)/4);
v_3r = 3.76*lambda*(1+(y-2*x)/4);

v_1p = 1;
v_2p = y/2;
v_3p = (lambda-1)*(1+(y-2*x)/4);
v_4p = 3.76*lambda*(1+(y-2*x)/4);

comp_f_CO2 = v_1p*M_CO2/(v_1p*M_CO2+v_2p*M_H2O+v_3p*M_O2+v_4p*M_N2);
comp_f_H2O = v_2p*M_H2O/(v_1p*M_CO2+v_2p*M_H2O+v_3p*M_O2+v_4p*M_N2);
comp_f_O2 = v_3p*M_O2/(v_1p*M_CO2+v_2p*M_H2O+v_3p*M_O2+v_4p*M_N2);
comp_f_N2 = v_4p*M_N2/(v_1p*M_CO2+v_2p*M_H2O+v_3p*M_O2+v_4p*M_N2);

R_f = (8.314/(comp_f_CO2*M_CO2+comp_f_H2O*M_H2O+comp_f_O2*M_O2+comp_f_N2*M_N2))/1000;

% T_4 iteration
step = 1;
t_4 = 1674;
t_4_iter = 0;

while abs(t_4-t_4_iter)>0.01
    t_4_iter = t_4;
    T_c3 = [300.*ones(floor(300-273.15)+1,1)' 301:step:t_3];
    T_c4 = [300.*ones(floor(300-273.15)+1,1)' 301:step:t_4];
    c_pa_moy3 = mean(janaf('CO2',T_c3)*comp_f_CO2 + janaf('H2O',T_c3)*comp_f_H2O...
              + janaf('O2',T_c3)*comp_f_O2 + janaf('N2',T_c3)*comp_f_N2);
    c_pa_moy4 = mean(janaf('CO2',T_c4)*comp_f_CO2 + janaf('H2O',T_c4)*comp_f_H2O...
              + janaf('O2',T_c4)*comp_f_O2 + janaf('N2',T_c4)*comp_f_N2);
    t_4 = t_3*r^(-eta_PiT*(R_f*((t_4-273.15)-T_3)/(c_pa_moy4*(t_4-273.15)-c_pa_moy3*T_3))); % cf. eq 3.23-3.24
end

T_4 = t_4-273.15;

% Computation of h_4, s_4, e_4
c_pa_moy34 = (c_pa_moy4*T_4-c_pa_moy3*T_3)/(T_4-T_3); % cf. eq 3.22
h_4 = h_3 + c_pa_moy34*(T_4-T_3); % formule enthalpie
s_4 = s_3 - (1-eta_PiT)/eta_PiT*c_pa_moy34*log(t_4/t_3); % cf. eq 3.16
e_4 = (h_4-h_1)-(T_0+273.15)*(s_4-s_1); % cf. eq 1.19

%%%%%%%%%%%%
% Massflow %
%%%%%%%%%%%%

m_a = P_e/((1-k_mec)*(m_fa*(h_3-h_4)-(h_2-h_1))); % cf. eq 3.30
m_f = m_fa*m_a;
m_c = m_a/m_ac;

%%%%%%%%%%%%%%
% Combustion %
%%%%%%%%%%%%%%

m_CO2 = comp_f_CO2*m_f;
m_H2O = comp_f_H2O*m_f;
m_O2 = comp_f_O2*m_f;
m_N2 = comp_f_N2*m_f;

hhv = 55695;
c_pc = 35.3/(1000*M_CH4);
S_CH4 = 183.1/(1000*M_CH4);
S_O2 = 202.8/(1000*M_O2);
S_CO2 = 210.4/(1000*M_CO2);
S_H2O = 69.5/(1000*M_H2O);

T_ec = 300.*ones(floor(t_0-273.15),1)';

e_c = hhv + 15*(c_pc + comp_f_O2/m_O2*mean(janaf('O2',T_ec))...
    -comp_f_CO2/m_CO2*mean(janaf('CO2',T_ec)))...
    -comp_f_H2O/m_H2O*mean(janaf('H2O',T_ec)) -t_0*(S_CH4+c_pc*log(t_0/273.15))...
    -t_0*comp_f_O2/m_O2*(S_O2+mean(janaf('O2',T_ec))*log(t_0/273.15)-(8.314/M_O2)*log(20.64/100))...
    +t_0*comp_f_CO2/m_CO2*(S_CO2+mean(janaf('CO2',T_ec))*log(t_0/273.15)-(8.314/M_CO2)*log(0.04/100))...
    +t_0*comp_f_H2O/m_H2O*(S_H2O+mean(janaf('H2O',T_ec))*log(t_0/273.15)); % cf. eq 1.52


c_pg = comp_f_CO2*janaf('CO2',400) + comp_f_H2O*janaf('H2O',400)...
     + comp_f_O2*janaf('O2',400) + comp_f_N2*janaf('N2',400);
 
%%%%%%%%%%%%%%%%%%
% Energetic loss %
%%%%%%%%%%%%%%%%%%

perte_mecen = k_mec*(m_f*(h_3-h_4)+m_a*(h_2-h_1)); % cf. eq 3.29
perte_echen = m_f*(h_4-h_1); 
W_m = m_fa*(h_3-h_4)-(h_2-h_1); % cf. eq 3.10

%%%%%%%%%%%%%%%%%%
% Exergetic loss %
%%%%%%%%%%%%%%%%%%

perte_mecex = perte_mecen;
perte_rotex = m_f*(e_3-e_4)-m_a*(e_2-e_1)-m_a*W_m; % cf. eq 3.34
perte_combex = m_c*e_c-(m_f*e_3-m_a*e_2); % cf. eq 3.36
perte_echex  = m_f*(e_4-e_1);

%%%%%%%%%%%%%%%%
% Efficiencies %
%%%%%%%%%%%%%%%%

eta_cyclen = W_m/Q_comb; % cf. eq 3.28
eta_mec = 1-(perte_mecen/(P_e+perte_mecen)); % cf. eq 3.30
eta_toten = eta_mec*eta_cyclen; % cf. eq 3.31
eta_cyclex = m_a*W_m/(m_f*e_3-m_a*e_2); % cf. eq 3.32
eta_rotex = m_a*W_m/(m_f*(e_3-e_4)-m_a*(e_2-e_1)); % cf. eq 3.34
eta_combex = (m_f*e_3-m_a*e_2)/(m_c*e_c); % cf. eq 3.36

%%%%%%%%%%%%%%%%%%%%%%%%
% Allocation solutions %
%%%%%%%%%%%%%%%%%%%%%%%%

% ETA vector
ETA(1) = eta_cyclen; ETA(2) = eta_toten; ETA(3) = eta_cyclex; 
ETA(4) = eta_rotex; ETA(5) = eta_combex; 

% DATEN vector
DATEN(1) = perte_mecen; DATEN(2) = perte_echen;

% DATEX vector
DATEX(1) = perte_mecex; DATEX(2) = perte_rotex; DATEX(3) = perte_combex; 
DATEX(4) = perte_echex; 

% DAT Vector
DAT(1,1) = T_1; DAT(1,2) = T_2; DAT(1,3) = T_3; DAT(1,4) = T_4;
DAT(2,1) = p_1; DAT(2,2) = p_2; DAT(2,3) = p_3; DAT(2,4) = p_4; 
DAT(3,1) = h_1; DAT(3,2) = h_2; DAT(3,3) = h_3; DAT(3,4) = h_4; 
DAT(4,1) = s_1; DAT(4,2) = s_2; DAT(4,3) = s_3; DAT(4,4) = s_4; 
DAT(5,1) = e_1; DAT(5,2) = e_2; DAT(5,3) = e_3; DAT(5,4) = e_4; 

% MASSFLOW vector
MASSFLOW(1) = m_a; MASSFLOW(2) = m_c; MASSFLOW(3) = m_f; 

% COMBUSTION structure
COMBUSTION.LHV = l_hv; COMBUSTION.e_c = e_c; COMBUSTION.lambda = lambda;
COMBUSTION.Cp_g = c_pg; COMBUSTION.fum(1) = m_O2; COMBUSTION.fum(2) = m_N2;
COMBUSTION.fum(3) = m_CO2; COMBUSTION.fum(4) = m_H2O;

% Display
ETA
DATEN
DATEX
DAT
MASSFLOW
COMBUSTION
end
