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
%   -combustion.e_c    : the combuistible exergie         [kJ/kg]
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
           P_e = 230e3;%100MW
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

% Cycle State 1 - Before Compression
T_1 = T_ext;
t_1 = T_ext + 273.15;
p_1 = 1;
c_pa = 1.006;

h_1 = c_pa*T_1;
s_1 = 0.058;
e_1 = 0;

% Cycle State 2 - After Compression
p_2 = r*p_1;
R_sa = 287.058/1000;

step = 1;
t_2 = 700;
t_2_iter = 0;
T1 = 300.*ones(301,1)';

while abs(t_2-t_2_iter)>0.01
    t_2_iter = t_2;
    T_c1 = [T1 301:step:t_1];
    T_c2 = [T1 301:step:t_2];
    c_pa_moy1 = (mean(janaf('N2',T_c1))*0.79 + mean(janaf('O2',T_c1))*0.21);
    c_pa_moy2 = (mean(janaf('N2',T_c2))*0.79 + mean(janaf('O2',T_c2)*0.21));
    t_2 = t_1*r^((1/eta_PiC)*(R_sa*((t_2-273.15)- T_1)/(c_pa_moy2*(t_2-273.15)-c_pa_moy1*T_1)));
end


T_2 = t_2-273.15;

c_pa_moy12 = (c_pa_moy2*T_2-c_pa_moy1*T_1)/(T_2-T_1);
h_2 = h_1 + c_pa_moy12*(T_2-T_1);
s_2 = s_1 + (1-eta_PiC)*c_pa_moy12*log(t_2/t_1);
e_2 = (h_2-h_1)-(T_0+273.15)*(s_2-s_1);

% Cycle State 3 - After Combustion
t_3 = T_3 + 273.15;
p_3 = k_cc*p_2;

M_O2 = 31.998000000000000/1000;
M_CO2 = 44.008000000000000/1000;
M_H2 = 2.015940000000000/1000;
M_H2O = 18.014940000000000/1000;
M_N2 = 28.014000000000000/1000;
M_CO = 28.0090000000000/1000;
M_CH4 = 16.0430000000/1000;

T1 = 300.*ones(301,1)';
T_react = [T1 301:step:t_2];
T_prod = [T1 301:step:t_3];
x = 0;
y = 4;
l_hv = 51.5*10^3;

T_c1 = [T1 301:step:t_1];
T_c2 = [T1 301:step:t_2];
T_c3 = [T1 301:step:t_3];

c_pa_mat = ones(2,6);

c_pa_mat(1,1) = (mean(janaf('N2',T_c1)));
c_pa_mat(2,1) = (mean(janaf('N2',T_c2)));
c_pa_mat(1,2) = (mean(janaf('O2',T_c1)));
c_pa_mat(2,2) = (mean(janaf('O2',T_c2)));
c_pa_mat(1,3) = (mean(janaf('CO2',T_c1)));
c_pa_mat(2,3) = (mean(janaf('CO2',T_c3)));
c_pa_mat(1,4) = (mean(janaf('H2O',T_c1)));
c_pa_mat(2,4) = (mean(janaf('H2O',T_c3)));
c_pa_mat(1,5) = (mean(janaf('O2',T_c1)));
c_pa_mat(2,5) = (mean(janaf('O2',T_c3)));
c_pa_mat(1,6) = (mean(janaf('N2',T_c1)));
c_pa_mat(2,6) = (mean(janaf('N2',T_c3)));

f = @(lambda) -l_hv - t_2*(lambda*(1+(y-2*x)/4)*(c_pa_mat(1,2)*T_2-c_pa_mat(2,2)*T_1)/(T_2-T_1)...
-   lambda*(1+(y-2*x)/4)*3.76*(c_pa_mat(1,1)*T_2-c_pa_mat(2,1)*T_1)/(T_2-T_1)) + t_3*((c_pa_mat(1,3)*T_3-c_pa_mat(2,3)*T_1)/(T_3-T_1)...
+ y/2*(c_pa_mat(1,4)*T_3-c_pa_mat(2,4)*T_1)/(T_3-T_1) + (1+(y-2*x)/4)*(lambda-1)*(c_pa_mat(1,5)*T_3-c_pa_mat(2,5)*T_1)/(T_3-T_1)...
+ (1+(y-2*x)/4)*3.76*lambda*(c_pa_mat(1,6)*T_3-c_pa_mat(2,6)*T_1)/(T_3-T_1));
    
lambda = fsolve(f,1);

M_a = (32+3.76*28)*(1+(y-2*x)/4)/(12+y+16*x);
M_s = 1 + M_a;

m_ac = M_a*lambda;
m_fa = (1+ 1/(m_ac));
Q_comb = l_hv/m_ac;

h_3 = (Q_comb+h_2)*(1/m_fa);
%s_3 = ?
%e_3 = (h_3-h_3)-(T_0+273.15)*(s_3-s_2);

% Cycle State 4 - After decompression
p_4 = p_3/(k_cc*r);

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

step = 1;
t_4 = 700;
t_4_iter = 0;
T1 = 300.*ones(301,1)';

while abs(t_4-t_4_iter)>0.01
    t_4_iter = t_4;
    T_c3 = [T1 301:step:t_3];
    T_c4 = [T1 301:step:t_4];
    c_pa_moy3 = mean(janaf('CO2',T_c3)*comp_f_CO2 + janaf('H2O',T_c3)*comp_f_H2O...
              + janaf('O2',T_c3)*comp_f_O2 + janaf('N2',T_c3)*comp_f_N2);
    c_pa_moy4 = mean(janaf('CO2',T_c4)*comp_f_CO2 + janaf('H2O',T_c4)*comp_f_H2O...
              + janaf('O2',T_c4)*comp_f_O2 + janaf('N2',T_c4)*comp_f_N2);
    t_4 = t_3*r^(-eta_PiT*(R_f*(T_3-(t_4-273.15))/(c_pa_moy3*T_3-c_pa_moy4*(t_4-273.15))));
end

T_4 = t_4-273.15;

c_pa_moy34 = (c_pa_moy3*T_3-c_pa_moy4*T_4)/(T_3-T_4);
h_4 = h_3 + c_pa_moy34*(T_4-T_3);
% s_4 = s_3 - (1-eta_PiT)/eta_PiT*c_pa_moy34*log(t_4/t_3);
% e_4 = (h_4-h_3)-(T_0+273.15)*(s_4-s_3);

% Allocation solutions
DAT(1,1) = T_1;
DAT(1,2) = T_2;
DAT(1,3) = T_3;
DAT(1,4) = T_4;
DAT(2,1) = p_1;
DAT(2,2) = p_2;
DAT(2,3) = p_3;
DAT(2,4) = p_4;
DAT(3,1) = h_1;
DAT(3,2) = h_2;
DAT(3,3) = h_3;
DAT(3,4) = h_4
DAT(4,1) = s_1;
DAT(4,2) = s_2;
DAT(4,3) = 0;
DAT(4,4) = 0;

end
