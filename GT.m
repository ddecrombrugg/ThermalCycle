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


%% Parameters Verification

if nargin < 3
   display = 0;
   if nargin < 2
       options = struct();
       if nargin < 1
           P_e = 100e3; % 100[MW]
       end
   end
end

if isfield(options,'T_0') == 0
    options.T_0 = 273.15;
else
    options.T_0 = options.T_0 + 273.15;
end

if isfield(options,'T_ext') == 0
    options.T_ext = 288.15;
else
    options.T_ext = options.T_ext + 273.15;
end

if isfield(options,'r') == 0
    options.r = 18;
end

if isfield(options,'T_3') == 0
    options.T_3 = 1400 + 273.15;
else
    options.T_3 = options.T_3 + 273.15;
end

options.eta_PiC = .9;
options.eta_PiT = .9;
options.eta_kcc = .95;


%% Other parameters

p_ext = 100e3; % [Pa]

R = 8.314472; % The ideal gas's constant [J/mol/K]
R_O2  = R / 31.99800e-3; % [J/(kg*K)]
R_N2  = R / 28.01400e-3; % [J/(kg*K)]
R_CO2 = R / 44.00800e-3; % [J/(kg*K)]
R_H2O = R / 18.01494e-3; % [J/(kg*K)]
R_air = .21*R_O2 + .79*R_N2; % [J/(kg*K)]

T_vect = [ones(1,300)*300,300:5000];
Cp_O2  = janaf('O2',T_vect) *1e3; % [J/(kg*K)]
Cp_N2  = janaf('N2',T_vect) *1e3;
Cp_CO2 = janaf('CO2',T_vect) *1e3;
Cp_H2O = janaf('H2O',T_vect) *1e3;
Cp_air = .21*Cp_O2 + .79*Cp_N2;
T_vect = 0:5000;

gamma  = Cp_air./(Cp_air - R_air); % Laplace coefficient of a diatomic atom [-]

%% Bases for the calculations of the combustion
% In the combustion chamber where the compressed air is coming at (T_2,
% p_2,h_2,s_2,e_2) at a certain flow m_air, we inject CH_4 (T_CH4,m_CH4)
% so that the following transformation can happen to create energy :
% CH_4 + w*(O_2 + 3.76*N_2) -> CO_2 + a*O_2 + 2*H_2O + 3.76*w*N_2

%sym a w
%sol = solve([w - a - 2; options.r*w - a] == [0;0],[a,w]);
%a = sol.a; w = sol.w;
%lambda = w/2;


%% Calculations of all states such that
% 1 -> 2 : polytropic compression
% 2 -> 3 : isobaric warming
% 3 -> 4 : polytropic relaxation
% 4 -> 1 : isobaric cooling


% First state
p_1 = p_ext;
T_1 = options.T_ext;
h_1 = Cp_air(int16(T_1)) * (T_1 - options.T_0);
s_1 = Cp_air(int16(T_1)) * log(T_1/options.T_0);

% Second state
p_2 = p_1 * options.r;

err = 10; iter = 0;
try_1 = 500; try_2 = 1000;
while abs(err) > 1e-3 && iter < 100
   exp_1 = R_air / mean(Cp_air(int16(T_1):int16(try_1)));
   T_21  = T_1 * options.r ^(exp_1/options.eta_PiC);
   err_1 = try_1 - T_21;
   
   exp_2 = R_air / mean(Cp_air(int16(T_1):int16(try_2)));
   T_22  = T_1 * options.r ^(exp_2/options.eta_PiC);
   err_2 = try_2 - T_22;
   
   if err_1 < 0 && err_2 < 0
      if err_1 > err_2
         try_2 = try_1 *2;
         T_2 = T_21;
         err = err_1;
      else
         try_1 = try_2 *2;
         T_2 = T_22;
         err = err_2;
      end
   elseif err_1 > 0 && err_2 > 0
      if err_1 < err_2
         try_2 = try_1 /2;
         T_2 = T_21;
         err = err_1;
      else
         try_1 = try_2 /2;
         T_2 = T_22;
         err = err_2;
      end
   elseif abs(err_1) < abs(err_2)
      try_2 = (try_1 + try_2) /2;
      T_2 = T_21;
      err = err_1;
   else
      try_1 = (try_1 + try_2) /2;
      T_2 = T_22;
      err = err_2;
   end
   iter = iter +1;
end

h_2 = Cp_air(int16(T_2)) * (T_2 - options.T_0);
s_2 = Cp_air(int16(T_2)) * log(T_2/options.T_0) - R_air * log(p_2/p_ext);

% Third state
T_3 = options.T_3;
p_3 = p_2 * options.k_cc;
h_3 = Cp_air(int16(T_3)) * (T_3 - options.T_0);
s_3 = Cp_air(int16(T_3)) * log(T_3/options.T_0) - R_air * log(p_3/p_ext);

% Fourth state
p_4 = p_ext;

err = 10; iter = 0;
try_1 = 500; try_2 = 1000;
while abs(err) > 1e-3 && iter < 100
   exp_1 = R_air / mean(Cp_air(int16(try_1):int16(T_3)));
   T_41  = T_3 * (1/options.k_cc/options.r) ^(exp_1*options.eta_PiT);
   err_1 = try_1 - T_41;

   exp_2 = R_air / mean(Cp_air(int16(try_2):int16(T_3)));
   T_42  = T_3 * (1/options.k_cc/options.r) ^(exp_2*options.eta_PiT);
   err_2 = try_2 - T_42;

   if err_1 < 0 && err_2 < 0
      if err_1 > err_2
         try_2 = try_1 *2;
         T_4 = T_41;
         err = err_1;
      else
         try_1 = try_2 *2;
         T_4 = T_42;
         err = err_2;
      end
   elseif err_1 > 0 && err_2 > 0
      if err_1 < err_2
         try_2 = try_1 /2;
         T_4 = T_41;
         err = err_1;
      else
         try_1 = try_2 /2;
         T_4 = T_42;
         err = err_2;
      end
   elseif abs(err_1) < abs(err_2)
      try_2 = (try_1 + try_2) /2;
      T_4 = T_41;
      err = err_1;
   else
      try_1 = (try_1 + try_2) /2;
      T_4 = T_42;
      err = err_2;
   end
   iter = iter +1;
end

h_4 = Cp_air(int16(T_4)) * (T_4 - options.T_0);
s_4 = Cp_air(int16(T_4)) * log(T_4/options.T_0);


DAT = struct();

DAT.T_1 = T_1 -273.15;
DAT.T_2 = T_2 -273.15;
DAT.T_3 = T_3 -273.15;
DAT.T_4 = T_4 -273.15;

DAT.p_1 = p_1 *1e-5;
DAT.p_2 = p_2 *1e-5;
DAT.p_3 = p_3 *1e-5;
DAT.p_4 = p_4 *1e-5;

DAT.h_1 = h_1 *1e-3;
DAT.h_2 = h_2 *1e-3;
DAT.h_3 = h_3 *1e-3;
DAT.h_4 = h_4 *1e-3;

DAT.s_1 = s_1 *1e-3;
DAT.s_2 = s_2 *1e-3;
DAT.s_3 = s_3 *1e-3;
DAT.s_4 = s_4 *1e-3;

DAT.e_1 = DAT.h_1 - options.T_0*DAT.s_1;
DAT.e_2 = DAT.h_2 - options.T_0*DAT.s_2;
DAT.e_3 = DAT.h_3 - options.T_0*DAT.s_3;
DAT.e_4 = DAT.h_4 - options.T_0*DAT.s_4;

%% Figures

FIG = [];
end
