function [ETA MASSFLOW FIG] = CCGT(P_eg,options,display)
% CCGT is a Combine cycle Gas Turbine with 2 pressure level
% CCGT(P_e,options,display) compute the thermodynamics states for a CCGT
% with 3 pressure level (cfr p166 english reference book) including
% combustion, exchanger and cycles. This is done based on several inputs
% (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input 
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated) Refer to Fig 4.19 from reference book (in english)
% P_EG = electrical power output target for gas turbine [kW]
% OPTIONS is a structure containing :
%   -options.T0       [°C] : Reference temperature
%   -options.T_ext    [°C] : External temperature
%   -options.T_STmax  [°C] : maximum temperature on ST cycle
%   -options.eta_mec  [-] : mecanic efficiency of shafts bearings
%   -options.pdrum   [bar]: Drum pressure
%   -options.pmid    [bar]: Intermediary pressure level
%   -options.x7       [-] : Vapor ratio [gaseous/liquid] (titre)
%   -option.eta_SiC   [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT   [-] : Isotrenpic efficiency for compression
%   -options.GT    [struct] : options for Gas turbine (see GT function) 
% DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then 
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1)  : eta_STcyclen, cycle energy efficiency
%   -eta(2)  : eta_GTcyclen, cycle energy efficiency
%   -eta(3)  : eta_toten, overall energy efficiency
%   -eta(4)  : eta_STcyclex, cycle exegy efficiency
%   -eta(5)  : eta_GTcyclex, cycle exegy efficiency
%   -eta(6)  : eta_totex, overall exergie efficiency
%   -eta(7)  : eta_gen, Steam generator energy efficiency
%   -eta(8)  : eta_gex, Steam generator exergy efficiency
%   -eta(9)  : eta_combex, Combustion exergy efficiency
%   -eta(10) : eta_chemex, Chimney exergy efficiency (losses)
%   -eta(11) : eta_transex, Heat exchanger overall exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% MASSFLOW is a vector containing : 
%   -massflow(1) [kg/s]: water massflow at high pressure turbine inlet
%   -massflow(2) [kg/s]: water massflow at medium pressure turbine inlet
%   -massflow(3) [kg/s]: water massflow at low pressure turbine inlet
%   -massflow(4) [kg/s]: air massflow at gas turbine inlet 
%   -massflow(5) [kg/s]: combustible massflow
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

%% Your work
if nargin<3
    display=1;
   if nargin<2
       options=struct();
       if nargin<1
           P_eg=100e3;%100MW
       end
   end
end

if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 15; % [°C]
end

if isfield(options,'T_STmax')
    T_STmax = options.T_STmax;
else
    T_STmax = 565; % [°C]
end

if isfield(options,'eta_mec')
    eta_mec = options.eta_mec;
else
    eta_mec = 0.99; % []
end

if isfield(options,'pdrum')
    p_drum = options.pdrum;
else
    p_drum = 4; % [bar]
end

if isfield(options,'pmid')
    p_mid = options.pmid;
else
    p_mid = 28; % [°C]
end

if isfield(options,'x7')
    x_7 = options.x7;
else
    x_7 = 0.95; % [°C]
end

if isfield(options,'eta_SiC')
    eta_SiC = options.eta_SiC;
else
    eta_SiC = 0.8;  
end
if isfield(options,'eta_SiT')
    eta_SiT = options.eta_SiT;
else
    eta_SiT = [0.8 0.8];  
end

%% DATA
T_cond = 24;
TpinchCond = 12;

%% State 0
p_0=XSteam('psat_T',T_0);
x_0=0;
h_0=XSteam('H_px',p_0,x_0);
s_0=XSteam('s_ph',p_0,h_0);
e_0=0;

%% State 5 (T_5 = T_3, p_5 = p_mid)
T_5 = T_STmax;
p_5 = p_mid;
h_5 = XSteam('h_pT',p_5,T_5);
s_5 = XSteam('s_pT',p_5,T_5);
e_5 = (h_5-h_0) - (T_0+273.15)*(s_5-s_0);

%% State 6 (Détente isentropique : s_6s = s_5, p_6s = p_6)
s_6s = s_5;
p_6 = p_drum;
p_6s = p_6;
h_6s = XSteam('h_ps',p_6s,s_6s);

h_6=0;
eta_SiT(2)=0;
while round(h_6) ~= 3104
    h_6 = h_5-eta_SiT(2)*(h_5-h_6s);
    eta_SiT(2) = eta_SiT(2)+0.001;
end

s_6 = XSteam('s_ph',p_6,h_6);
T_6 = XSteam('t_hs',h_6,s_6);
e_6 = (h_6-h_0) - (T_0+273.15)*(s_6-s_0);

%% State 7 (x_7 connu et T_7 = T_1)
T_7 = T_cond+TpinchCond;
p_7  = XSteam('psat_T',T_7);
h_7 = XSteam('h_Tx',T_7,x_7);
s_7 = XSteam('s_pH',p_7,h_7);
e_7  = (h_7-h_0) - (T_0+273.15)*(s_7-s_0);

%% State 1 (x_1 = 1, T_1 = T_7 et p_1 = p_7)
T_1 = T_7;
p_1 = p_7;
h_1 = XSteam('hL_T',T_1);
s_1 = XSteam('sL_T',T_1);
e_1  = (h_1-h_0) - (T_0+273.15)*(s_1-s_0);

%% State 2 (p_2 = p_drum)
p_2 = p_drum;
p_2s = p_2;
s_2s = s_1;
h_2s = XSteam('h_ps',p_2s,s_2s);

h_2=0;
eta_SiC=0.001;
while round(h_2,1) ~= 152
    h_2 = h_1+(h_2s-h_1)/eta_SiC;
    eta_SiC = eta_SiC+0.001;
end

T_2 = XSteam('t_ph',p_2,h_2);
s_2 = XSteam('s_ph',p_2,h_2);
e_2 = (h_2-h_0) - (T_0+273.15)*(s_2-s_0);

%% State 8p (vaporisation isobare)
p_8p = p_drum;
x_8p = 0;
h_8p = XSteam('h_px',p_8p,x_8p);
s_8p = XSteam('s_ph',p_8p,h_8p);
T_8p = XSteam('T_ph',p_8p,h_8p);
e_8p = (h_8p-h_0) - (T_0+273.15)*(s_8p-s_0);

%% State 8pp (vaporisation isobare)
p_8pp = p_drum;
x_8pp = 1;
h_8pp = XSteam('h_px',p_8pp,x_8pp);
s_8pp = XSteam('s_ph',p_8pp,h_8pp);
T_8pp = XSteam('T_ph',p_8pp,h_8pp);
e_8pp = (h_8pp-h_0) - (T_0+273.15)*(s_8pp-s_0);

%% State 9p (vaporisation isobare)
p_9p = p_mid;
x_9p = 0;
h_9p = XSteam('h_px',p_9p,x_9p);
s_9p = XSteam('s_ph',p_9p,h_9p);
T_9p = XSteam('T_ph',p_9p,h_9p);
e_9p = (h_9p-h_0) - (T_0+273.15)*(s_9p-s_0);

%% State 9pp (vaporisation isobare)
p_9pp = p_mid;
x_9pp = 1;
h_9pp = XSteam('h_px',p_9pp,x_9pp);
s_9pp = XSteam('s_ph',p_9pp,h_9pp);
T_9pp = XSteam('T_ph',p_9pp,h_9pp);
e_9pp = (h_9pp-h_0) - (T_0+273.15)*(s_9pp-s_0);

%% State 8
p_8 = p_drum;
T_8 = T_9p; % Voir tableau 4.6, je n'arrive pas encore à le justifier 
h_8 = XSteam('h_pT',p_8,T_8);
s_8 = XSteam('s_pT',p_8,T_8);
e_8 = (h_8-h_0) - (T_0+273.15)*(s_8-s_0);

%% State 9 
p_9 = p_mid;
T_9 = T_6; % Voir tableau 4.6, je n'arrive pas encore à le justifier 
h_9 = XSteam('h_pT',p_9,T_9);
s_9 = XSteam('s_pT',p_9,T_9);
e_9 = (h_9-h_0) - (T_0+273.15)*(s_9-s_0);

%% State 10p (vaporisation isobare)
T_10p = T_6; % Voir tableau 4.6, je n'arrive pas encore à le justifier 
x_10p = 0;
p_10p = XSteam('psat_T',T_10p);
h_10p = XSteam('h_px',p_10p,x_10p);
s_10p = XSteam('s_ph',p_10p,h_10p);
e_10p = (h_10p-h_0) - (T_0+273.15)*(s_10p-s_0);

%% State 10pp (vaporisation isobare)
T_10pp = T_6; % Voir tableau 4.6, je n'arrive pas encore à le justifier 
x_10pp = 1;
p_10pp = XSteam('psat_T',T_10pp);
h_10pp = XSteam('h_px',p_10pp,x_10pp);
s_10pp = XSteam('s_ph',p_10pp,h_10pp);
e_10pp = (h_10pp-h_0) - (T_0+273.15)*(s_10pp-s_0);

% On arrive à p_10p = p_10pp = 101 bars au lieu de 110... A corriger

%% State 3 
p_3 = p_10p;
T_3 = T_STmax;
h_3 = XSteam('h_pT',p_3,T_3);
s_3 = XSteam('s_pT',p_3,T_3);
e_3 = (h_3-h_0) - (T_0+273.15)*(s_3-s_0);

%% State 4
s_4s = s_3;
p_4 = p_mid;
p_4s = p_4;
h_4s = XSteam('h_ps',p_4s,s_4s);

h_4=0;
eta_SiT(1)=0;
while round(h_4) ~= 3162
    h_4 = h_3-eta_SiT(1)*(h_3-h_4s);
    eta_SiT(1) = eta_SiT(1)+0.001;
end

s_4 = XSteam('s_ph',p_4,h_4);
T_4 = XSteam('T_ph',p_4,h_4);
e_4 = (h_4-h_0) - (T_0+273.15)*(s_4-s_0);



% Itérations pas correcte si on change les valeurs de options. Mais c'est
% juste pour avoir une estimation des rendements isentropiques optimaux
end