function [ETA MASSFLOW FIG] = CCGT3P(P_eg,options,display)
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
           P_eg=225e3;%100MW
       end
   end
end

if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 15; 
end
if isfield(options,'T_STmax')
    T_STmax = options.T_STmax;
else
    T_STmax = 565; 
end
if isfield(options,'eta_mec')
    eta_mec = options.eta_mec;
else
    eta_mec = 0.99; 
end
if isfield(options,'pdrum')
    p_drum = options.pdrum;
else
    p_drum = 4;
end
if isfield(options,'pmid')
    p_mid = options.pmid;
else
    p_mid = 28; 
end

if isfield(options,'x7')
    x_7 = options.x7;
else
    x_7 = 0.95;
end
if isfield(options,'eta_SiC')
    eta_SiC = options.eta_SiC;
else
    eta_SiC = 0.8;  
end
if isfield(options,'eta_SiT')
    eta_SiT = options.eta_SiT;
else
    eta_SiT = [0.9 0.9];  
end

%% DATA
T_cond = 35;
Tpinchf = 10;

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

% h_6=0;
% eta_SiT(2)=0;
% while round(h_6) ~= 3104
%     h_6 = h_5-eta_SiT(2)*(h_5-h_6s);
%     eta_SiT(2) = eta_SiT(2)+0.001;
% end
h_6 = h_5-eta_SiT(2)*(h_5-h_6s);
s_6 = XSteam('s_ph',p_6,h_6);
T_6 = XSteam('t_hs',h_6,s_6);
e_6 = (h_6-h_0) - (T_0+273.15)*(s_6-s_0);

%% State 7 (x_7 connu et T_7 = T_1)
T_7 = T_cond; % Condensation isobare
p_7  = XSteam('psat_T',T_7);
h_7 = XSteam('h_Tx',T_7,x_7);
s_7 = XSteam('s_pH',p_7,h_7);
e_7  = (h_7-h_0) - (T_0+273.15)*(s_7-s_0);

%% State 1 ( T_1 = T_7 et p_1 = p_7)
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

% h_2=0;
% eta_SiC=0;
% while round(h_2,1) ~= 152
%     eta_SiC = eta_SiC+0.001;
%     h_2 = h_1+(h_2s-h_1)/eta_SiC;
% end
h_2 = h_1+(h_2s-h_1)/eta_SiC;
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
T_8 = T_9p;
h_8 = XSteam('h_pT',p_8,T_8);
s_8 = XSteam('s_pT',p_8,T_8);
e_8 = (h_8-h_0) - (T_0+273.15)*(s_8-s_0);

%% State 9 
p_9 = p_mid;
T_9 = T_6;
h_9 = XSteam('h_pT',p_9,T_9);
s_9 = XSteam('s_pT',p_9,T_9);
e_9 = (h_9-h_0) - (T_0+273.15)*(s_9-s_0);

%% State 10p (vaporisation isobare)
T_10p = T_9;
x_10p = 0;
p_10p = XSteam('psat_T',T_10p);
h_10p = XSteam('h_px',p_10p,x_10p);
s_10p = XSteam('s_ph',p_10p,h_10p);
e_10p = (h_10p-h_0) - (T_0+273.15)*(s_10p-s_0);

%% State 10pp (vaporisation isobare)
p_10pp = p_10p;
x_10pp = 1;
T_10pp = XSteam('Tsat_p',p_10pp);
h_10pp = XSteam('h_px',p_10pp,x_10pp);
s_10pp = XSteam('s_ph',p_10pp,h_10pp);
e_10pp = (h_10pp-h_0) - (T_0+273.15)*(s_10pp-s_0);

%% State after pump LP-IP (for efficiency computing)
p_8ppp = p_9p;
p_8ppps = p_8ppp;
s_8ppps = s_8p;
h_8ppps = XSteam('h_ps',p_8ppps,s_8ppps);
h_8ppp = (h_8ppps-h_8p)/eta_SiC+h_8p;

%% State after pump IP-HP
p_9ppp = p_10p;
p_9ppps = p_9ppp;
s_9ppps = s_9p;
h_9ppps = XSteam('h_ps',p_9ppps,s_9ppps);
h_9ppp = (h_9ppps-h_9p)/eta_SiC+h_9p;

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

% h_4=0;
% eta_SiT(1)=0;
% while round(h_4) ~= 3162
%     eta_SiT(1) = eta_SiT(1)+0.001;
%     h_4 = h_3-eta_SiT(1)*(h_3-h_4s);
% end
h_4 = h_3-eta_SiT(1)*(h_3-h_4s);
s_4 = XSteam('s_ph',p_4,h_4);
T_4 = XSteam('T_ph',p_4,h_4);
e_4 = (h_4-h_0) - (T_0+273.15)*(s_4-s_0);

%% Gas Turbine
[ETA_GT,DATEN_GT,DATEX_GT,DAT_GT,MASSFLOW_GT,COMBUSTION_GT,Cp_g,FIG] = GT_CCGT(P_eg,options,display);

p_4g = DAT_GT(2,4);
T_4g = DAT_GT(1,4)+273.15; 
h_4f = DAT_GT(3,4);
s_4f = DAT_GT(4,4);
e_4f = DAT_GT(5,4);

p_1f = DAT_GT(2,1);
T_1f = DAT_GT(1,1)+273.15;
h_1f = DAT_GT(3,1);
s_1f = DAT_GT(4,1);
e_1f = DAT_GT(5,1);

LHV = COMBUSTION_GT.LHV;
e_c = COMBUSTION_GT.e_c; 
compf = COMBUSTION_GT.compf;
R_f = COMBUSTION_GT.R_f;

m_a = MASSFLOW_GT(1);
m_c = MASSFLOW_GT(2); 
m_f = MASSFLOW_GT(3); 

%% Mass flow
T_HPf = T_10p+Tpinchf; 
T_MPf = T_9p+Tpinchf;
T_LPf = T_8p+Tpinchf;

h_HPf = cp_fumes(300,(T_HPf+273.15),compf)*((T_HPf+273.15)-300);
h_MPf = cp_fumes(300,(T_MPf+273.15),compf)*((T_MPf+273.15)-300);
h_LPf = cp_fumes(300,(T_LPf+273.15),compf)*((T_LPf+273.15)-300); 

A = [h_8-h_8p h_9p-h_8p h_9p-h_8p;h_6-h_8 h_9-h_9p h_10p-h_9p;0 h_5-h_9 h_3-h_10p+h_5-h_9];
B = [m_f*(h_MPf-h_LPf);m_f*(h_HPf-h_MPf);m_f*(h_4f-h_HPf)];

m_v = linsolve(A,B);
m_vLP = m_v(1);
m_vMP = m_v(2);
m_vHP = m_v(3);

MASSFLOW = [m_vHP m_vHP+m_vMP m_vHP+m_vMP+m_vLP m_a m_c];

%% State 5g
p_5f = 1;
f = @(h_5f) m_f*(h_LPf-h_5f)-(m_vLP+m_vMP+m_vHP)*(h_8p-h_2); 
h_5f = fsolve(f,1);

T_5f = 1000; T_it = 0;
while abs((T_5f-T_it))>1e-5
    T_it = T_5f;
    T_5f = (h_5f/cp_fumes(300,T_5f+273.15,COMBUSTION_GT.compf))+300;
end
s_5f = cp_fumes(300,T_5f,COMBUSTION_GT.compf)*log(T_5f/300)-R_f*log(p_5f/1);
e_5f = e_1f+(h_5f-h_1f)-273.15*(s_5f-s_1f);

%%%%%%%%%%%%%%%%
% Efficiencies %
%%%%%%%%%%%%%%%%
P_op_LMP = MASSFLOW(2)*(h_8ppp-h_8p);
P_op_MHP = MASSFLOW(1)*(h_9ppp-h_9p);
P_op_pump = MASSFLOW(3)*(h_2-h_1);
P_op = P_op_LMP+P_op_MHP+P_op_pump;

P_mov_HP = MASSFLOW(1)*(h_3-h_4);
P_mov_IP = MASSFLOW(2)*(h_5-h_6);
P_mov_LP = MASSFLOW(3)*(h_6-h_7);
P_mov = P_mov_HP+P_mov_IP+P_mov_LP;

P_mcy_ST = (P_mov-P_op);
P_elec_GT = P_eg;
P_elec_ST = P_mcy_ST*eta_mec;

QI = m_vHP*(h_3-h_2)+m_vHP*(h_5-h_4)+m_vMP*(h_5-h_2)+m_vLP*(h_6-h_2);
QI_ex = m_vHP*(e_3-e_2)+m_vHP*(e_5-e_4)+m_vMP*(e_5-e_2)+m_vLP*(e_6-e_2);

ETA(1) = P_mcy_ST/QI;
ETA(2) = ETA_GT(1);
ETA(3) = (P_elec_ST*10^3+P_elec_GT*10^3)/(m_c*LHV);
ETA(4) = P_mcy_ST/QI_ex;
ETA(5) = ETA_GT(3);
ETA(6) = (P_elec_ST+P_elec_GT)/(m_c*e_c);
ETA(7) = QI/(m_c*LHV/1000);
ETA(8) = QI_ex/(m_c*e_c);
ETA(9) = ETA_GT(6);
ETA(10) = (e_4f-e_5f)/e_4f;
ETA(11) = QI_ex/(m_f*(e_4f-e_5f));

%% FIGURES
DAT_plot = [T_2, T_3, T_4, T_5, T_6, T_7, T_8, T_8p, T_8pp, T_9, T_9p, T_9pp, T_10p, T_10pp;...
            h_2, h_3, h_4, h_5, h_6, h_7, h_8, h_8p, h_8pp, h_9, h_9p, h_9pp, h_10p, h_10pp;
            s_2, s_3, s_4, s_5, s_6, s_7, s_8, s_8p, s_8pp, s_9, s_9p, s_9pp, s_10p, s_10pp];
FIG = [];
if display    
    Tplot = linspace(0.01,373.9458,1000);
    for i=1:1000
        sL_plot(i) = XSteam('sL_T',Tplot(i));
        sV_plot(i) = XSteam('sV_T',Tplot(i));
        hL_plot(i) = XSteam('hL_T',Tplot(i));
        hV_plot(i) = XSteam('hV_T',Tplot(i));
    end
    
    FIG(1) = figure;
    plot(sV_plot,Tplot,'b','LineWidth',2);
    hold on
    plot(sL_plot,Tplot,'b','LineWidth',2);
    hold on;
    plotData = [1,p_3,p_4,s_3,h_3,h_4,eta_SiT(2);
                1,p_5,p_7,s_5,h_5,h_7,eta_SiT(2);
                2,p_2,p_3,s_2,h_2,h_3,eta_SiT(1);
                2,p_2,p_5,s_2,h_2,h_5,eta_SiT(1);
                2,p_2,p_6,s_2,h_2,h_6,eta_SiT(1)];
    for i=1:length(plotData(:,1))
        hold on;
        [vect vech vecs] = CCGTPlot(plotData(i,:));
        plot(vecs,vect,'r','Linewidth',2);
    end
    hold on
    plot([s_2 s_7],[T_2 T_7],'r','Linewidth',2);
    hold on
    scatter(DAT_plot(3,:),DAT_plot(1,:),'k','filled','Linewidth',2);
    labels = {'1=2','3','4','5','6','7','8','8''','8''''','9','9''','9''''','10''','10'''''};
    text(DAT_plot(3,:),DAT_plot(1,:),labels,'FontSize',12,'VerticalAlignment','bottom','HorizontalAlignment','right')
    xlabel('Entropy [kJ/kgK]');
    ylabel('Temperature [°C]');
    title('T-S Diagram of CCGT');
    
    FIG(2) = figure;
    plot(sV_plot,hV_plot,'b','LineWidth',2);
    hold on
    plot(sL_plot,hL_plot,'b','LineWidth',2);
    hold on;
    for i=1:length(plotData(:,1))
        hold on;
        [vect vech vecs] = CCGTPlot(plotData(i,:));
        plot(vecs,vech,'r','Linewidth',2);
        scatter(vecs(end),vech(end),'k','filled','Linewidth',2);
        scatter(vecs(1),vech(1),'k','filled','Linewidth',2);
    end
    plot([s_2 s_7],[h_2 h_7],'r','Linewidth',2);
    scatter([s_2 s_7],[h_2 h_7],'k','filled','Linewidth',2);
    xlabel('Entropy [kJ/kgK]');
    ylabel('Enthalpy [kJ/kg]');
    title('H-S Diagram of CCGT'); 
end

end