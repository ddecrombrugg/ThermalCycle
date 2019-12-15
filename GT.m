function [ETA, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION, Cp_g, FIG] = ...
    GT(P_e,options,display)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS VERIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3
       display = 1;
       if nargin < 2
           options = struct();
           if nargin < 1
               P_e = 230; % 100[MW]
           end
       end
    end

    if isfield(options,'T_0')
        T_0 = options.T_0 +273.15;
    else
        T_0 = 15 +273.15;
    end

    if isfield(options,'T_ext')
        T_ext = options.T_ext +273.15;
    else
        T_ext = 15 +273.15;
    end

    if isfield(options,'r')
        r = options.r;
    else
        r = 18;
    end

    if isfield(options,'T_3')
        T_3 = options.T_3 +273.15;
    else
        T_3 = 1400 +273.15;
    end

    if isfield(options,'eta_PiC')
        eta_PiC = options.eta_PiC;
    else
        eta_PiC = .9;
    end

    if isfield(options,'eta_PiT')
        eta_PiT = options.eta_PiT;
    else
        eta_PiT = .9;
    end

    if isfield(options,'k_cc')
        k_cc = options.k_cc;
    else
        k_cc = .95;
    end

    if isfield(options,'k_mec')
        k_mec = options.k_mec;
    else
        k_mec = .015;
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OTHER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M_O2  = 31.99800e-3; % [kg/mol]
    M_N2  = 28.01400e-3;
    M_CO2 = 44.00800e-3;
    M_H2O = 18.01494e-3;
    M_air = .21*M_O2 + .79*M_N2;

    p_ext = 100e3; % [Pa]

    R = 8.314472e-3; % The ideal gas's constant [kJ/mol/K]
    R_O2  = R / M_O2; % [kJ/(kg*K)]
    R_N2  = R / M_N2; % [kJ/(kg*K)]
    R_CO2 = R / M_CO2; % [kJ/(kg*K)]
    R_H2O = R / M_H2O; % [kJ/(kg*K)]
    R_air = 287.058e-3; % [kJ/(kg*K)]

    T_vect = [ones(1,300)*300,301:5000];
    Cp_O2  = janaf('O2',T_vect); % [J/(kg*K)]
    Cp_N2  = janaf('N2',T_vect);
    Cp_CO2 = janaf('CO2',T_vect);
    Cp_H2O = janaf('H2O',T_vect);

        function [X] = Cp_air(T)
            T(T<300) = 300; T(T>5000) = 5000;
            X = .21*janaf('O2',T) + .79*janaf('N2',T);
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% REFERENCE  STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIRST  STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p_1 = p_ext;
    T_1 = T_ext;
    h_1 = 1.006 * (T_1 - 273.15);
    s_1 = 0.054;
    e_1 = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECOND STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p_2 = p_1 *r;

    function [T,Cp_moy] = Compression(ratio,T)
        T_it = 0;
        while abs(T-T_it) > 1e-5
            T_it = T;
            Cp_moy = integral(@Cp_air,T_1,T)/(T-T_1);
            T = T_1 * ratio ^(R_air/Cp_moy/eta_PiC); % cf. eq 3.19-3.22
        end
    end

    [T_2,Cp_moy] = Compression(r,1674);
    h_2 = h_1 + Cp_moy*(T_2 - T_1);
    s_2 = s_1 + Cp_moy*log(T_2/T_1) - R_air*log(r);
    e_2 = (h_2-h_1) - T_0*(s_2-s_1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THIRD  STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMBUSTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p_3 = p_2 * k_cc;
    x = 0; y = 4; LHV = 50.15e3; T_CH4 = 298.15;
    coeff1 = [-.703029 108.4773 -42.52157 5.862788 .678565];
    Cp_CH4 = integral(@(T) coeff1(1)+coeff1(2)*(T./1000)...
        +coeff1(3)*((T./1000).^2)+coeff1(4)*((T./1000).^3)...
        +(coeff1(5)./((T./1000).^2)),273.15,T_CH4)/(T_CH4-273.15) *1e-3;

    %Cp_CH4 = 35.639; % Mass heat of methane at standart temperature [J/(mol*K)]
    M_c   = (12.01+1.01*y+16*x)*1e-3; % Molar mass of methane [kg/mol]

    a = @(lam) (lam-1)*(1+y/4-x/2); b = y/2;
    w = @(lam) lam*(1+y/4-x/2); % Stoechiometric coefficients

    t = 273:ceil(T_3);
    fun = @(lam) (T_3 -273.15) * (mean(Cp_CO2(t))*M_CO2 + b*mean(Cp_H2O(t))*M_H2O...
        + a(lam)*mean(Cp_O2(t))*M_O2 + 3.76*w(lam)*mean(Cp_N2(t))*M_N2) ...
        - w(lam)*integral(@Cp_air,273,T_2)*M_air/.21 - 25*Cp_CH4 - LHV*M_c;
    % Bilan d'enthalpie sur la combustion
    opt = optimset('Display','off');
    lambda = fsolve(fun,1,opt); % Exces d'air [mol_air/mol_c]

    a = double(a(lambda)); w = double(w(lambda)); % [mol]

    comp_f_tot = M_CO2 + b*M_H2O + a*M_O2 + 3.76*w*M_N2; % [kg]
    comp_f_CO2 = M_CO2 / comp_f_tot; % [-]
    comp_f_H2O = b*M_H2O / comp_f_tot; % [-]
    comp_f_O2  = a*M_O2 / comp_f_tot; % [-]
    comp_f_N2  = 3.76*w*M_N2 / comp_f_tot; % [-]
    R_f = R / comp_f_tot * (1+b+a+3.76*w); % [kJ/(kg*K)]

    function [X] = Cp_f(T)
        T(T<300) = 300; T(T>5000) = 5000;
        X = comp_f_CO2*janaf('CO2',T) + comp_f_H2O*janaf('H2O',T) + ...
            comp_f_O2*janaf('O2',T) + comp_f_N2*janaf('N2',T); % [J/(kg*K)]
    end

    m_a1 = (1+y/4-x/2) * (M_air/.21)/M_c; % [mol]
    m_ac = lambda * m_a1 ; % [-]
    m_ag = (1 + 1/m_ac)^(-1); % [-]
    Q_comb = LHV / m_ac;

    h_3 = (Q_comb + h_2) * m_ag;
    s_3 = s_2 + integral(@Cp_f,T_2,T_3)*log(T_3/T_2)/(T_3-T_2) - R_f*log(k_cc);
    e_3 = (h_3-h_1) - T_0*(s_3-s_1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOURTH STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p_4 = p_ext;

    function [T,Cp_moy] = Detente(ratio,T)
        T_it = 0;
        while abs(T - T_it) > 1e-5
            T_it = T;
            Cp_moy = integral(@Cp_f,T_3,T)/(T-T_3);
            T = T_3 * ratio^(R_f/Cp_moy*eta_PiT);
        end
    end

    [T_4,Cp_moy] = Detente(1/k_cc/r,500);
    h_4 = h_3 + Cp_moy * (T_4 - T_3);
    s_4 = s_1 + integral(@Cp_f,T_1,T_4)*log(T_4/T_1)/(T_4-T_0);
    e_4 = (h_4-h_1) - T_0*(s_4-s_1);

    DAT = [[T_1, T_2, T_3, T_4]-273.15;... [°C]
            [p_1, p_2, p_3, p_4] *1e-5;...  [bar]
            [h_1, h_2, h_3, h_4];...  [kJ/kg]
            [s_1, s_2, s_3, s_4];...  [kJ/(kg*K)]
            [e_1, e_2, e_3, e_4]];   % [kJ/kg]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EFFICIENCY, MASSFLOWS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    P_e = P_e *1e3; % [kW]

    eta_mec = 1 - k_mec*((h_3-h_4)/m_ag + (h_2-h_1))/((h_3-h_4)/m_ag - (h_2-h_1));
    W_m = (h_3 - h_4)/m_ag - (h_2 - h_1); % [J/kg_air]
    P_m = P_e / eta_mec; % [kW]

    A = [(h_1-h_2), 0, (h_3-h_4);
         -1, lambda*m_a1, 0;
         (1 + lambda*m_a1), 0, -lambda*m_a1];
    b = [P_m;0;0];
    m = A\b;

    m_a = m(1); m_c = m(2); m_g = m(3); % [kg/s]
    MASSFLOW = [m_a, m_c, m_g]; % [kg/s]

    m_CO2f = comp_f_CO2 * m_g; m_H2Of = comp_f_H2O * m_g;
    m_O2f  = comp_f_O2  * m_g; m_N2f  = comp_f_N2  * m_g;

    PCS = 55695;
    e_c = PCS + 15 * (Cp_CH4/M_c + comp_f_O2*Cp_O2(288) - comp_f_CO2*Cp_CO2(288) - comp_f_H2O*Cp_H2O(288)) ...
        - 288.15 * (183.1e-3/M_c + Cp_CH4/M_c*log(288.15/273.15)) ...
        - 288.15 * (202.8e-3/M_O2 + Cp_O2(288)*log(288.15/273.15) - R_O2*log(.2064)) * comp_f_O2 ...
        + 288.15 * (210.4e-3/M_CO2 + Cp_CO2(288)*log(288.15/273.15) - R_CO2*log(.0003)) * comp_f_CO2 ...
        + 288.15 * (69.5e-3/M_H2O + Cp_H2O(288)*log(288.15/273.15)) * comp_f_H2O;

    eta_cyclen = W_m / Q_comb;
    eta_toten  = eta_mec * eta_cyclen;
    eta_cyclex = P_m / (m_g*e_3 - m_a*e_2);
    eta_totex  = P_e / m_c / e_c;
    eta_rotex  = P_m / (m_g*(e_3 - e_4) - m_a*(e_2 - e_1));
    eta_combex = (m_g*e_3 - m_a*e_2) / (m_c*e_c);

    ETA = [eta_cyclen,eta_toten,eta_cyclex,eta_totex,eta_rotex,eta_combex]; % [%]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOSSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    P_1 = (h_1 * m_a) *1e-3;
    P_2 = (h_2 * m_a) *1e-3;
    P_3 = (h_3 * 450) *1e-3;
    P_4 = (h_4 * 450) *1e-3;

    perte_mecen  = (P_m - P_e) *1e-3;
    perte_echen  = m_g * (h_4-h_1) *1e-3;

    perte_mecex  = perte_mecen;
    perte_rotex  = (m_g*(e_3 - e_4) - m_a*(e_2 - e_1) - P_m) *1e-3;
    perte_combex = (m_c*e_c-(m_g*e_3-m_a*e_2)) *1e-3;
    perte_echex  = m_g * (e_4-e_1) *1e-3;

    DATEN = [perte_mecen, perte_echen]; % [kW]
    DATEX = [perte_mecex, perte_rotex, perte_combex, perte_echex]; % [kW]

    COMBUSTION = struct();
    COMBUSTION.LHV    = LHV; % [kJ/kg]
    COMBUSTION.e_c    = e_c; % [kJ/kg]
    COMBUSTION.lambda = lambda; % [mol_air/mol_c]
    COMBUSTION.Cp_g   = Cp_f(400); % [kJ/(kg*K)]
    COMBUSTION.fum = [m_CO2f, m_H2Of, m_O2f, m_N2f]; % [kg/s]

    Cp_g = COMBUSTION.Cp_g;
    P_prim = P_e + (perte_mecex + perte_rotex + perte_combex + perte_echex)*1e3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FIG = [];

if display
    samp = 25;
    p_12 = linspace(p_1,p_2,samp+1); T_12 = [T_1,zeros(1,samp)];
    h_12 = [h_1,zeros(1,samp)]; s_12 = [s_1,zeros(1,samp)];

    p_23 = linspace(p_2,p_3,samp+1); T_23 = linspace(T_2,T_3,samp+1);
    p_23 = p_23(2:end); T_23 = T_23(2:end);
    h_23 = zeros(1,samp); s_23 = zeros(1,samp);

    p_34 = linspace(p_3,p_4,samp+1); T_34 = [T_3,zeros(1,samp)];
    h_34 = [h_3,zeros(1,samp)]; s_34 = [s_3,zeros(1,samp)];

    p_41 = p_1*ones(1,samp+2); T_41 = linspace(T_4,T_1,samp+2);
    h_41 = [h_4,zeros(1,samp),h_1]; s_41 = [s_4,zeros(1,samp),s_1];
    for i = 1:samp
        [T_12(i+1),Cp] = Compression(p_12(i+1)/p_1,T_12(i)*1.5);
        h_12(i+1) = h_1 + Cp*(T_12(i+1) - T_1); % formule enthalpie
        s_12(i+1) = s_1 + (1-eta_PiC) * Cp*log(T_12(i+1)/T_1); % cf. eq 3.15

        Cp = integral(@Cp_air,T_2,T_23(i))*(1-i/samp) + integral(@Cp_f,T_2,T_23(i))*i/samp;
        h_23(i) = h_2 + Cp;%*(T_23(i+1) - T_2);
        s_23(i) = s_2 + Cp*log(T_23(i)/T_2)/(T_23(i)-T_2)...
            - (R_air*(1-i/samp)+R_f*i/samp)*log(p_23(i)/p_2);

        [T_34(i+1),Cp] = Detente(p_34(i+1)/p_3,T_34(i)/1.5);
        h_34(i+1) = h_3 + Cp*(T_34(i+1) - T_3); % formule enthalpie
        s_34(i+1) = s_3 + Cp*log(T_34(i+1)/T_3) - R_f*log(p_34(i+1)/p_3);

        if T_41(i) < 300
            h_41(i+1) = 1.006e3 * (T_41(i+1) - T_0);
            s_41(i+1) = 1.006e3 * log(T_41(i+1)/T_0);
        else 
            Cp = integral(@Cp_air,T_4,T_41(i+1))*i/samp + integral(@Cp_f,T_4,T_41(i+1))*(1-i/samp);
            h_41(i+1) = h_4 + Cp;
            s_41(i+1) = s_4 + Cp*log(T_41(i+1)/T_4) /(T_41(i+1) - T_4);
        end
    end
    T_12 = T_12-273.15; T_23 = T_23-273.15;
    T_34 = T_34-273.15; T_41 = T_41-273.15;

    % HS - GRAPH
        FIG(1) = figure('Color','w'); grid; hold on;
        labels = {'1','2','3','4'};
        X = DAT(4,:); Y = DAT(3,:);
        my_fig = plot(X,Y,'LineStyle','none');
        my_fig.Marker = '.'; my_fig.MarkerSize = 10;
        my_fig.Color = [0 0.4470 0.7410];
        text(X,Y,labels,'FontSize',10,'VerticalAlignment','bottom',...
            'HorizontalAlignment','right');
        my_fig = plot([s_12,s_23,s_34,s_41],[h_12,h_23,h_34,h_41]);
        my_fig.Color = [0 0.4470 0.7410]; my_fig.LineWidth = .75;
        %axis([0 2 0 1800]);
        ylabel('enthalpie, $h$ [kJ/kg]','Interpreter','latex');
        xlabel('entropie, $s$ [kJ/(kg K)]','Interpreter','latex');
        %title('H-S Graph');
    
    % TS - GRAPH
        FIG(2) = figure('Color','w'); grid; hold on;
        labels = {'1','2','3','4'};
        X = DAT(4,:); Y = DAT(1,:);
        my_fig = plot(X,Y,'LineStyle','none');
        my_fig.Marker = '.'; my_fig.MarkerSize = 10;
        my_fig.Color = [0 0.4470 0.7410];
        text(X,Y,labels,'FontSize',12,'VerticalAlignment','bottom','HorizontalAlignment','right')
        my_fig = plot([s_12,s_23,s_34,s_41],[T_12,T_23,T_34,T_41]);
        my_fig.Color = [0 0.4470 0.7410]; my_fig.LineWidth = .75;
        %axis([0 2 0 1800]);
        ylabel('temperature, $T$ [$^\circ$C]','Interpreter','latex');
        xlabel('entropie, $s$ [kJ/(kg K)]','Interpreter','latex');
        %title('T-S Graph');

    % pie chart - energy
        FIG(3) = figure('Color','w'); grid;
        X = [P_e*1e-3 DATEN];
        pie(X,{'','',''});
        legend(sprintf('Puissance effective : %g [MW]',P_e*1e-3),...
            sprintf('Pertes mecaniques : %.3g [MW]',perte_mecen),...
            sprintf('Pertes echappement : %.5g [MW]',perte_echen),...
            'Location','southoutside');
    
    % pie chart - exergy
        FIG(4) = figure('Color','w');
        eta_turb = (h_3-h_4)/(e_3-e_4); pertes_turb = m_g*(e_3-e_4)*(1-eta_turb)*1e-3;
        eta_comp = (e_2-e_1)/(h_2-h_1); pertes_comp = m_a*(e_2-e_1)*(1-eta_comp)*1e-3;
        X = [P_e*1e-3,DATEX(1),pertes_turb,pertes_comp,DATEX(3:4)];
        pie(X,{'','','','','',''});
        legend(sprintf('Puissance effective : %g [MW]', P_e*1e-3),...
            sprintf('Pertes mecaniques : %.3g [MW]', perte_mecex),...
            sprintf('Irr. a la turbine : %.4g [MW]', pertes_turb),...
            sprintf('Irr. au compresseur : %.4g [MW]', pertes_comp),...
            sprintf('Irr. combustion : %.5g [MW]', perte_combex),...
            sprintf('Pertes a la cheminee : %.5g [MW]', perte_echex),...
            'Location','southoutside');
        %title(sprintf('Primary Exergy Flux (%g MW)', P_prim*1e-3));
end


end

