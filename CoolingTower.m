function [DAT_WATER DAT_AIR MASSFLOW] = CoolingTower(P_w,options)
% COOLINGTOWER is a cooling tower 0D modelisation
% COOLINGTOWER(P_w,options) compute the thermodynamics states for a Cooling
% tower based on several inputs (given in OPTION) and based on a given 
% water power to dissipate P_w.
% It returns the main results. 
%
% INPUTS :
% P_W = Heat power output at the condenser [kW]
% OPTIONS is a structure containing :
%   -options.Tcond  [°C]: Temperature in the condenser
%   -options.Tpinch [°C]: Minimum tempearture pinch between Tw_out and the
%                         condenser temperature.
%   -options.Tw_out [°C]: Cooling water temperature at the condenser outlet
%   -options.Tw_in  [°C]: Cooling water temperature at the condenser inlet
%   -options.Triver [°C]: River temperature 
%   -options.Ta_in  [°C]: Atmospheric air temperature 
%   -options.Ta_out [°C]: Air outlet temperature of cooling tower 
%   -options.Phi_atm [-]: Relative humidity of atmospheric air
%   -options.Phi_out [-]: Maximum relative humidity of air at the cooling 
%                         tower outlet.
%
% OUTPUT :
% MassFlow [kg/s]: Vector containing the different massflow :
%   -massflow(1) : water massflow at the condenser
%   -massflow(2) : additionnal water massflow = water flow evaporated
%   -massflow(3) : air massflow at the cooling tower   
%
%  dat_water = [T_e1       , T_e2       , T_e3       , T_e4;  %[°C]
%               h_e1       , h_e2       , h_e3       , h_e4;  %[kJ/kg]
%               m_e1       , m_e2       , m_e3       , m_e4]; %[kg/s]
% 
%  dat_air   = [Ta_in       , Ta_out  ;  %[°C]
%               ha_in       , ha_out  ;  %[kJ/kg]
%               xa_in       , xa_out  ;  %[kg_water/kg_dry_air]
%               Phia_in     , Phia_out]; %[-] relative humidity
%  
%
% ADDITIONNAL INFORMATIONS
% Water points : 
%       1 : water outlet of cooling tower
%       2 : water just before condenser
%       3 : water just after  condenser
%       4 : water from the river (coming between 1 & 2)
%
% Air points :
%       a_in : air at the cooling tower inlet
%       a_out : air at the cooling tower outlet
%

%% YOUR WORK

if nargin<2
    options=struct();
    if nargin<1
        P_w=200e3;%200MW_heat
    end
end

if isfield(options,'Tcond')
    Tcond = options.Tcond;
else
    Tcond = 35; %[K]
end

if isfield(options,'Tpinch')
    Tpinch = options.Tpinch;
else
    Tpinch = 4; %[K]
end

if isfield(options,'Tw_out')
    Tw_out = options.Tw_out;
else
    Tw_out = 44; %[K]
end

if isfield(options,'Tw_in')
    Tw_in = options.Tw_in;
else
    Tw_in = 30; %[K]
end

if isfield(options,'Triver')
    Triver = options.Triver;
else
    Triver = 15; %[K]
end

if isfield(options,'Ta_in')
    Ta_in = options.Ta_in;
else
    Ta_in = 15; %[K]
end

if isfield(options,'Ta_out')
    Ta_out = options.Ta_out;
else
    Ta_out = 25; %[K]
end

if isfield(options,'Phi_atm')
    Phi_atm = options.Phi_atm;
else
    Phi_atm = 0.8; %[K]
end

if isfield(options,'Phi_out')
    Phi_out = options.Phi_out;
else
    Phi_out = 1; %[K]
end

if(Tw_out-Tcond)<Tpinch
    disp('Problem');
end

%% DAT_AIR
[Tdb_in, w_in, phi_in, h_in, Tdp_in, v_in, Twb_in] = Psychrometrics('Tdb',Ta_in,'phi',Phi_atm*100);
[Tdb_out, w_out, phi_out, h_out, Tdp_out, v_out, Twb_out] = Psychrometrics('Tdb',Ta_out,'phi',Phi_out*100);
DAT_AIR = [Tdb_in, Tdb_out; h_in/1000, h_out/1000; w_in, w_out; phi_in/100, phi_out/100]; 

%% MASSFLOW
cp_e = 4.18;
MASSFLOW(1) = P_w/(cp_e*(Tw_out-Tw_in));
MASSFLOW(3) = (MASSFLOW(1)*(XSteam('hL_T',Tw_out)-XSteam('hL_T',Tw_in)))...
            /((DAT_AIR(2,2)-DAT_AIR(2,1))-(DAT_AIR(3,2)-DAT_AIR(3,1))*XSteam('h_pT',1,Tw_in));
MASSFLOW(2) = MASSFLOW(3)*(DAT_AIR(3,2)-DAT_AIR(3,1));

%% DAT_WATER

f = @(t_e)(MASSFLOW(1)-MASSFLOW(2))*cp_e*t_e -MASSFLOW(1)*cp_e*Tw_in+MASSFLOW(2)*cp_e*Triver;

DAT_WATER(1,1) = fsolve(f,1);
DAT_WATER(2,1) = XSteam('hL_T',DAT_WATER(1,1));
DAT_WATER(3,1) = MASSFLOW(1)-MASSFLOW(2);

DAT_WATER(1,2) = Tw_in;
DAT_WATER(2,2) = XSteam('hL_T',DAT_WATER(1,2));
DAT_WATER(3,2) = MASSFLOW(1);

DAT_WATER(1,3) = Tw_out;
DAT_WATER(2,3) = XSteam('hL_T',DAT_WATER(1,3));
DAT_WATER(3,3) = MASSFLOW(1);

DAT_WATER(1,4) = Triver;
DAT_WATER(2,4) = XSteam('hL_T',DAT_WATER(1,4));
DAT_WATER(3,4) = MASSFLOW(2);
end