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
        P_W=444e3;%200MW_heat
    end
end

if isfield(options,'Tcond')
    Tcond = options.Tcond;
else
    Tcond = 30; %[K]
end

if isfield(options,'Tpinch')
    Tpinch = options.Tpinch;
else
    Tpinch = 4; %[K]
end

if isfield(options,'Tw_out')
    Tw_out = options.Tw_out;
else
    Tw_out = 43; %[K]
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

%% dat_air
[Tdb_in, w_in, phi_in, h_in, Tdp_in, v_in, Twb_in] = Psychrometrics('Tdb',Ta_in,'phi',Phi_atm*100);
[Tdb_out, w_out, phi_out, h_out, Tdp_out, v_out, Twb_out] = Psychrometrics('Tdb',Ta_out,'phi',Phi_out*100);
dat_air = [Tdb_in, Tdb_out; h_in/1000, h_out/1000; w_in, w_out; phi_in/100, phi_out/100]; 

%% Massflows
cp_e = 4.18;
massflow(1) = P_W/(cp_e*(Tw_out-Tw_in));
massflow(3) = (massflow(1)*(XSteam('h_pT',1,Tw_out)-XSteam('h_pT',1,Tw_in)))...
            /((dat_air(2,2)-dat_air(2,1))-(dat_air(3,2)-dat_air(3,1))*XSteam('h_pT',1,Tw_in));
massflow(2) = massflow(3)*(dat_air(3,2)-dat_air(3,1));

%% dat_water
f = @(t_e) massflow(3)*((dat_air(2,2)-dat_air(2,1))-(dat_air(3,2)-dat_air(3,1))*cp_e*t_e)...
           -massflow(1)*cp_e*(Tw_out-t_e);
       
dat_water(1,1) = fsolve(f,1);
dat_water(2,1) = XSteam('h_pT',1,dat_water(1,1));
dat_water(3,1) = massflow(1)-massflow(2);

dat_water(1,2) = Tw_in;
dat_water(2,2) = XSteam('h_pT',1,dat_water(1,2));
dat_water(3,2) = massflow(1);

dat_water(1,3) = Tw_out;
dat_water(2,3) = XSteam('h_pT',1,dat_water(1,3));
dat_water(3,3) = massflow(1);

dat_water(1,4) = Triver;
dat_water(2,4) = XSteam('h_pT',1,dat_water(1,4));
dat_water(3,4) = massflow(2);
end