options = struct();

options.k_mec = .7; % [-] : Shaft losses 
options.T_0 = 0; % [°C] : Reference temperature
options.T_ext = 15; % [°C] : External temperature
options.r = 10; % [-] : Compression ratio
options.k_cc = % [-] : Coefficient of pressure losses due to combustion
%                   chamber
options.T_3 = 1050; %[°C] : Temperature after combustion (before turbine)
option.eta_PiC = .871; %[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for compression
option.eta_PiT = .850; %[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for expansion