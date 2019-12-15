% Fichier Test
% Etude de l'effet de certains parametres sur le rendements des cycles modelises


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GAS TURBINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eta_toten pour différentes températures en fonction de r
	%{
	options = struct();
	P_e = 230;
	display = 0;

	options.eta_PiC = .90;
	options.eta_PiT = .90;
	options.k_cc = .95;
	options.T_ext = 15;

	r = 3:3:69; L = length(r)
	eta_cyclen = zeros(5,L);
	for i = 1:L
		options.r = r(i);i

		options.T_3 = 900;
		ETA = GT(P_e,options,display);
		eta_cyclen(1,i) = ETA(2);

		options.T_3 = 1050;
		ETA = GT(P_e,options,display);
		eta_cyclen(2,i) = ETA(2);

		options.T_3 = 1200;
		ETA = GT(P_e,options,display);
		eta_cyclen(3,i) = ETA(2);

		options.T_3 = 1350;
		ETA = GT(P_e,options,display);
		eta_cyclen(4,i) = ETA(2);

		options.T_3 = 1500;
		ETA = GT(P_e,options,display);
		eta_cyclen(5,i) = ETA(2);
	end

	figure('Color','w'); grid; hold on;
	my_fig = plot(r,eta_cyclen(1,:)); my_fig.Color = [0 0 153]/255;
	my_fig = plot(r,eta_cyclen(2,:)); my_fig.Color = [76 0 153]/255;
	my_fig = plot(r,eta_cyclen(3,:)); my_fig.Color = [153 0 153]/255;
	my_fig = plot(r,eta_cyclen(4,:)); my_fig.Color = [153 0 76]/255;
	my_fig = plot(r,eta_cyclen(5,:)); my_fig.Color = [153 0 0]/255;

	m = find(eta_cyclen(1,:) == max(eta_cyclen(1,:)));
	my_fig = plot(r(m),eta_cyclen(1,m)); my_fig.Color = [0 0 153]/255;
	my_fig.Marker = '.'; my_fig.MarkerSize = 10;
	m = find(eta_cyclen(2,:) == max(eta_cyclen(2,:)));
	my_fig = plot(r(m),eta_cyclen(2,m)); my_fig.Color = [76 0 153]/255;
	my_fig.Marker = '.'; my_fig.MarkerSize = 10;
	m = find(eta_cyclen(3,:) == max(eta_cyclen(3,:)));
	my_fig = plot(r(m),eta_cyclen(3,m)); my_fig.Color = [153 0 153]/255;
	my_fig.Marker = '.'; my_fig.MarkerSize = 10;
	m = find(eta_cyclen(4,:) == max(eta_cyclen(4,:)));
	my_fig = plot(r(m),eta_cyclen(4,m)); my_fig.Color = [153 0 76]/255;
	my_fig.Marker = '.'; my_fig.MarkerSize = 10;

	legend('900 \circC','1050 \circC','1200 \circC','1350 \circC','1500 \circC');
	axis([0 82 0 .6]);
	xlabel('$p_2/p_1$ [-]','Interpreter','latex');
	ylabel('$\eta_{toten}$ [-]','Interpreter','latex');
	%}

% Diagramme TS pour differents ratio de compression
	%{
	options = struct();
	P_e = 230;
	display = 0;

	r = [3,10,30,50,100]; L = length(r)
	S = []; T = [];
	for i = 1:L
		options.r = r(i);

		[~, ~, ~, ~, ~, ~, ~, FIG] = GT(P_e,options,display);
		S = [S;FIG(1,:)]; T = [T;FIG(2,:)];
	end

	figure('Color','w'); grid; hold on;
	my_fig = plot(S(1,:),T(1,:)); my_fig.Color = [0 0 255]/255;
	my_fig.LineWidth = .7;
	my_fig = plot(S(2,:),T(2,:)); my_fig.Color = [0 128 255]/255;
	my_fig.LineWidth = .7;
	my_fig = plot(S(3,:),T(3,:)); my_fig.Color = [0 204 102]/255;
	my_fig.LineWidth = .7;
	my_fig = plot(S(4,:),T(4,:)); my_fig.Color = [0 255 0]/255;
	my_fig.LineWidth = .7;
	my_fig = plot(S(5,:),T(5,:)); my_fig.Color = [153 255 51]/255;
	my_fig.LineWidth = .7;

	legend('3','10','30','50','100');
	xlabel('entropie, $s$ [kJ/(kg $^\circ$K)]','Interpreter','latex');
	ylabel('temperature, $T$ [$^\circ$C]','Interpreter','latex');
	%}

% Diagramme TS pour differentes temperatures maximales
	%{
	options = struct();
	P_e = 230;
	display = 0;

	T_max = [7,9,11,13,15]*1e2; L = length(T_max)
	S = []; T = [];
	for i = 1:L
		options.T_3 = T_max(i);

		[~, ~, ~, ~, ~, ~, ~, FIG] = GT(P_e,options,display);
		S = [S;FIG(1,:)]; T = [T;FIG(2,:)];
	end

	figure('Color','w'); grid; hold on;
	my_fig = plot(S(1,:),T(1,:)); my_fig.Color = [0 0 153]/255;
	my_fig.LineWidth = .7;
	my_fig = plot(S(2,:),T(2,:)); my_fig.Color = [76 0 153]/255;
	my_fig.LineWidth = .7;
	my_fig = plot(S(3,:),T(3,:)); my_fig.Color = [153 0 153]/255;
	my_fig.LineWidth = .7;
	my_fig = plot(S(4,:),T(4,:)); my_fig.Color = [153 0 76]/255;
	my_fig.LineWidth = .7;
	my_fig = plot(S(5,:),T(5,:)); my_fig.Color = [153 0 0]/255;
	my_fig.LineWidth = .7;

	legend('700','900','1100','1300','1500');
	xlabel('entropie, $s$ [kJ/(kg $^\circ$K)]','Interpreter','latex');
	ylabel('temperature, $T$ [$^\circ$C]','Interpreter','latex');
	%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STEAM TURBINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eta_toten en fonction de la pression maximale
	%{
	options = struct();
	P_e = 288e3;
	display = 0;

	options.T_max = 560;
	p_max = 100:20:400; eta_totex = zeros(length(p_max),1);
	for i = 1:length(p_max)
		options.p3_hp = p_max(i)
		ETA = ST(P_e,options,display);
		eta_totex(i) = ETA(4);
	end

	figure('Color','w'); grid;
	my_fig = plot(p_max,eta_totex);
	my_fig.LineWidth = .7;
	legend('560 [\circC]');
	xlabel('$p_{max}$ [bar]','Interpreter','latex');
	ylabel('$\eta_{totex}$ [\%]','Interpreter','latex');
	%}

% eta_toten en fonction de la pression maximale
	% {
	options = struct();
	P_e = 288e3;
	display = 1;
	options.p3_hp = 200;
	T_max = 300:10:600;
	options.lam_chaud = 1;
	options.lam_exch = 1;
	options.nsout = 4;
	options.drumFlag = 0;
	eta_totex = zeros(length(T_max),1);
	for i = 1:length(T_max)
		options.T_max = T_max(i)
		ETA = ST(P_e,options,display);
		eta_totex(i) = ETA(4);
	end

	figure('Color','w'); grid;
	my_fig = plot(T_max,eta_totex);
	my_fig.LineWidth = .7; 
	legend('200 [bar]');
	xlabel('$T_{max}$ ^[\circC]','Interpreter','latex');
	ylabel('$\eta_{totex}$ [\%]','Interpreter','latex');
	%}

% Rendement exergétique total en fonction du nombre de soutirages
	%{
	options = struct();
	P_e = 288e3;
	display = 0;
	options.p3_hp = 400;
	options.T_max = 600;
	options.nsout = 15;
	options.TpinchEx = 0;
	options.TpinchSub = 3;
	options.eta_SiT = .927;
	options.lam_exch = 1;

	reheat = 1:10;
	eta_totex = zeros(length(reheat),2);
	for i = 1:length(reheat)
		options.reheat = reheat(i);
		options.drumFlag = 1 
		ETA = ST(P_e,options,display);
		eta_totex(i,1) = ETA(4);
		%options.drumFlag = 0
		%ETA = ST(P_e,options,display);
		%eta_totex(i,2) = ETA(4);
	end

	figure('Color','w'); grid; hold on;
	p = plot(reheat,eta_totex(:,1),'.-'); p.MarkerSize = 10;
	%p = plot(reheat,eta_totex(:,2),'.-'); p.MarkerSize = 10;
	legend('avec tambour')%,'sans tambour');
	xlabel('Nombre de rechauffes, $n_{reheat}$','Interpreter','latex');
	ylabel('$\eta_{totex}$ [\%]','Interpreter','latex');
	%}

	%{
	options = struct();
	P_e = 288e3;
	display = 0;
	options.reheat = 1;
	options.TpinchEx = 3;

	nsout = 2:20; eta_totex = zeros(length(nsout),2);
	for i = 1:length(nsout)
		options.nsout = nsout(i);
		options.drumFlag = 1
		ETA = ST(P_e,options,display);
		eta_totex(i,1) = ETA(4);
		options.drumFlag = 0
		ETA = ST(P_e,options,display);
		eta_totex(i,2) = ETA(4);
	end

	figure('Color','w'); grid; hold on;
	p = plot(nsout,eta_totex(:,1),'.-'); p.MarkerSize = 10;
	p = plot(nsout,eta_totex(:,2),'.-'); p.MarkerSize = 10;
	legend('avec tambour','sans tambour');
	xlabel('Nombre de soutirages, $n_{sout}$','Interpreter','latex');
	ylabel('$\eta_{totex}$ [\%]','Interpreter','latex');
	%}