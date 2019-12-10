function Correction_4students(path)


if nargin<2
    noRef=1;% Pas de reference
    if nargin<1
        path='.\GrXX\';%here
        TestGT     = 1;
        TestST     = 1;
        TestCCGT2P = 1;
        TestCCGT3P = 1;
        TestCT     = 1;
    end
end

addpath(path)
fprintf('\n ----------------------------------------------------------------')
fprintf('\n Group : %s',path)
fprintf('\n ----------------------------------------------------------------')

ScoreGT=0;
ScoreST=0;
ScoreCCGT2P=0;
ScoreCCGT3P=0;
ScoreCT=0;
Score=0;  


                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %%%%%%%%%%                     %%%%%%%%%%%
                      %%           Test GT (X/100)            %%
                      %%%%%%%%%%                     %%%%%%%%%%%
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TestGT
    P_e=50e3;%50MW
    options=optionsGT();
    display=1;
    
  
    try % I test the default case:
        [ETA DATEN DATEX DAT MASSFLOW COMBUSTION Cp_g  FIG] = GT();
    catch
        disp('Basic GT failed');
    end
    
    if display
        adress = sprintf('%s%s',path,'RESULTS/GT_');
        for k=1:length(FIG)
            saveas(FIG(k),sprintf('%sFigure_%d',adress,k)); % Using export_fig instead of saveas.
        end
    end
    close all;
    
    
    % Test GT pour different coef de pression (X/100)
    P_e=50e6;
    options=optionsGT(); % YOUR JOB to define a setup here
    options.r=15;
    display=0;
        
    try
        [ETA DATEN DATEX DAT MASSFLOW COMBUSTION Cp_g ] = GT(P_e,options,display);
        %    I compute your score : 
    catch
        disp('GT with assistant inputs failed');
    end
    
    
end

                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %%%%%%%%%%                     %%%%%%%%%%%
                      %%           Test ST (YY/100)           %%
                      %%%%%%%%%%                     %%%%%%%%%%%
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TestST
    P_e=250e3;
    options=optionsST();
    display=1;
    
    try % <=> compare with student results
        [ETA XMASSFLOW DATEN DATEX DAT MASSFLOW COMBUSTION FIG] = ST(P_e,options,display);        
    catch
        disp('Basic ST failed');
        display=0;
    end
    Score=0;
    
    if display
        adress = sprintf('%s%s',path,'RESULTS/ST_');
        for k=1:length(FIG)
            saveas(FIG(k),sprintf('%sFigure_%d',adress,k)); % Using export_fig instead of saveas.
        end
    end
    close all;
    
    %% Another pressure level
    options=optionsST();
    options.p3_hp = 200;%200 bars
    display=1;
    
    try % <=> compare with student results
        [ETA XMASSFLOW DATEN DATEX DAT MASSFLOW COMBUSTION FIG] = ST(P_e,options,display);
    catch
        disp('Basic ST failed');
        display=0;
    end
    
    Score=0;
    
    if display
        adress = sprintf('%s%s',path,'RESULTS/ST_x4_');
        for k=1:length(FIG)
            saveas(FIG(k),sprintf('%sFigure_%d',adress,k)); % Using export_fig instead of saveas.
        end
    end
    close all;
    
   
    
end
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %%%%%%%%%%                     %%%%%%%%%%%
                      %%         Test CCGT3P (15/100)         %%
                      %%%%%%%%%%                     %%%%%%%%%%%
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TestCCGT3P
    P_eg=250e3;%50[MW]
    options=optionsCCGT3P();
    display=1;
    
    try % <=> compare with student results
        [ETA MASSFLOW FIG] = CCGT3P(P_eg,struct(),1);
    catch
        disp('Basic CCGT3P failed');
    end
   
    if display
        adress = sprintf('%s%s',path,'RESULTS/CCGT3P_');
        for k=1:length(FIG)
            saveas(FIG(k),sprintf('%sFigure_%d',adress,k)); % Using export_fig instead of saveas.
        end
    end
    close all;

end



                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %%%%%%%%%%                     %%%%%%%%%%%
                      %%           Test CT (10/100)           %%
                      %%%%%%%%%%                     %%%%%%%%%%%
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if TestCT
    options=optionsCT();
    P_w=444e3;
   
    try
        [DAT_WATER DAT_AIR MASSFLOW] = CoolingTower(P_w,options);
    catch
        disp('Basic CT failed');
    end
end


rmpath(path)

%% Fonctions auxilliaires

    function options=optionsGT()
        options.T_0    =15; %°C Put the value you want 
        % Define your inputs here. They don't need to be all defined.
    end

    function options=optionsST()
        options.nsout       = 7; %Put the value you want 
        % Define your inputs here. They don't need to be all defined.
    end

    function options=optionsCCGT3P()
        options.GT      = optionsGT();
        options.T0      =  15; %°C Put the value you want 
        % Define your inputs here. They don't need to be all defined.
    end
    
    function options=optionsCT()
        options.Tcond  =33; %°C Put the value you want 
        % Define your inputs here. They don't need to be all defined.
    end

end