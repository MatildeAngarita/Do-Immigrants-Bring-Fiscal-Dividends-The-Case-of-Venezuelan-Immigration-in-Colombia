%%------------------------------------------------------%%
%--------------------------------------------------------%
%    Do immigrants  Bring Fiscal Dividends? The          %
%    Case of Venexuelan Immigration in Colombia          %
%                                                        %
%     Matilde Angarita  -  Marcela De Castro             %
%     Juan Camilo Santaella - Oscar Valencia             %
%%------------------------------------------------------%%
%--------------------------------------------------------%
%%% This code calls and run the three scenarios
    % Note: Dynare is needed, we use version 4.5.7
   
%% --------------------  
    %   Low Immigration Scenario
    clear all
    
dynare init_scenario.mod
A0 = oo_.steady_state;
save A_P0 A0;
clear 

dynare Low_scenario.mod
load A_P0 A0;
AR1 = oo_.steady_state;

Low_results = simult_(A0(:,1),oo_.dr,zeros(100,M_.exo_nbr),1);
save results_1 Low_results 
%% --------------------  
        %  Medium Immigration Scenario  
         clear all
         
dynare init_scenario.mod
A0 = oo_.steady_state;
save A_P0 A0;
clear 

dynare medium_scenario.mod
load A_P0 A0;
AR2 = oo_.steady_state;

Med_results = simult_(A0(:,1),oo_.dr,zeros(100,M_.exo_nbr),1);
save results_2 Med_results 
%% --------------------  
        %  High Immigration Scenario  
         clear all
        
dynare init_scenario.mod
A0 = oo_.steady_state;
save A_P0 A0;
clear 

dynare High_scenario.mod
load A_P0 A0;
AR3 = oo_.steady_state;

Hig_results = simult_(A0(:,1),oo_.dr,zeros(100,M_.exo_nbr),1);
save results_3 Hig_results  
