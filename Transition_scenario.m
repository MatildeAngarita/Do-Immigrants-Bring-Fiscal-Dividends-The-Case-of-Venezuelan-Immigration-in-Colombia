%%------------------------------------------------------%%
%--------------------------------------------------------%
%                                                        %
%               ----- DSGE: MADURO -----                 %
%                     Junio 2019                         %
%                     DGPM - MHCP                        %
%                                                        %
%%% Documento para la Transición y uso de escenarios   %%%



%______________________________
% Escenario Base 

dynare MADURO_0.mod
A0 = oo_.steady_state;
save A_P0 A0;
%

dynare MADURO_base1.mod
load A_P0 A0;
A1 = oo_.steady_state;

B0 = simult_(A0(:,1),oo_.dr,zeros(100,M_.exo_nbr),1);

%% ______________________________
% Escenario Medio

dynare MADURO_0.mod
A0 = oo_.steady_state;
save A_P0 A0;
%clear 

dynare MADURO_medio1.mod
load A_P0 A0;
A1 = oo_.steady_state;

B2 = simult_(A0(:,1),oo_.dr,zeros(100,M_.exo_nbr),1);

%% ______________________________
% Escenario Alto 

dynare MADURO_0.mod
A0 = oo_.steady_state;
save A_P0 A0;
%clear 

dynare MADURO_alto1.mod
load A_P0 A0;
A1 = oo_.steady_state;

B3 = simult_(A0(:,1),oo_.dr,zeros(100,M_.exo_nbr),1);


 
 
 
 
 
 
% irf_periods=10; %IRF should have 20 periods
% drop_periods=20; %drop 200 periods in simulation as burnin
% irf_replication=1000; %take GIRF average over 1000 periods
% 
% 
% impulse_vec=zeros(1,M_.exo_nbr); %initialize impulse vector to 0
% impulse_vec(strmatch('eps_g',M_.exo_names,'exact'))=1; %impulse of 1 percent of GDP shock to g
% 
% starting_point=A0(:,1); %define starting point of simulations; here: start at steady state and use 200 periods burnin to get to ergodic distribution
% 
% %% initialize shock matrices to 0
% shocks_baseline = zeros(irf_periods+drop_periods,M_.exo_nbr); %baseline
% shocks_impulse = shocks_baseline;
% 
% B = simult_(starting_point,oo_.dr,shocks_impulse,options_.order); %simulation with shock
% 
% 
% 
% 
% 
% 
% 
% 
% %%  eliminate shocks with 0 variance
% i_exo_var = setdiff([1:M_.exo_nbr],find(diag(M_.Sigma_e) == 0 )); %finds shocks with 0 variance
% nxs = length(i_exo_var); %number of those shocks
% chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));%get Cholesky of covariance matrix to generate random numbers
% 
%     shocks_baseline(:,i_exo_var) = randn(irf_periods+drop_periods,nxs)*chol_S; %generate baseline shocks
%     shocks_impulse = shocks_baseline; %use same shocks in impulse simulation
%     shocks_impulse(drop_periods+1,:) = shocks_impulse(drop_periods+1,:)+impulse_vec; %add deterministic impulse
%     y_baseline = simult_(starting_point,oo_.dr,shocks_baseline,options_.order); %baseline simulation
%     y_shock = simult_(starting_point,oo_.dr,shocks_impulse,options_.order); %simulation with shock
% 
