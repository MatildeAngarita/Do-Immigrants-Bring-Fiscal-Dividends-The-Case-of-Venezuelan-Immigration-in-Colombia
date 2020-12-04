%%------------------------------------------------------%%
%--------------------------------------------------------%
%    Do immigrants  Bring Fiscal Dividends? The          %
%    Case of Venexuelan Immigration in Colombia          %
%                                                        %
%     Matilde Angarita  -  Marcela De Castro             %
%     Juan Camilo Santaella - Oscar Valencia             %
%%------------------------------------------------------%%
%--------------------------------------------------------%
%%% Model code for the Low Immigration scenario

%% CALIBRATION %%
%%%%%%%%%%%%%%%%%

init= 5*ones(4,1);
[xs,fval] = fsolve(@solve_model_init,init); % Call solver

pmg_s_cali = xs(1);
pmg_u_cali = xs(2);
c_h_cali   = xs(3);
c_f_cali   = xs(4);

calibra_model;
load q_h_s_cali q_h_u_cali q_f_s_cali q_f_u_cali;


options = optimoptions('fsolve','Display','iter','TolFun',1e-6,'TolX',1e-8);
initend= 1*ones(8,1);

[xend,fval] = fsolve(@solve_model_endval_high,initend); % Call solver


pmg_s_end      = real(xend(1));
pmg_u_end      = real(xend(2));
c_h_end        = real(xend(3));
c_f_end        = real(xend(4));
miu_h_s_end    = real(xend(5));
miu_h_u_end    = real(xend(6));
miu_f_s_end    = real(xend(7));
miu_f_u_end    = real(xend(8));
% 

%%


%--------------------%
%      Variables     %
%--------------------%


    var 
         g_h         g_h_s       g_h_u      g_f_s       N_h_s      N_h_u       N_f_s        N_f_u       g_f_u         N_h
         g_f          g          N_f         N          mig        def          E            A          H_h_s        H_h_u     
         H_f_s       H_f_u       H_u         H_s        H_h        H_f        miu_h_s      miu_h_u     miu_f_s      miu_f_u 
        theta_h_s   theta_h_u   theta_f_s   theta_f_u   V_h_s      V_h_u       V_f_s        V_f_u       w_h_s        w_h_u   
         w_f_s       w_f_u       r_h         r_f        C_h        C_f           K            C           I            Y    
         h_h_s       h_h_u       h_f_s       h_f_u      h_h        h_u          h_s          h_f         c_h          c_f       
           c          k           B          B_h        B_f        pi          L_h_s         L_h_u       L_f_s       L_f_u
         L_u         L_s          i          y         pmg_k       v_h_s       v_h_u         v_f_s       v_f_u       pmg_s   
         pmg_u       phi_h      phi_f       omega       L_h        L_f       nmono_h_s     nmono_h_u    nmono_f_s   nmono_f_u   
          u          u_h         u_f        u_h_s      u_h_u       u_f_s       u_f_u        recaudo      g_k          Ak     
        miubar_s    miubar_u    U_h_s       U_h_u      U_f_s       U_f_u        U_h           U_f        U_s          U_u        
          U          u_s        u_u       recaudo_c   recaudo_h  deficit_y  gasto_flujo  gasto_flujo_y  Y_anual     def_anual  
      gasto_anual gasto_flujo_y_anual  d_y_anual d_c_anual  d_i_anual  d_recaudo_anual ratio_s   ratio_u ratio_h ratio_f
              recaudo_h_mig recaudo_h_loc recaudo_c_mig recaudo_c_locales
         recaudo_c_y 
     recaudo_h_y

      ;
         
         
  % exogenous shocks
    varexo e_ptf e_ak e_shock e_mig e_smiu e_umiu ;
    
    
    
%--------------------%
%      Parameters    %
%--------------------%    
    
   parameters 
   
   %Desutilidad del trabajo 
      psi 
    %Impuestos  
      tau_c
      tau_h_s 
      tau_h_u
      tau_f_s 
      tau_f_u
    % Subsidios
      b_h_s 
      b_h_u 
      b_f_s 
      b_f_u 
  
    %Depreciasión
      delta 
      kappa
    % Tasa de separación  
      x     
    % Participación del Unskilled en la función de producción  
      sigma 
      alpha 
    % Participación del kapital 
      rho 
    % Elasticidad de sustitución entre capital y trabajo calificado  
      ipsilon 
    % Productividad del capital   
      beta
    % Elasticidad del mercado laboral  
      eta 
    %Elasticidad de Frish
      gamma 
    % Negociación (poder de la firma)
      lambda_h_s
      lambda_h_u 
      lambda_f_s
      lambda_f_u 
    % crecimientos
      g_h_s_bar
      g_h_u_bar 
      g_f_s_bar 
     % g_f_u_bar
      g_h_bar
      
      rhopob
     N_h_s_bar   
     N_h_u_bar  
     N_f_s_bar 
     N_f_u_bar  
     
     N_h_s_end   
     N_h_u_end  
     N_f_s_end 
     N_f_u_end  
    % Migración
      rhomig
      migbar 
   
      q_h_s
      q_h_u 
      q_f_s
      q_f_u
      
      pmg_s_bar
      pmg_u_bar
      c_h_bar
      c_f_bar
   
   
      %Choques
      edu
      rho_edu
      
      Abar
      za 
      ze 
      
      Akbar
      rho_ak
      
      miu_ss_s
      rho_smiu
      
      miu_ss_u
      rho_umiu
      b_flujo
      phi_g
     G
   ;

%--------------------%
%    Calibration     %
%--------------------%   
         %Desutilidad del trabajo 
      psi = 0.15; 
    %Impuestos  
      tau_c = 0.08;
      tau_h_s = 0.04;
      tau_h_u = 0.04;
      tau_f_s = 0.04;
      tau_f_u = 0.04;
    % Subsidios
%       b_h_s = 0.88;
%       b_h_u = 0.88;
%       b_f_s = 0.88;
%       b_f_u = 0.88;
      b_h_s = 0;
      b_h_u = 0;
      b_f_s = 0;
      b_f_u = 0;

    %Depreciasión
      delta = 0.025;
      kappa = 8;
    % Tasa de separación  
      x     = 0.11;
    % Participación del Unskilled en la función de producción  
      sigma = 0.4047;
      alpha = -0.2;
    % Participación del kapital 
      rho = 0.7331;
    % Elasticidad de sustitución entre capital y trabajo calificado  
      ipsilon = -0.2;
    % Productividad del capital   
      beta = 0.985;
    % Elasticidad del mercado laboral  
      eta = 0.5;
    %Elasticidad de Frish
      gamma = 0.2;
    % Negociación (poder de la firma)
      lambda_h_s = 0.85;
      lambda_h_u = 0.85;
      lambda_f_s = 0.5;
      lambda_f_u = 0.5; 
    % crecimientos
      g_h_s_bar = 1;
      g_h_u_bar = 1;
      g_f_s_bar = 1;
      g_h_bar = 1;
      
      rhopob=0.23;
            
     
     N_h_s_end   =  10554.128;   
     N_h_u_end   =  36661.677;  
     N_f_s_end   =  152.717; 
     N_f_u_end   =  1022.026 + 4877.550;   

     N_h_s_bar   =  10554.128;   
     N_h_u_bar   =  36661.677;  
     N_f_s_bar   =  152.717; 
     N_f_u_bar   =  N_f_u_end; 
     
   % Migración
      rhomig = 0.8;
      migbar = 0;
     
    % Educación 
      edu = 1;
      rho_edu = 0.2;
          %Shocks
      Abar       = 0.6;
      za         = 0.8;
      ze         = 1;
      
      Akbar      =0.3;
      rho_ak     =0.8;

      miu_ss_s   =0.6;
      rho_smiu   =0.8;
      
      miu_ss_u   =0.6; 
      rho_umiu   =0.8;
          %Shocks
          
      gasto_ven = 0;    

%       b_flujo = 950;
       b_flujo = 0.299;

      phi_g   = 0.9;  
       G=800;
    
      
      %LLAMAR !NO OLVIDAR!
     q_h_s = q_h_s_cali;
     q_h_u = q_h_u_cali;
     q_f_s = q_f_s_cali;
     q_f_u = q_f_u_cali;
     
     pmg_s_bar = pmg_s_cali;
     pmg_u_bar = pmg_u_cali;
     c_h_bar = c_h_cali;
     c_f_bar = c_f_cali;

     
  
      
     
%%      
%----------------------------%
%            MODEL           %
%----------------------------%    
  
model;


% Crecimiento poblacional
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
   % Ec. 1 - Total  
     g_h = N_h/N_h(-1); 
     
   % Ec. 2 - tasa de crecimiento home skilled
     g_h_s = N_h_s/N_h_s(-1);    
      
   % Ec. 3 - tasa de crecimiento home unskilled
     g_h_u = N_h_u/N_h_u(-1); 
      
   % Ec. 4 - tasa de crecimiento home
     g_f_s = N_f_s/N_f_s(-1);
  
   % Ec. 5 - Home- Skilled
     N_h_s=(1-rhopob)*N_h_s(-1)+rhopob*N_h_s_bar;
     
   % Ec. 6 - Home - Unskilled
     N_h_u=(1-rhopob)*N_h_u(-1)+rhopob*N_h_u_bar;
     
   % Ec. 7 - Foreign - Skilled
     N_f_s=(1-rhopob)*N_f_s(-1)+rhopob*N_f_s_bar;
      
   % Ec. 8 - Crecimiento de Población extranjera
     g_f_u = N_f_u/N_f_u(-1);
     
   % Ec. 9 - Foreign - Unskilled
     N_f_u=(1-rhopob)*N_f_u(-1)+rhopob*N_f_u_end+mig;    
            
   % Ec. 10 - Home
     N_h = N_h_s + N_h_u;   
     
   % Ec. 11 - Foreign
     g_f = N_f/N_f(-1); 
     
   % Ec. 12 - Total
     g = N/N(-1); 
     
   % Ec. 13 - Foreign agregado
     N_f = N_f_s + N_f_u;
   
   % Ec. 14 - total agregado
     N= N_h + N_f;
  
   % Ec. 15 - Choque
     mig = rhomig*migbar + (1-rhomig)*mig(-1) + e_mig;

     
     
     
% Gobierno
%------------------------------------------------------------------------------------------------------------------------------------------------
    % Ec. 16 - Déficit
     def = tau_c*C + tau_h_s*w_h_s*H_h_s(-1) + tau_h_u*w_h_u*H_h_u(-1) + tau_f_s*w_f_s*H_f_s(-1) + tau_f_u*w_f_u*H_f_u(-1) -  b_h_s*(N_h_s-H_h_s)- b_h_u*(N_h_u-H_h_u(-1))- b_f_s*(N_f_s-H_f_s(-1))- b_f_u*(N_f_u-H_f_u(-1))- gasto_flujo-G;
     
      
    % Ec. 17 Educación
      E = (1-rho_edu)*E(-1) + edu*g_f*E(-1)*e_shock;
      
    % Ec. 18 PTF
      A = (1-za)*Abar + za*A(-1) +e_ptf;


 
% Mercado Laboral
%------------------------------------------------------------------------------------------------------------------------------------------------      
    % Ec. 19 - Dinámica del empleo (Home-skilled)
        H_h_s = (1-x)*H_h_s(-1) + miu_h_s*V_h_s;
        
    % Ec. 20 - Dinámica del empleo (Home- unskilled)
        H_h_u = (1-x)*H_h_u(-1) + miu_h_u*V_h_u;
        
    % Ec. 21 - Dinámica del empleo (Foreign - skilled)
        H_f_s = (1-x)*H_f_s(-1) + miu_f_s*V_f_s;
        
    % Ec. 22 - Dinámica del empleo (Foreign - unskilled)
       H_f_u = (1-x)*H_f_u(-1) + miu_f_u*V_f_u;
       
    % Ec. 23 - trabajo unskilled
      H_u= H_h_u + H_f_u;

    % Ec. 24 - trabajo unskilled
      H_s= H_h_s + H_f_s;  
     
    % Ec. 25 - trabajo home
      H_h= H_h_s + H_h_u;

    % Ec. 26 - trabajo foreign
      H_f= H_f_s + H_f_u;      
 
    % Ec. 27 - Euler de las firmas (Home-skilled)
      q_h_s/miu_h_s  =  beta*(pmg_s(+1)-w_h_s(+1)+(q_h_s/miu_h_s(+1))*(1-x));
     
    % Ec. 28 - Euler de las firmas (Home- unskilled)
      q_h_u/miu_h_u  =  beta*(pmg_u(+1)-w_h_u(+1)+(q_h_u/miu_h_u(+1))*(1-x));
      
    % Ec. 29 - Euler de las firmas (Foreign - skilled)
      q_f_s/miu_f_s  = beta*(pmg_s(+1)-w_f_s(+1)+(q_f_s/miu_f_s(+1))*(1-x));
    
    % Ec. 30 - Euler de las firmas (Foreign - unskilled)
      q_f_u/miu_f_u  = beta *(pmg_u(+1)-w_f_u(+1)+(q_f_u/miu_f_u(+1))*(1-x));
     
    % Ec. 31 - Probabilidad de usar una vacante  (Home-skilled)
      miu_h_s = miubar_s*(theta_h_s)^(-eta);

    % Ec. 32 - Probabilidad de usar una vacante  (Home- unskilled)
      miu_h_u = miubar_u*(theta_h_u)^(-eta);
     
    % Ec. 33 - Probabilidad de usar una vacante  (Foreign - skilled)
      miu_f_s = miubar_s*(theta_f_s)^(-eta);
        
    % Ec. 34 - Probabilidad de usar una vacante  (Foreign - unskilled)
      miu_f_u = miubar_u*(theta_f_u)^(-eta);
      
    % Ec. 35 - Holgura del mercado (Home-skilled)
      theta_h_s = V_h_s/(N_h_s - H_h_s(-1));

    % Ec. 36 - Holgura del mercado (Home- unskilled)
      theta_h_u = V_h_u/(N_h_u - H_h_u(-1));
      
    % Ec. 37 - Holgura del mercado (Foreign - skilled)
      theta_f_s = V_f_s/(N_f_s - H_f_s(-1));
    
    % Ec. 38 - Holgura del mercado (Foreign - unskilled)
      theta_f_u = V_f_u/(N_f_u - H_f_u(-1));
     
    % Ec. 39 - Salarios de equilibrio (Home- Skilled)
      w_h_s = lambda_h_s*(pmg_s +(theta_h_s*q_h_s))+ (1-lambda_h_s)*((psi*((H_h(-1)/N_h)^gamma))*(C_h/N_h)*(1+tau_c)+b_h_s)/(1-tau_h_s);
   
    % Ec. 40 - Salarios de equilibrio (Home- Unskilled)
      w_h_u = lambda_h_u*(pmg_u +(theta_h_u*q_h_u))+ (1-lambda_h_u)*((psi*((H_h(-1)/N_h)^gamma))*(C_h/N_h)*(1+tau_c)+b_h_u)/(1-tau_h_u);
 
    % Ec. 41 - Salarios de equilibrio (Foreign- Skilled)
      w_f_s = lambda_f_s*(pmg_s +(theta_f_s*q_f_s))+ (1-lambda_f_s)*((psi*((H_f(-1)/N_f)^gamma))*(C_f/N_f)*(1+tau_c)+b_f_s)/(1-tau_f_s);

    % Ec. 42 - Salarios de equilibrio (Foreign- Unskilled)
      w_f_u = lambda_f_u*(pmg_u +(theta_f_u*q_f_u))+ (1-lambda_f_u)*((psi*((H_f(-1)/N_f)^gamma))*(C_f/N_f)*(1+tau_c)+b_f_u)/(1-tau_f_u);

     
% Hogares
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    
   % Ec. 43 - Euler (locales)
     c_h(+1)/c_h = beta*(r_h(+1)/g_h(+1));
     
   % Ec. 44 - Euler (foráneos)
     c_f(+1)/c_f = beta*(r_f(+1)/g_f(+1));
    
   % Capital
   %-----------------------------
  
     
   % Ec. 45 -  Euler del capital
     (1-(kappa*((i(+1)/k(+1))-delta)))/(1-kappa*((i/k)-delta)) = beta*((1/g(+1))*((pmg_k(+1)*(1-kappa*((i(+1)/k(+1))-delta))) + (1-delta) + kappa*((i(+1)/k(+1))-delta)*(i(+1)/k(+1))-(kappa/2)*(((i(+1)/k(+1))-delta)^2)));
     
   % Ec. 46 - Producto marginal del capital
     pmg_k =((y^(1-alpha))*(A^alpha))*(1-sigma)*((rho*((Ak*k)^ipsilon) + (1-rho)*(L_s^ipsilon))^((alpha/ipsilon)-1))*rho*(Ak^ipsilon)*((k)^(ipsilon-1));        
 
   % Ec. 47 - Restricción presupuestal (Home)
     C_h*(1+tau_c) + B_h = (w_h_s*H_h_s(-1))*(1-tau_h_s) + (w_h_u*H_h_u(-1)*(1-tau_h_u)) + b_h_s*(N_h_s-H_h_s(-1)) +  b_h_u*(N_h_u-H_h_u(-1)) + r_h*B_h(-1) + pi;
     
   % Ec. 48 - Restricción presupuestal (Foreign)
     C_f*(1+tau_c) + B_f = (w_f_s*H_f_s(-1))*(1-tau_f_s) + (w_f_u*H_f_u(-1)*(1-tau_f_u)) + b_f_s*(N_f_s-H_f_s(-1)) +  b_f_u*(N_f_u-H_f_u(-1)) + r_f*B_f(-1) + gasto_flujo;
     
   % Ec. 49 - Agregación de los bonos
%      B = B_f + B_h;
   
   % Ec. 50 - Agregación del consumo
     C = C_h + C_f;

   % Ec. 51 - Bonos Home
      B_h = phi_g*B;
      
   % Ec. 52 - Bonos Home
      B_f = (1-phi_g)*B;
   
   % Ec. 53 - Dinámica capital Home
      K= I - ((kappa/2)*(((I/K(-1))-delta)^2)*K(-1)) +(1-delta)*K(-1);
    
   % Ec. 54 - Beneficios de las firmas
      pi = Y - w_h_s*H_h_s - w_h_u*H_h_u - w_f_s*H_f_s - w_f_u*H_f_u - V_h_s*q_h_s - V_h_u*q_h_u - V_f_s*q_f_s - V_f_u*q_f_u - I;

           
    
% Firmas 
%------------------------------------------------------------------------------------------------------------------------------------------------

    % Ec. 55 - Función de producción
      y   =  A*(sigma*((L_u)^alpha) + (1-sigma)*(rho*((Ak*k)^ipsilon)+(1-rho)*((L_s)^ipsilon))^(alpha/ipsilon))^(1/alpha) ;

      
   % Variables percapita
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
    % Ec. 56 - trabajadores skilled home 
      h_h_s = H_h_s(-1)/N_h;   

    % Ec. 57 - trabajadores unskilled home 
      h_h_u = H_h_u(-1)/N_h;
      
    % Ec. 58 - trabajadores skilled foreign 
      h_f_s = H_f_s(-1)/N_f;       
        
    % Ec. 59 - trabajadores unskilled foreign 
      h_f_u = H_f_u(-1)/N_f;  
      
    % Ec. 60 - Empleo per cápita
      h_h = H_h(-1) /N_h;     

    % Ec. 61 - Empleo per cápita
      h_f = H_f(-1) /N_f;
            
    % Ec. 62 - Empleo per cápita
      h_s = H_s(-1) /N;
      
    % Ec. 63 - Empleo per cápita
      h_u = H_u(-1) /N;
        
    % Ec. 64 - Consumo home
      c_h = C_h/N_h;   
     
    % Ec. 65 - Consumo foreign
      c_f = C_f/N_f;   
       
    % Ec. 66 - Consumo
      c = C/N; 

    % Ec. 67 - Capital
      k = K(-1)/N;  
   
    % Ec. 68 - Bonos
      B = -def +r_h*B_h(-1) + r_f*B_f(-1);

    % Ec. 69 - Recaudo de impuesto al consumo
      recaudo_c = tau_c*C;
        
    % Ec. 70 - Recaudo de impuesto al trabajo
      recaudo_h = tau_h_s*H_h_s(-1)*w_h_s + tau_h_u*H_h_u(-1)*w_h_u +tau_f_s*H_f_s(-1)*w_f_s +tau_f_u*H_f_u(-1)*w_f_u; 

    % Ec- 71 - Recaudo del gobierno
          recaudo = tau_c*C + tau_h_s*H_h_s(-1)*w_h_s + tau_h_u*H_h_u(-1)*w_h_u +tau_f_s*H_f_s(-1)*w_f_s +tau_f_u*H_f_u(-1)*w_f_u; 
      
    % Ec. 72 - trabajadores skilled home (percápita)
      L_h_s = H_h_s(-1)/N;   

    % Ec. 73 - trabajadores unskilled home (percápita)
      L_h_u = H_h_u(-1)/N;
      
    % Ec. 74 - trabajadores skilled foreign (percápita)
      L_f_s = H_f_s(-1)/N;       
        
    % Ec. 75 - trabajadores unskilled foreign (percápita)
      L_f_u = H_f_u(-1)/N;   
      
    % Ec. 76 - trabajadores unskilled (percápita)
      L_u = H_u(-1)/N;   

    % Ec. 77 - trabajadores skilled(percápita)
      L_s = H_s(-1)/N;   

    % Ec. 78 - Inversión
      i = I/N;
    
    % Ec. 79 - Producto per cápita
      y = Y/N;
       
    % Ec. 80 - definición de vacante percápita
      v_h_s = V_h_s/N;
     
    % Ec. 81 - definición de vacante percápita
      v_h_u = V_h_u/N;
      
    % Ec. 82 - definición de vacante percápita
      v_f_s = V_f_s/N;
      
    % Ec. 83 - definición de vacante percápita
      v_f_u = V_f_u/N;

      
% Ecuaciones auxiliares
%-----------------------------------------------------------------------------------------------------------------------------------------------------------
       
    % Ec. 84 - Producto marginal skilled
      pmg_s = ((y^(1-alpha))*(A^alpha))*((1-sigma)*((rho*((Ak*k)^ipsilon)+(1-rho)*(L_s^ipsilon))^((alpha/ipsilon)-1))*(1-rho)*(L_s^(ipsilon-1)));

    % Ec. 85 - Producto marginal unskilled
      pmg_u = ((y^(1-alpha))*(A^alpha))*sigma*(L_u^(alpha -1));     
      
    % Ec. 86 - proporción de skilled
      phi_h = N_h_s/ N_h;

    % Ec. 87 - proporción de unskilled
      phi_f = N_f_s/ N_f;
            
    % Ec. 88 - proporción de home
      omega = N_h/N;
      
    % Ec. 89  
      L_h = L_h_s + L_h_u;
      
    % Ec. 90  
      L_f = L_f_s + L_f_u; 
              
    % Ec. 91 - Población fijada en Estado Estacionario
      nmono_h_s = N_h_s/N;

    % Ec. 92 - Población fijada en Estado Estacionario
      nmono_h_u = N_h_u/N;
      
    % Ec. 93 - Población fijada en Estado Estacionario
      nmono_f_s = N_f_s/N;
      
    % Ec. 94 - Población fijada en Estado Estacionario
      nmono_f_u = N_f_u/N;  
    
    % Ec. 95 Tasa de desempleo agregada
      u = U/N;
    
    % Ec. 96 Tasa de desempleo home
      u_h = U_h/N_h;
      
    % Ec. 97 Tasa de desempleo foreign
      u_f = U_f/N_f;
    
    % Ec. 98 Tasa de desempleo Home Skilled
      u_h_s = U_h_s/N_h_s;
    
    % Ec. 98 Tasa de desempleo Home Unkilled
      u_h_u = U_h_u/N_h_u;
    
    % Ec. 100 Tasa de desempleo Foreign Skilled
      u_f_s = U_f_s/N_f_s;
    
    % Ec. 101 Tasa de desempleo Foreign Unkilled
      u_f_u = U_f_u/N_f_u;

    % Ec- 102 - Shock de productividad del capital
      Ak = (1-rho_ak)*Akbar + rho_ak*Ak(-1) +e_ak;

    % Ec- 103 - Shock de productividad matching Skilled
      miubar_s = (1-rho_smiu)*miu_ss_s + rho_smiu*miubar_s(-1) +e_smiu;

    % Ec- 104 - Shock de productividad matching Unskilled
      miubar_u = (1-rho_umiu)*miu_ss_u + rho_umiu*miubar_u(-1) +e_umiu;

    % Ec. 105 Desempleo Home Skilled
      U_h_s = N_h_s - H_h_s(-1);
    
    % Ec. 106 Desempleo Home Unkilled
      U_h_u = N_h_u - H_h_u(-1);
    
    % Ec. 107 Desempleo Foreign Skilled
      U_f_s = N_f_s - H_f_s(-1);
    
    % Ec. 108 Desempleo Foreign Unkilled
      U_f_u = N_f_u - H_f_u(-1);
    
    % Ec. 109 Desempleo Home 
      U_h   = U_h_s + U_h_u;
    
    % Ec. 110 Desempleo Foreign 
      U_f   = U_f_s + U_f_u;
    
    % Ec. 111 Desempleo Skilled
      U_s = U_h_s + U_f_s;
    
    % Ec. 112 Desempleo Unskilled
      U_u = U_h_u + U_f_u;
    
    % Ec. 113 Desempleo Total
      U = U_h + U_f;
     
    % Ec. 114 Tasa de desempleo Skilled
      u_s = U_s/(N_h_s + N_f_s);
    
    % Ec. 115 Tasa de desempleo Unskilled
      u_u = U_u/(N_h_u + N_f_u);
           
% ---- Déficit y gasto 

    % Ec. 117 - Déficit como % del PIB
%       deficit_y = def/Y;
       deficit_y = def_anual/Y_anual;

 
%    gasto_flujo = b_flujo*(g_f_u-1);
   gasto_flujo = b_flujo*(((N_f_u + N_f_u(-1)+ N_f_u(-2)+ N_f_u(-3))/4)-((N_f_u(-1) + N_f_u(-2)+ N_f_u(-3)+ N_f_u(-4))/4));

   gasto_flujo_y = gasto_flujo/Y;
   
   %Ec. 20
   Y_anual = Y + Y(-1) + Y(-2) + Y(-3);   
  

   
   def_anual = def + def(-1) + def(-2) +def(-3);
   gasto_anual = gasto_flujo + gasto_flujo(-1) + gasto_flujo(-2) +gasto_flujo(-3);
   
   gasto_flujo_y_anual = gasto_anual/Y_anual;
   
     d_y_anual = (Y_anual/Y_anual(-4))-1;
     d_c_anual = ((C + C(-1) + C(-2) + C(-3))/(C(-4) + C(-5) + C(-6) + C(-7)))-1;
  
     
   %Ec. 26
     d_i_anual = (I + I(-1) + I(-2) + I(-3))/(I(-4) + I(-5) + I(-6) + I(-7))-1; 
     d_recaudo_anual = (recaudo + recaudo(-1) + recaudo(-2) + recaudo(-3))/(recaudo(-4) + recaudo(-5) + recaudo(-6) + recaudo(-7))-1;
     
     ratio_s = V_f_s/V_h_s;
     ratio_u = V_f_u/V_h_u;
     ratio_f = V_f_s/V_f_u;
     ratio_h = V_h_s/V_h_u; 
         
    % Ec. 132 - Tasa de crecimiento del capital
      g_k = k/k(-1);   
     
     recaudo_h_mig = (tau_f_u*H_f_u*w_f_u +tau_f_s*H_f_s*w_f_s)/Y ;
     recaudo_h_loc = (tau_h_s*H_h_s*w_h_s + tau_h_u*H_h_u*w_h_u)/Y;
     recaudo_c_mig = (tau_c*C_f)/Y;
     recaudo_c_locales = (tau_c*C_h)/Y;


     recaudo_c_y = recaudo_c/Y;
     recaudo_h_y = recaudo_h/Y;

     end;
%----------------------------------------%
%           FINAL STEADY STATE           %
%----------------------------------------%

initval;

     
    pmg_s       = pmg_s_end;
    pmg_u       = pmg_u_end;
    c_h         = c_h_end;
    c_f         = c_f_end;
    miu_h_s     = miu_h_s_end;
    miu_h_u     = miu_h_u_end;
    miu_f_s     = miu_f_s_end;
    miu_f_u     = miu_f_u_end;
 
    
 
% Tasas de crecimiento
     g            = 1; %1
     g_h          = g_h_bar; %2
     g_f          = 1; %3
     g_h_s        = g_h_s_bar; %4
     g_h_u        = g_h_u_bar; %5
     g_f_s        = g_f_s_bar; %6
     g_f_u        = 1; %7

% Población (Exógena)
     N_h_s        =  N_h_s_end ;   
     N_h_u        =  N_h_u_end ;  
     N_f_s        =  N_f_s_end ;  
     N_f_u        =  N_f_u_end ;  
     
     N_h          =  N_h_s + N_h_u;   %12
     N_f          =  N_f_s + N_f_u ;  %13
     N            =  N_h + N_f ;  %14

     mig          =  migbar;  %15
    % Choques 
     E            = 0 ;  %17
     A            = Abar;  %18
     Ak           = Akbar;
     miubar_s     = miu_ss_s;
     miubar_u     = miu_ss_u; 
%Capital 
     r_h          = (1/beta);  %19
     r_f          = r_h; %20

     V_h_s   = N_h_s/(((miu_h_s/miubar_s)^(1/eta))+(miu_h_s/x));
     V_h_u   = N_h_u/(((miu_h_u/miubar_u)^(1/eta))+(miu_h_u/x));
     V_f_s   = N_f_s/(((miu_f_s/miubar_s)^(1/eta))+(miu_f_s/x));
     V_f_u   = N_f_u/(((miu_f_u/miubar_u)^(1/eta))+(miu_f_u/x));

      v_h_s = V_h_s/N;
      v_h_u = V_h_u/N;
      v_f_s = V_f_s/N;
      v_f_u = V_f_u/N;    
     

      
     H_h_s  = (miu_h_s*V_h_s)/x;
     H_h_u  = (miu_h_u*V_h_u)/x;
     H_f_s  = (miu_f_s*V_f_s)/x;
     H_f_u  = (miu_f_u*V_f_u)/x;

     L_h_s        = H_h_s/N;
     L_h_u        = H_h_u/N;
     L_f_s        = H_f_s/N;
     L_f_u        = H_f_u/N;
     
%      H_h_s        = L_h_s*N;
%      H_h_u        = L_h_u*N ;
%      H_f_s        = L_f_s*N;
%      H_f_u        = L_f_u*N ;
%    
     H_u          = H_h_u + H_f_u;
     H_s          = H_h_s + H_f_s; 
     H_h          = H_h_s + H_h_u;
     H_f          = H_f_s + H_f_u;         
     
     L_s          = L_f_s +L_h_s ;
     L_u          = L_f_u +L_h_u  ;
     L_h          = L_h_s +L_h_u;
     L_f          = L_f_s +L_f_u  ;


     h_h_u        = H_h_u/N_h;
     h_h_s        = H_h_s/N_h;   
     h_f_u        = H_f_u/N_f;  
     h_f_s        = H_f_s/N_f;  
     h_s          = H_s/N;
     h_u          = H_u/N;
     h_h          = H_h/N_h;
     h_f          = H_f/N_f;

% Mercado de trabajo

     nmono_h_s    = N_h_s/N;
     nmono_h_u    = N_h_u/N;
     nmono_f_s    = N_f_s/N;
     nmono_f_u    = N_f_u/N;
     

     theta_h_s    =  v_h_s/(nmono_h_s-L_h_s); %11
     theta_h_u    =  v_h_u/(nmono_h_u-L_h_u); %12
     theta_f_s    =  v_f_s/(nmono_f_s-L_f_s); %13
     theta_f_u    =  v_f_u/(nmono_f_u-L_f_u); %14


     % Ponderadores
     phi_h        = N_h_s/ N_h;
     phi_f        = N_f_s/ N_f;
     omega        = N_h/N;
 
     
% Firmas     
     k            = (((1-rho)*(r_h-(1-delta))/(rho*(Ak^ipsilon)*pmg_s))^(1/(ipsilon-1)))*(L_s) ;
     i            = delta*k; %29
     y            = A*(sigma*((L_u)^alpha) + (1-sigma)*(rho*((Ak*k)^ipsilon)+(1-rho)*((L_s)^ipsilon))^(alpha/ipsilon))^(1/alpha) ; %28
  
   
% Variables percápita
  
      K     = k*N;
      I     = i*N;
      Y     = y*N;

    C_h = c_h*N_h;
    C_f = c_f*N_f;
    C = C_h + C_f;
    c = C/N;

    
 % Salarios  
     w_h_s = lambda_h_s*(pmg_s + (theta_h_s*q_h_s))+ (1-lambda_h_s)*((psi*((H_h/N_h)^gamma))*(C_h/N_h)*(1+tau_c)+b_h_s)/(1-tau_h_s);

 w_h_u = lambda_h_u*(pmg_u +(theta_h_u*q_h_u))+ (1-lambda_h_u)*((psi*((H_h/N_h)^gamma))*(C_h/N_h)*(1+tau_c)+b_h_u)/(1-tau_h_u);
 
 w_f_s = lambda_f_s*(pmg_s +(theta_f_s*q_f_s))+ (1-lambda_f_s)*((psi*((H_f/N_f)^gamma))*(C_f/N_f)*(1+tau_c)+b_f_s)/(1-tau_f_s);
 
  w_f_u = lambda_f_u*(pmg_u +(theta_f_u*q_f_u))+ (1-lambda_f_u)*((psi*((H_f/N_f)^gamma))*(C_f/N_f)*(1+tau_c)+b_f_u)/(1-tau_f_u);



 %Desempleo
    %Home
     U_h_s = N_h_s - H_h_s;
     U_h_u = N_h_u - H_h_u;
     U_h   = U_h_s + U_h_u;
     
    %Foreign 
     U_f_s = N_f_s - H_f_s;
     U_f_u = N_f_u - H_f_u;
     U_f   = U_f_s + U_f_u;
     
    %Skilled y Unskilled
    U_s = U_h_s + U_f_s;
    U_u = U_h_u + U_f_u;
    
    %Agregada
    U = U_h + U_f;
    
 % Tasas de desempleo
 
     %Home
     u_h_s = U_h_s/N_h_s;
     u_h_u = U_h_u/N_h_u;
     u_h = U_h/N_h;
     
    %Foreign 
     u_f_s = U_f_s/N_f_s;
     u_f_u = U_f_u/N_f_u;
     u_f = U_f/N_f;
     
    %Skilled y Unskilled
    u_s = U_s/(N_h_s + N_f_s);
    u_u = U_u/(N_h_u + N_f_u);
    
    %Agregada
    u = U/N;
    
     recaudo_c = tau_c*C;
     recaudo_h = tau_h_s*H_h_s*w_h_s + tau_h_u*H_h_u*w_h_u +tau_f_s*H_f_s*w_f_s +tau_f_u*H_f_u*w_f_u; 
     recaudo_c_y = recaudo_c/Y;
     recaudo_h_y = recaudo_h/Y;

     recaudo = recaudo_c + recaudo_h;
     gasto_flujo = 0;
     def = tau_c*C + tau_h_s*w_h_s*H_h_s + tau_h_u*w_h_u*H_h_u + tau_f_s*w_f_s*H_f_s + tau_f_u*w_f_u*H_f_u -  b_h_s*(N_h_s-H_h_s)- b_h_u*(N_h_u-H_h_u)- b_f_s*(N_f_s-H_f_s)- b_f_u*(N_f_u-H_f_u)-G;
     B = (def)/(r_h-1);
     B_h = phi_g*B;
     B_f = (1-phi_g)*B;
     pi = Y - w_h_s*H_h_s - w_h_u*H_h_u - w_f_s*H_f_s - w_f_u*H_f_u - V_h_s*q_h_s - V_h_u*q_h_u - V_f_s*q_f_s - V_f_u*q_f_u - I;
     pmg_k =((y^(1-alpha))*(A^alpha))*(1-sigma)*((rho*((Ak*k)^ipsilon) + (1-rho)*(L_s^ipsilon))^((alpha/ipsilon)-1))*rho*(Ak^ipsilon)*((k)^(ipsilon-1));        
     g_k = 1;
     
     deficit_y = def/Y;
      
     gasto_flujo_y = 0;
     Y_anual = Y*4;    
     def_anual = def*4;

     gasto_anual = gasto_flujo*4;
     gasto_flujo_y_anual = gasto_flujo/Y;
     d_y_anual = 0;
     d_c_anual = 0;
     d_i_anual = 0; 
     d_recaudo_anual = 0;
     ratio_s = V_f_s/V_h_s;
     ratio_u = V_f_u/V_h_u;
     ratio_f = V_f_s/V_f_u;
     ratio_h = V_h_s/V_h_u; 
     
     recaudo_h_mig = (tau_f_u*H_f_u*w_f_u +tau_f_s*H_f_s*w_f_s)/Y ;
     recaudo_h_loc = (tau_h_s*H_h_s*w_h_s + tau_h_u*H_h_u*w_h_u)/Y;
     recaudo_c_mig = (tau_c*C_f)/Y;
     recaudo_c_locales = (tau_c*C_h)/Y;



end;

resid;
check;
steady;

     shocks;
     
  var   e_mig=1000;

  end;
  

 

  

stoch_simul(order=1, irf=40, nograph);

