
init= 5*ones(4,1);
[xs,fval] = fsolve(@solve_maduro_init,init); % Call solver

pmg_s = (xs(1));
pmg_u = (xs(2));
c_h   = (xs(3));
c_f   = (xs(4));
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
      g_f_u_bar =1;
      g_h_bar = 1;
      
%       rhopob=0.5;
      
%      N_h_s_end   =  10554.128;   
%      N_h_u_end   =  36661.677;  
%      N_f_s_end   =  152.717; 
%      N_f_u_end   =  1022.026 + 1143.962;  

     N_h_s_bar   =  10554.128;   
     N_h_u_bar   =  36661.677;  
     N_f_s_bar   =  152.717; 
     N_f_u_bar   =  1022.026; 

    % Migración
%       rhomig = 0.01;
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


      phi_g   = 0.9;   
        G=800;
%__________________________      
% Tasas de crecimiento
     g            = 1; %1
     g_h          = g_h_bar; %2
     g_f          = 1; %3
     g_h_s        = g_h_s_bar; %4
     g_h_u        = g_h_u_bar; %5
     g_f_s        = g_f_s_bar; %6
     g_f_u        = 1; %7

% Población (Exógena)
     N_h_s        =  N_h_s_bar ;   
     N_h_u        =  N_h_u_bar ;  
     N_f_s        =  N_f_s_bar ;  
     N_f_u        =  N_f_u_bar ;  
 
     N_h          =  N_h_s + N_h_u;   %12
     N_f          =  N_f_s + N_f_u ;  %13
     N            =  N_h + N_f ;  %14

    C_h = c_h*N_h;
    C_f = c_f*N_f;
    C = C_h + C_f;
      
     mig          =  migbar;  %15
    % Choques 
     E            = 0 ;  %17
     A            = Abar;  %18
     Ak           = Akbar;
     miubar_s     = miu_ss_s;
     miubar_u     =miu_ss_u; 
%Capital 
     r_h          = (1/beta) ; %19
     r_f          = r_h; %20
     
% Empleo (TARGET)
    
     L_h_s        = 0.1981 ;
     L_h_u        = 0.6883 ;
     L_f_s        = 0.0026 ;
     L_f_u        = 0.0173 ;
     
     H_h_s        = L_h_s*N;
     H_h_u        = L_h_u*N ;
     H_f_s        = L_f_s*N;
     H_f_u        = L_f_u*N ;
   
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
     
     v_h_s        = (x*L_h_s/(miubar_s*((nmono_h_s-L_h_s)^eta)))^(1/(1-eta)); 
     v_h_u        = (x*L_h_u/(miubar_u*((nmono_h_u-L_h_u)^eta)))^(1/(1-eta));
     v_f_s        = (x*L_f_s/(miubar_s*((nmono_f_s-L_f_s)^eta)))^(1/(1-eta)); 
     v_f_u        = (x*L_f_u/(miubar_u*((nmono_f_u-L_f_u)^eta)))^(1/(1-eta));

     theta_h_s    =  v_h_s/(nmono_h_s-L_h_s); %11
     theta_h_u    =  v_h_u/(nmono_h_u-L_h_u); %12
     theta_f_s    =  v_f_s/(nmono_f_s-L_f_s); %13
     theta_f_u    =  v_f_u/(nmono_f_u-L_f_u); %14

     miu_h_s      = miubar_s*(theta_h_s^(-eta)); %19
     miu_h_u      = miubar_u*(theta_h_u^(-eta));%20
     miu_f_s      = miubar_s*(theta_f_s^(-eta)); %21
     miu_f_u      = miubar_u*(theta_f_u^(-eta)); %22     
  
     % Ponderadores
     phi_h        = N_h_s/ N_h;
     phi_f        = N_f_s/ N_f;
     omega        = N_h/N;
 
     
% Firmas     
     k            = (((1-rho)*(r_h-(1-delta))/(rho*(Ak^ipsilon)*pmg_s))^(1/(ipsilon-1)))*(L_s) ;
     i            = delta*k; %29
     y            = A*(sigma*((L_u)^alpha) + (1-sigma)*(rho*((Ak*k)^ipsilon)+(1-rho)*((L_s)^ipsilon))^(alpha/ipsilon))^(1/alpha) ; %28
   
      V_h_s = v_h_s*N;
      V_h_u = v_h_u*N;
      V_f_s = v_f_s*N;
      V_f_u = v_f_u*N;    
   
% Variables percápita

      K     = k*N;
      I     = i*N;
      Y     = y*N;

      q_h_s = (miu_h_s*beta*(1-lambda_h_s)*(pmg_s-((psi*(h_h^gamma)*c_h*(1+tau_c) + b_h_s)/(1-tau_h_s))))/(1+beta*theta_h_s*lambda_h_s*miu_h_s -beta*(1-x));
      q_h_u = (miu_h_u*beta*(1-lambda_h_u)*(pmg_u-((psi*(h_h^gamma)*c_h*(1+tau_c) + b_h_u)/(1-tau_h_u))))/(1+beta*theta_h_u*lambda_h_u*miu_h_u -beta*(1-x));
      q_f_s = (miu_f_s*beta*(1-lambda_f_s)*(pmg_s-((psi*(h_f^gamma)*c_f*(1+tau_c) + b_f_s)/(1-tau_f_s))))/(1+beta*theta_f_s*lambda_f_s*miu_f_s -beta*(1-x));
      q_f_u = (miu_f_u*beta*(1-lambda_f_u)*(pmg_u-((psi*(h_f^gamma)*c_f*(1+tau_c) + b_f_u)/(1-tau_f_u))))/(1+beta*theta_f_u*lambda_f_u*miu_f_u -beta*(1-x));
      
    C_h = c_h*N_h;
    C_f = c_f*N_f;
    C = C_h + C_f;
    c = c_f*(1-omega) + c_h*omega;

 
 % Salarios  
 
    w_h_s = lambda_h_s*(pmg_s + (theta_h_s*q_h_s))+ (1-lambda_h_s)*((psi*(h_h^gamma))*c_h*(1+tau_c)+b_h_s)/(1-tau_h_s);
    w_h_u = lambda_h_u*(pmg_u + (theta_h_u*q_h_u))+ (1-lambda_h_u)*((psi*(h_h^gamma))*c_h*(1+tau_c)+b_h_u)/(1-tau_h_u);
    w_f_s = lambda_f_s*(pmg_s + (theta_f_s*q_f_s))+ (1-lambda_f_s)*((psi*(h_f^gamma))*c_f*(1+tau_c)+b_f_s)/(1-tau_f_s);
    w_f_u = lambda_f_u*(pmg_u + (theta_f_u*q_f_u))+ (1-lambda_f_u)*((psi*(h_f^gamma))*c_f*(1+tau_c)+b_f_u)/(1-tau_f_u);

    
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
   recaudo = recaudo_c + recaudo_h;
   gasto_flujo = 0;
   def = tau_c*C + tau_h_s*w_h_s*H_h_s + tau_h_u*w_h_u*H_h_u + tau_f_s*w_f_s*H_f_s + tau_f_u*w_f_u*H_f_u -  b_h_s*(N_h_s-H_h_s)- b_h_u*(N_h_u-H_h_u)- b_f_s*(N_f_s-H_f_s)- b_f_u*(N_f_u-H_f_u)-G;
  
     B = (def)/(r_h-1);
     B_h = phi_g*B;
     B_f = (1-phi_g)*B;
 

q_h_s_cali = q_h_s;
q_h_u_cali = q_h_u;
q_f_s_cali = q_f_s;
q_f_u_cali = q_f_u;



    
save q_h_s_cali;
save q_h_u_cali;
save q_f_s_cali;
save q_f_u_cali;
save Y_init;


