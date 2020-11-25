function f = solve_maduro_endval_alto(xe, q_h_s_cali, q_h_u_cali, q_f_s_cali, q_f_u_cali)

load q_h_s_cali;
load q_h_u_cali;
load q_f_s_cali;
load q_f_u_cali;


q_h_s = q_h_s_cali;
q_h_u = q_h_u_cali;
q_f_s = q_f_s_cali;
q_f_u = q_f_u_cali;



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
      kappa = 1;
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
      
      rhopob=0.5;
      
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

      b_flujo = 0.31;
      phi_g   = 0.9;   
      G=800;
 % % ------------------ % % 
% % Solución del solver % % 
%_________________________

    pmg_s       = xe(1);
    pmg_u       = xe(2);
    c_h         = xe(3);
    c_f         = xe(4);
    miu_h_s     = xe(5);
    miu_h_u     = xe(6);
    miu_f_s     = xe(7);
    miu_f_u     = xe(8);
        
    
    
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

    % Choques 
     mig          =  migbar;  %15
     E            = 0 ;    %17
     A            = Abar;  %18
     Ak           = Akbar;
     miubar_s     = miu_ss_s;
     miubar_u     = miu_ss_u; 
%Capital 
     r_h          = (1/beta); %19
     r_f          = r_h; %20
     
% Vacantes
%      v_h_s        = (x*L_h_s/(miubar_s*((nmono_h_s-L_h_s)^eta)))^(1/(1-eta)); 
%      v_h_u        = (x*L_h_u/(miubar_u*((nmono_h_u-L_h_u)^eta)))^(1/(1-eta));
%      v_f_s        = (x*L_f_s/(miubar_s*((nmono_f_s-L_f_s)^eta)))^(1/(1-eta)); 
%      v_f_u        = (x*L_f_u/(miubar_u*((nmono_f_u-L_f_u)^eta)))^(1/(1-eta));

    V_h_s   = N_h_s/(((miu_h_s/miubar_s)^(1/eta))+(miu_h_s/x));
     V_h_u   = N_h_u/(((miu_h_u/miubar_u)^(1/eta))+(miu_h_u/x));
     V_f_s   = N_f_s/(((miu_f_s/miubar_s)^(1/eta))+(miu_f_s/x));
     V_f_u   = N_f_u/(((miu_f_u/miubar_u)^(1/eta))+(miu_f_u/x));

      v_h_s = V_h_s/N;
      v_h_u = V_h_u/N;
      v_f_s = V_f_s/N;
      v_f_u = V_f_u/N;    
       
     
%      V_h_s = v_h_s*N;
%       V_h_u = v_h_u*N;
%       V_f_s = v_f_s*N;
%       V_f_u = v_f_u*N;    
%      
     H_h_s  = (miu_h_s*V_h_s)/x;
     H_h_u  = (miu_h_u*V_h_u)/x;
     H_f_s  = (miu_f_s*V_f_s)/x;
     H_f_u  = (miu_f_u*V_f_u)/x;
     
     L_h_s        = H_h_s/N;
     L_h_u        = H_h_u/N;
     L_f_s        = H_f_s/N;
     L_f_u        = H_f_u/N;
   
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
     


     theta_h_s    =  v_h_s/(nmono_h_s-L_h_s); 
     theta_h_u    =  v_h_u/(nmono_h_u-L_h_u); 
     theta_f_s    =  v_f_s/(nmono_f_s-L_f_s); 
     theta_f_u    =  v_f_u/(nmono_f_u-L_f_u); 
  
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
 
 %   w_h_s = lambda_h_s*(pmg_s + (theta_h_s*q_h_s))+ (1-lambda_h_s)*((psi*(h_h^gamma))*c_h*(1+tau_c)+b_h_s)/(1-tau_h_s);
     w_h_s = lambda_h_s*(pmg_s + (theta_h_s*q_h_s))+ (1-lambda_h_s)*((psi*((H_h/N_h)^gamma))*(C_h/N_h)*(1+tau_c)+b_h_s)/(1-tau_h_s);

  %  w_h_u = lambda_h_u*(pmg_u + (theta_h_u*q_h_u))+ (1-lambda_h_u)*((psi*(h_h^gamma))*c_h*(1+tau_c)+b_h_u)/(1-tau_h_u);
     w_h_u = lambda_h_u*(pmg_u +(theta_h_u*q_h_u))+ (1-lambda_h_u)*((psi*((H_h/N_h)^gamma))*(C_h/N_h)*(1+tau_c)+b_h_u)/(1-tau_h_u);
 
  %  w_f_s = lambda_f_s*(pmg_s + (theta_f_s*q_f_s))+ (1-lambda_f_s)*((psi*(h_f^gamma))*c_f*(1+tau_c)+b_f_s)/(1-tau_f_s);
     w_f_s = lambda_f_s*(pmg_s +(theta_f_s*q_f_s))+ (1-lambda_f_s)*((psi*((H_f/N_f)^gamma))*(C_f/N_f)*(1+tau_c)+b_f_s)/(1-tau_f_s);
 
  %  w_f_u = lambda_f_u*(pmg_u + (theta_f_u*q_f_u))+ (1-lambda_f_u)*((psi*(h_f^gamma))*c_f*(1+tau_c)+b_f_u)/(1-tau_f_u);
     w_f_u = lambda_f_u*(pmg_u +(theta_f_u*q_f_u))+ (1-lambda_f_u)*((psi*((H_f/N_f)^gamma))*(C_f/N_f)*(1+tau_c)+b_f_u)/(1-tau_f_u);


      u= 1-(H_h_s+H_h_u+H_f_s+H_f_u)/N;
      u_h = 1- h_h;
      u_f = 1- h_f;
      u_h_s = 1- h_h_s;
      u_h_u = 1- h_h_u;      
      u_f_s = 1- h_f_s;      
      u_f_u = 1- h_f_u;    
      

      
     recaudo_c = tau_c*C;
     recaudo_h = tau_h_s*H_h_s*w_h_s + tau_h_u*H_h_u*w_h_u +tau_f_s*H_f_s*w_f_s +tau_f_u*H_f_u*w_f_u; 
     recaudo = recaudo_c + recaudo_h;
     gasto_flujo = 0;
     def = tau_c*C + tau_h_s*w_h_s*H_h_s + tau_h_u*w_h_u*H_h_u + tau_f_s*w_f_s*H_f_s + tau_f_u*w_f_u*H_f_u -  b_h_s*(N_h_s-H_h_s)- b_h_u*(N_h_u-H_h_u)- b_f_s*(N_f_s-H_f_s)- b_f_u*(N_f_u-H_f_u)-G;
     B = (def)/(r_h-1);
     B_h = phi_g*B;
     B_f = (1-phi_g)*B;
     pi = Y - w_h_s*H_h_s - w_h_u*H_h_u - w_f_s*H_f_s - w_f_u*H_f_u - V_h_s*q_h_s - V_h_u*q_h_u - V_f_s*q_f_s - V_f_u*q_f_u - I;
     pmg_k =((y^(1-alpha))*A^alpha)*(1-sigma)*((rho*((Ak*k)^ipsilon) + (1-rho)*(L_s^ipsilon))^((alpha/ipsilon)-1))*rho*(Ak^ipsilon)*((k)^(ipsilon-1));        
     g_k = 1;
     
     
      f = [];

          f(1) =  pmg_s - (((y^(1-alpha))*(A^alpha))*((1-sigma)*((rho*((Ak*k)^ipsilon)+(1-rho)*(L_s^ipsilon))^((alpha/ipsilon)-1))*(1-rho)*(L_s^(ipsilon-1)))); 
          f(2) =  pmg_u - (((y^(1-alpha))*(A^alpha))*sigma*(L_u^(alpha -1)));            
          f(3) =  c_h   - ((w_h_s*h_h_s*(1-tau_h_s) + w_h_u*h_h_u*(1-tau_h_u) + b_h_s*(phi_h-h_h_s) +  b_h_u*((1-phi_h)-h_h_u) + (B_h/N_h)*(r_h-1) + (pi/N_h))/(1+tau_c)); 
          f(4) =  c_f   - ((w_f_s*h_f_s*(1-tau_f_s) + w_f_u*h_f_u*(1-tau_f_u) + b_f_s*(phi_f-h_f_s) +  b_f_u*((1-phi_f)-h_f_u) + (B_f/N_f)*(r_h-1))/(1+tau_c));
          f(5) =  miu_h_s - (r_h*(q_h_s*(1-beta*(1-x)))/((pmg_s - w_h_s)));  
          f(6) =  miu_h_u - (r_h*(q_h_u*(1-beta*(1-x)))/((pmg_u - w_h_u))); 
          f(7) =  miu_f_s - (r_f*(q_f_s*(1-beta*(1-x)))/((pmg_s - w_f_s))); 
          f(8) =  miu_f_u - (r_f*(q_f_u*(1-beta*(1-x)))/((pmg_u - w_f_u))); 

          f = f';  

end