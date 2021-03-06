% Least square fitting codes for IL-6, d-dimer and HSTC of survivors and
% non-survivors. The codes are used to perform bestfit and parameter
% uncertainty estimation
% Antiviral parameters are fixed. Please get these parameters from running
% "Surv_death_antiviral_bestfit.m".

% Load data and antiviral parameters of the interested group
load 5P_surv;
load surv_antiviral_pbest;

global T0 A0 I0 sd0 sh0 alpha V0 N beta gamma delta epsilon eta tau theta

% Please copy initial guess of parameters from Table S8, S7

T0= 0.84; 
A0= 0.13;
I0= 5.74;
sd0 = 0.25;
sh0 =0;

kappa = 1.8;
lambda = 2.32;
mu_d = 0.04;
nu_d = 0.04;
mu_h = 0.25;
nu_h = 0.42;

p0 = [kappa, lambda, mu_d, nu_d, mu_h, nu_h];

alpha = pbest(1); 
V0= 10^pbest(2);
N = pbest(3);% Incubation period
beta = pbest(4);
gamma = pbest(5);
% T
delta = pbest(6);
epsilon = pbest(7);
T02=T0;
% A
eta = pbest(8);
tau = pbest(9);
theta = pbest(10);

% Unmark the following section for survivors, and go to the last function
% section to unmark the corresponding section

% Import inflammation data to local variables
l1 = surv_5p_i(:,1);
l2 = surv_5p_i(:,2);
% n=1; % This line is used when estimating parameter uncertainty.  n=1 for
% IL-6, n=2 for d-dimer, n=3 for HSCT
% l2(n) = l2(n) - surv_5p_i(n,3); % This line is used when adding lower
% limit of data point for parameter uncertainty estimation. Same as below
% for all the marked lines with similar format
% l2(n) = l2(n) + surv_5p_i(n,4);
l3 = surv_5p_s_ddimer(:,1);
l4 = surv_5p_s_ddimer(:,2);
% l4(n) = l4(n) + surv_5p_s_ddimer(n,4);
l5 = surv_5p_s_hsct(:,1);
l6 = surv_5p_s_hsct(:,2);

% Unmark the following sections for non-survivors, and go to the last function
% section to unmark the corresponding section
% l1 = death_5p_i(:,1); 
% l2 = death_5p_i(:,2); 
% n = 4;
% % l2(n) = l2(n) - death_5p_i(n,3);
% % l2(n) = l2(n) + death_5p_i(n,4);
% l3 = death_5p_s_ddimer(:,1);
% l4 = death_5p_s_ddimer(:,2);
% % l4(n) = l4(n) - death_5p_s_ddimer(n,3);
% % l4(n) = l4(n) + death_5p_s_ddimer(n,4);
% l5 = death_5p_s_hsct(:,1);
% l6 = death_5p_s_hsct(:,2);
% % l6(n) = l6(n) - death_5p_s_hsct(n,3);
% l6(n) = l6(n) + death_5p_s_hsct(n,4);

% Fit

A = [];
b = [];
Aeq = [];
beq = [];
nonlon=[];

% Please copy upper and lower bound of parameter contraints to be varied
% from Table S8

lb = [ 0 0 0 0.4 0 0.4];
ub= [ 3 6 1 1 1 1];

options = optimoptions('fmincon','MaxFunctionEvaluations', 20000,'MaxIterations',2000);

% Fitting process

[infla_pbest, fval, exitflag, output, lambda, grad, hessian]  = fmincon(@(p) objective(p, l1, l2, l3, l4, l5, l6), p0, A,b,Aeq,beq, lb, ub,nonlon,options );

% Save for later use
save surv_infla_pbest infla_pbest

%% Objective function

function diff = objective(p, i_day, i_data, sd_day, sd_data, sh_day, sh_data)
global T0 A0 I0 sd0 sh0 alpha V0 N beta gamma delta epsilon eta tau theta

tspan = [N, 60];
lags = tau;
options = ddeset('InitialY', [V0, 0, A0, I0, sd0, sh0] );
sol = dde23(@(t,y,Z) ddefun(t,y,Z, alpha,beta,delta,epsilon,gamma,eta,theta,p(1),p(2),p(3),p(4),p(5),p(6)), lags, @history, tspan, options);

i_cal = deval(sol, i_day);
i_cal = i_cal(4,:);
sd_cal = deval(sol, sd_day);
sd_cal = sd_cal(5,:);
sh_cal = deval(sol, sh_day);
sh_cal = sh_cal(6,:);

% Please unmark the following section for non-survivors
% i_scale_factor = log10(max(i_data))^-1 ;
% i_diff = sum( (i_scale_factor* (log10(i_data) - log10(i_cal') ) ).^2 ) / length(i_data); 
% 
% sd_scale_factor = log10(max(sd_data))^-1 ;
% sd_diff = sum( (sd_scale_factor* (log10(sd_data) - log10(sd_cal') ) ).^2 ) / length(sd_data); 
% 
% sh_scale_factor = log10(max(sh_data))^-1 ;
% sh_diff = sum( (sh_scale_factor* (log10(sh_data) - log10(sh_cal') ) ).^2 ) / length(sh_data); 

% Please unmark the following section for survivors
 i_scale_factor = max(i_data)^-1;
 i_diff = sum( ( i_scale_factor* ( i_data - i_cal' ) ).^2 ) / length(i_data);
 
 sd_scale_factor = max(sd_data)^-1;
 sd_diff = sum( ( sd_scale_factor* ( sd_data - sd_cal' ) ).^2 ) / length(sd_data);
 
 sh_scale_factor = max(sh_data)^-1;
 sh_diff = sum( ( sh_scale_factor* ( sh_data - sh_cal' ) ).^2 ) / length(sh_data);
 
 diff = i_diff + sh_diff + sd_diff ;

end

function dydt = ddefun(t,y,Z,alpha,beta,delta,epsilon,gamma,eta,theta,kappa,lambda,mu1,nu1,mu2,nu2)
global A0;
  ylag1 = Z(:,1);
  dydt = [y(1)*(alpha - beta * y(2) - gamma*( y(3) - A0) )*(y(1)>10); 
          delta *1e-8* y(1) - epsilon *0.1 * y(2); 
          eta*1e-8*ylag1(1) - theta * 0.01 * (y(3) -A0 ) ; 
          kappa * (y(3) -A0 )  - lambda * ( y(4) - 5.4  );  % 5.4 is the upper limit of normal range
          mu1 * y(4) - nu1 * y(5);
          mu2 * y(4) - nu2 * y(6)];
end

function s = history(t)
global A0
global I0
global sd0;
global sh0;
global T0;

    s(1) =0;
	s(2) = T0;
	s(3) = A0;
    s(4) = I0;
    s(5) = sd0;
    s(6) = sh0;
end


