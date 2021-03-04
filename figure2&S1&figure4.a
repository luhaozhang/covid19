% Best fit of Virus, T cell, Antibody data
% With the codes you should be able to reproduce bestfit and parameter
% uncertainty in Table S4. Applied to individual patients.
clear

% Data import. Please select the mat file of the intersted patient.
load no902_3o;
% Set index that chooses data to be fitted. v: virus; a: antibody; t: t cell
% Please copy the selected index from Table S3
v_vector = [1:7]; t_vector =  1:16; a_vector = 1:7;
        l1 = p_5p_v(v_vector,1);
        l2 = p_5p_v(v_vector,2);
        l3 = p_5p_t(t_vector,1);
        l4 = p_5p_t(t_vector,2); % If the patients are p1/p3/p4/p5, then l4 = p_5p_t(t_vector,2) * 0.589;
        l5 = p_5p_a( a_vector ,1);
        l6 = p_5p_a( a_vector ,2);

% Initial guess. Please copy initial guess from TableS7
alpha = 1;
V0 = 1.5;
beta = 4.5;
gamma = 1.12;
delta = 0.5;
epsilon = 1.4;
eta = 1.8;
tau=4.5;
theta = 5.11;

global T0;
global A0; 
T0 = 1.06;
A0 = 0.11;
N = -10;

p0 = [alpha, V0, N, beta, gamma, delta, epsilon, eta, tau, theta];

% Fit

A = [];
b = [];
Aeq = [];
beq = [];
nonlon=[];

% Upper bound and lower bound of the parameter constraints to be optimized.
% Please copy them from Table S6
lb = [0.1, 1, -10, 1, 0.1, 0.01, 1.4, 0.02, 1, 5.11];
ub= [2.5, 8, -10, 10, 5, 5, 1.4, 5, 20, 5.11];

options = optimoptions('fmincon','MaxFunctionEvaluations', 20000,'MaxIterations',2000);

[pbest, fval, exitflag, output, lambda, grad, hessian]  = fmincon(@(p) objective(p, l1, l2, l3, l4, l5, l6), p0, A,b,Aeq,beq, lb, ub,nonlon,options );

% Save the fit either for figure production of bestfit or uncertainty
% estimation.
save no902_3o_pbest  pbest

%% Objective function

function diff = objective(p, v_day, v_data, t_day, t_data,a_day, a_data)
global A0
global T0
tspan = [p(3), 60];
lags = p(9);
options = ddeset('InitialY', [10^p(2), 0, A0] );
sol = dde23(@(t,y,Z) ddefun(t,y,Z, p(1), p(4), p(5), p(6), p(7),p(8),p(10)), lags, @history, tspan, options);
v_cal = deval(sol, v_day);
v_cal = v_cal(1,:);
t_cal = deval(sol, t_day);
t_cal = t_cal(2,:);
a_cal = deval(sol, a_day);
a_cal = a_cal(3,:);

v_scale_factor = log10(max(v_data))^-1 ;
v_diff = sum( (v_scale_factor* (log10(v_data) - log10(v_cal') ) ).^2 ) / length(v_data); 
t_scale_factor = max(T0 - t_data)^-1;
t_diff = sum( ( t_scale_factor * ( T0 - t_data - t_cal' ) ).^2 )/length(t_data); 
 a_scale_factor = max(a_data)^-1;
a_diff = sum( ( a_scale_factor* ( a_data - a_cal' ) ).^2 ) / length(a_data);

diff = v_diff + t_diff + a_diff;

end

function dydt = ddefun(t,y,Z,alpha, beta, gamma, delta, epsilon, eta, theta)
global A0;
  ylag1 = Z(:,1);
  dydt = [y(1)*(alpha - beta * y(2) - gamma* ( y(3) - A0) )*(y(1)>10); 
          delta *1e-8* y(1) - epsilon *0.1 * y(2); 
          eta*1e-8*ylag1(1) - theta * 0.01* (y(3) -A0 )];
end

function s = history(t)
global A0
    s(1) =0;
	s(2) = 0;
	s(3) = A0;
end


