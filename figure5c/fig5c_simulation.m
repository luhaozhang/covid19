% Plot fits of survivors and non-survivors as fig.3
% Virus, T cell, Antibody, IL-6, D-dimer, HSCT, V, T, A, I, S

global T0;
global A0; 
global I0; 
global sd0; 
global sh0; 

% Non-survivors:death people

for count = 1:8

T0 = 0.5*count;
A0 = 0.08;
I0 = 10.73;
sd0 = 0.25;
sh0 = 0;

load death_antiviral_pbest
load death_infla_pbest

% % equation parameters
alpha = pbest(1); 
V0= 10^pbest(2);
N = pbest(3);% Incubation period
beta = pbest(4);
gamma = pbest(5);
% T
delta = 0.25+21.67 * (T0-0.39);
epsilon = pbest(7);
T02=T0;
% A
eta = pbest(8);
tau = pbest(9);
theta = pbest(10);
% I
kappa = infla_pbest(1);
lambda = 4.62 * (T0 - 0.39);
% D-dimer
mu1 = infla_pbest(3);
nu1 = 0.889 * (T0 - 0.39);
% HSCT
mu2 = infla_pbest(5);
nu2 = 0.933 * (T0 - 0.39);

% non-survivors: simulation

tspan = [N, 50];
lags = tau;
options = ddeset('InitialY', [V0, 0, A0, I0, sd0, sh0] );
sol2 = dde23(@(t,y,Z) ddefun(t,y,Z, alpha, beta, delta, epsilon, gamma, eta, theta, kappa, lambda, mu1, nu1, mu2, nu2), lags, @history, tspan, options);

%% Plot

% fig_tlim = [-10 tspan(2)];
% t_vector = -10:1:tspan(2);
% h =figure;
% set(h,'DefaultAxesFontSize',22)
% % V
% subplot(2,3,1);
% v_cal = sol1.y(1,:);
% v_cal(v_cal<=10)=1;
% sol1.y(1,:) = v_cal;
% scatter(v1_total_day,v1_total_data,84, '^', 'DisplayName','Virus','MarkerEdgeColor', [0.75,0.00,0.75], 'MarkerFaceColor',[0.75,0.00,0.75]);
% hold on; plot( sol1.x, sol1.y(1,:), 'Linewidth', 3,'DisplayName','Fit','Color', [0.75,0.00,0.75]); % Virus
% hold on;
% scatter(v2_total_day,v2_total_data,84, 's', 'DisplayName','Virus','MarkerEdgeColor', [0 0 0],'MarkerFaceColor',[0 0 0]);
% v_cal = sol2.y(1,:);
% v_cal(v_cal<=10)=1;
% sol2.y(1,:) = v_cal;
% hold on; plot( sol2.x, sol2.y(1,:), 'Linewidth', 3,'DisplayName','Fit','Color', [0 0 0]); % Virus
% 
% title('Viral Load');
% hold on; 
% hlod=plot( t_vector ,100 * ones(size(t_vector)),'--', 'Color',[0.75, 0, 0.75],'DisplayName','LOD', 'Linewidth', 1); % LOD kelvin: 500, Francois��100 . take 100 for mixed data
% xlim(fig_tlim);
% yt = logspace(1,9,5);
% set(gca,'YTick',yt, 'Yscale','log')
% ylim([1e1,1e9])
% % xlabel('Days after symptom onset');
% ylabel('Copy/mL');
% set(gca,'YMinorTick','off','LineWidth',1.5)
% box on;
% 
% 
% % T
% subplot(2,3,2);
% errorbar(t1_total_day, t1_total_data, t1_neg, t1_pos,'^','MarkerSize',10,...
% 'MarkerEdgeColor',[0.75,0.00,0.75],'MarkerFaceColor',[0.75,0.00,0.75],'Color',[0.75,0.00,0.75])
% hold on; plot( sol1.x, T01 * ones(size(sol1.x)) - sol1.y(2,:),  'Linewidth', 3,'Color',[0.75,0.00,0.75]); % Activated innate immune cells
% hold on;
% errorbar(t2_total_day, t2_total_data, t2_neg, t2_pos,'s','MarkerSize',10,...
% 'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],'Color',[0,0,0])
% hold on; plot( sol2.x, T02 * ones(size(sol2.x)) - sol2.y(2,:),  'Linewidth', 3,'Color',[0,0,0]); % Activated innate immune cells
% hold on;
% plot( t_vector , 0.69 * ones(size(t_vector)),'--', 'Color',[0.25, 0.25, 0.25], 'Linewidth', 1);
% hold on;
% plot( t_vector , 2.54 * ones(size(t_vector)),'--', 'Color',[0.25, 0.25, 0.25], 'Linewidth', 1);
% title('T cell');
% xlim(fig_tlim);
% % xlabel('Days after symptom onset');
% ylabel('10^9/L');
% ylim([0, 2])
% set(gca,'YMinorTick','off','LineWidth',1.5)
% box on;
% 
% % A
% subplot(2,3,3);
% scatter(a1_total_day,a1_total_data, 84, '^','DisplayName','Virus','MarkerEdgeColor', [0.75,0.00,0.75],'MarkerFaceColor',[0.75,0.00,0.75]);
% hold on; plot( sol1.x, sol1.y(3,:),  'Linewidth', 3,'Color',[0.75,0.00,0.75]); % Activated innate immune cells
% hold on;
% scatter(a2_total_day,a2_total_data, 84, 's','DisplayName','Virus','MarkerEdgeColor', [0,0,0],'MarkerFaceColor',[0,0,0]);
% hold on; plot( sol2.x, sol2.y(3,:),  'Linewidth', 3,'Color',[0,0,0]); % Activated innate immune cells
% title('IgG');
% xlim(fig_tlim);
% % xlabel('Days after symptom onset');
% ylabel('OD or other');
% ylim([0, 6])
% set(gca,'YMinorTick','off','LineWidth',1.5)
% box on;
% 
% % I
% subplot(2,3,4);
% errorbar(i1_total_day, i1_total_data, i1_neg, i1_pos,'^','MarkerSize',10,...
% 'MarkerEdgeColor',[0.75,0.00,0.75],'MarkerFaceColor',[0.75,0.00,0.75],'Color', [0.75,0.00,0.75])
% hold on; plot( sol1.x, sol1.y(4,:),  'Linewidth', 2, 'DisplayName','Fit','Color',[0.75,0.00,0.75], 'Linewidth', 3); 
% hold on;
% errorbar(i2_total_day, i2_total_data, i2_neg, i2_pos,'s','MarkerSize',10,...
% 'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],'Color', [0,0,0])
% hold on; plot( sol2.x, sol2.y(4,:),  'Linewidth', 2, 'DisplayName','Fit','Color',[0,0,0], 'Linewidth', 3); 
% % normal range
% hold on;
% plot( t_vector ,5.4 * ones(size(t_vector)),'--', 'Color',[0.25, 0.25, 0.25], 'Linewidth', 1);   % Lu's extended data
% title('IL-6');
% % title('Fever');
% xlim(fig_tlim);
% % xlabel('Days after symptom onset');
% ylabel('pg/mL')
% ylim([0 60]);
% set(gca,'Yscale','log');
% set(gca,'YMinorTick','off','LineWidth',1.5)
% box on;

%% D-dimer
hold on; plot( sol2.x, sol2.y(5,:), 'Linewidth', 3,'DisplayName',['{\itT}_0 = ',num2str(T0)],'Color',[1 - 0.1 * count,0.1 * count,0.1* count]);
% normal range
hold on;
% plot( t_vector ,0.5 * ones(size(t_vector)),'--', 'Color',[0.25, 0.25, 0.25], 'Linewidth', 1);  % Lu's extended data
title('D-dimer');
xlim([0 22]);
set(gca,'Yscale','log')
ylim([0.1 100]);
xlabel('Days after symptom onset');
ylabel('\mug/L');
set(gca,'YMinorTick','off','LineWidth',1.5)
box on;
legend('box','off');

% % HSCT
% subplot(2,3,6);
% % Unmark the following line for survivors
% scatter(s1_h_total_day,s1_h_total_data, 84, '^','DisplayName','Virus','MarkerEdgeColor', [0.75,0.00,0.75],'MarkerFaceColor',[0.75,0.00,0.75]);
% hold on; 
% plot( sol1.x, sol1.y(6,:), 'Linewidth', 3,'DisplayName','Fit','Color',[0.75,0.00,0.75]); % Virus
% hold on;
% errorbar(s2_h_total_day, s2_h_total_data, s2_h_neg, s2_h_pos,'s','MarkerSize',10,...
% 'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],'Color', [0,0,0])
% hold on; plot( sol2.x, sol2.y(6,:), 'Linewidth', 3,'DisplayName','Fit','Color',[0,0,0]); 
% hold on;
% plot( t_vector , 25 * ones(size(t_vector)),'--', 'Color',[0.25, 0.25, 0.25], 'Linewidth', 1); % Average of women and men from ﻿David A. McAllister's paper
% title('HSCT');
% xlim(fig_tlim);
% set(gca,'Yscale','log')
% ylim([0.1 10000]);
% % xlabel('Days after symptom onset');
% ylabel('ng/mL');
% set(gca,'YMinorTick','off','LineWidth',1.5)
% box on;

set(findall(gcf,'-property','FontSize'),'FontSize',24);

end

function dydt = ddefun(t,y,Z,alpha,beta,delta,epsilon,gamma,eta,theta,kappa,lambda,mu1,nu1,mu2,nu2)
global A0;
global I0;
  ylag1 = Z(:,1);
  dydt = [y(1)*(alpha - beta * y(2) - gamma*( y(3) - A0) )*(y(1)>10); 
          delta *1e-8* y(1) - epsilon *0.1 *y(2); 
          eta*1e-8*ylag1(1) - theta *0.01*(y(3) -A0 ) ; 
          kappa * (y(3) -A0 )  - lambda * ( y(4) - 5.4  ); 
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


