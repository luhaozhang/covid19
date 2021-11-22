% This script is used for plotting integrated figure of antiviral dynamics,
% fig.1 and fig.S1. Each time you run the script meaning you add one
% patient subplot to the total plot

clear;

n = 3; % The n^th patient/group. Each time you go to the next patient please set the current n to n+1

% Load data and best of the interested patient/group
load VTA_mild;
load VTA_mild_pbest_1211;
% Set index that chooses data to be fitted. v: virus; a: antibody; t: t cell
% Please copy the selected index from Table S5
% v_vector = [1,4,11,13,14,15,16:19,21,22]; t_vector = [1,5,7]; a_vector = [1:3];

% Best fit parameters. Please copy from Table S6
global T0;
global A0;

T0 = 0.89;
A0 = 11.53;

% Import parameters from pbest or other varible the reader set
alpha = pbest(1);
V0= 10^pbest(2);
N = pbest(3);
beta = pbest(4);
gamma = pbest(5);
% T
delta = pbest(6); 
epsilon = pbest(7);
% A
eta = pbest(8);
tau = pbest(9);
theta = pbest(10);

lags = tau;

options = ddeset('InitialY', [V0, 0, A0] );

% Time can be shorten or elongated
tspan = [N 80];

sol = dde23(@(t,y,Z) ddefun(t,y,Z,alpha,beta,delta,epsilon,gamma,eta,theta), lags, @history, tspan, options);

% Calculation of goodness of fit
% Please unmark the following section for individuals
        
%
% Please unmark the following section for mild and severe groups
% 
%         l1 = Deng_v( v_vector,1);
%         l2 = Deng_v( v_vector,2);
%         l3 = Lu_t(:,1);
%         l4 = Lu_t(:,2);
%         l5 = Deng_a( a_vector ,1);
%         l6 = Deng_a( a_vector ,2);

% Please unmark the following section for survivors and non-survivors
% 
%         l1 = surv_5p_v(v_vector,1);
%         l2 = surv_5p_v(v_vector,2);
%         l3 = surv_5p_t(t_vector,1);
%         l4 = surv_5p_t(t_vector,2) * 0.589;
%         l5 = surv_5p_a( a_vector ,1);
%         l6 = surv_5p_a( a_vector ,2);
% 
% % Calculate the corresponding simulated value from numerical result        
%         v_cal = deval(sol, l1);
%         v_cal = v_cal(1,:);
%         t_cal = deval(sol, l3);
%         t_cal = t_cal(2,:);
%         a_cal = deval(sol, l5);
%         a_cal = a_cal(3,:);

% Calculation of goodness of fit
%         v_square = 1 - 1/log10(max (sol.y(1,:))) * sqrt( sum ( ( log10(v_cal') - log10(l2) ).^2 ) / length(l1) );
%         t_square = 1 - 1/max (T0 - sol.y(2,:)) * sqrt( sum ( ( (T0 - t_cal)' - l4).^2 ) / length(l3) );
%         a_square = 1 - 1/max (sol.y(3,:)) * sqrt( sum ( ( a_cal' - l6).^2 ) / length(l5) );
        
% goodness = [v_square t_square a_square];
% disp('goodness = '); disp(goodness);


% Calculation of T cell accumulated amount of killed virus, tv;
for td=2:length(sol.x)
tv(td) = trapz(sol.x(1:td), sol.y(1,1:td).*sol.y(2,1:td)*beta );
end

% Calculation of antibody accumulated amount of killed virus, av;
for td=2:length(sol.x)
av(td) = trapz(sol.x(1:td), sol.y(1,1:td).*(sol.y(3,1:td)-A0)*gamma );
end

% Calculation of Ft
% Ft = tv(end)/(av(end) + tv(end)); disp(['Ft = ', num2str(Ft)]);

% Calculation of maximum of viral load

disp('max of antibody = ')
disp( (max (sol.y(3,:))) );

%% Plot Virus-T cell-Antibody along with bestfit
% Import data

% Please unmark the following section for individual patients. 
%         v1_total_day = p_5p_v(:,1);
%         v1_total_data = p_5p_v(:,2);
%         % For p1-p5, t1_total_data = p_5p_t(:,2) * 0.589;
%         t1_total_day = p_5p_t(:,1);
%         t1_total_data = p_5p_t(:,2);
% % For p1-p5, please mark "t1_neg" "t1_neg" lines.
%         t1_neg = p_5p_t(:,3); 
%         t1_pos = p_5p_t(:,4); 
%         a1_total_day = p_5p_a(:,1);
%         a1_total_data = p_5p_a(:,2);
        
% Please unmark the following section for mild or severe group and go to
% the function section in the end of the code to change the cofficient
% % after eta to 1e-7
        v1_total_day = Deng_v(:,1);
        v1_total_data = Deng_v(:,2);
        t1_total_day = Lu_t(:,1);
        t1_total_data = Lu_t(:,2);
        t1_neg = Lu_t(:,3);
        t1_pos = Lu_t(:,4);
        a1_total_day = Deng_a( : ,1);
        a1_total_data = Deng_a( : ,2);


fig_tlim = [-10 60];
t_vector = -20:1:tspan(2);
hold on;

% You can choose how to arrange all the patients' figure in one integrated plot
subplot(3,3,n)
% V
v_cal = sol.y(1,:);
v_cal(v_cal<=10)=1;
sol.y(1,:) = v_cal;
scatter(v1_total_day,v1_total_data,72, 'o', 'DisplayName','Virus','MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor',[1, 1, 1],'LineWidth',2);
hold on; plot( sol.x, sol.y(1,:), 'Linewidth', 2,'DisplayName','Fit','Color',[0,0.7,0.18]); % Virus
hold on;
% Limit of detection, for individuals it is 500 copy/mL, for mild and
% severe group it is 38000
hlod=plot( t_vector , 500 * ones(size(t_vector)),'--', 'Color','red','DisplayName','LOD', 'Linewidth', 2);%3.8E4 for Deng
xlim(fig_tlim);
yt = logspace(0,10,6);
set(gca,'YTick',yt, 'Yscale','log')
ylim([1e0,1e10])
%xlabel('Days after symptom onset');
%ylabel('Copy/mL');
set(gca,'YMinorTick','off','LineWidth',1.5)
set(gca,'FontSize',20)
box on;

% Antibody
subplot(3,3,3+n)
scatter(a1_total_day,a1_total_data, 84, 's','DisplayName','Virus','MarkerEdgeColor', [0, 0, 0],'MarkerFaceColor',[1, 1, 1],'LineWidth',2);
hold on; plot( sol.x, sol.y(3,:), '-', 'Linewidth', 2,'Color',[0,0.7,0.18]); % Activated innate immune cells
ylim([0 6]); 
xlim(fig_tlim);
set(gca,'YMinorTick','off','LineWidth',1.5)
set(gca,'FontSize',20)
box on;
% ylabel('OD'); or other units of the antibody


% T cell
subplot(3,3,6+n)
% Please mark the following line for p1-p5
errorbar(t1_total_day, T0 - t1_total_data, t1_neg, t1_pos,'s','MarkerSize',6,...
'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0],'Color',[0,0,0])
hold on;

scatter(t1_total_day, T0 - t1_total_data,72, '^', 'DisplayName','Virus','MarkerEdgeColor', [0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',2);

hold on; plot( sol.x, sol.y(2,:),  'Linewidth', 2,'Color',[0,0.7,0.18]); % Activated innate immune cells
% plot( t_vector ,0.69 * ones(size(t_vector)),'--', 'Color','black', 'Linewidth', 1); % Normal range for CD3+ from HongzhouLu's publication
% hold on;
% plot( t_vector ,2.54 * ones(size(t_vector)),'--', 'Color','black', 'Linewidth', 1); % Normal range for CD3+ from HongzhouLu's publication
% title('CD3+');
xlim(fig_tlim);
%xlabel('Days after symptom onset');
%ylabel('10^9/L');
set(gca,'FontSize',20)
ylim([0, 1])
set(gca,'YMinorTick','off','LineWidth',1.5)
box on;

%%

function dydt = ddefun(t,y,Z,alpha,beta,delta,epsilon,gamma,eta,theta)
global T0;
global A0;
  ylag1 = Z(:,1);
  dydt = [y(1)*(alpha - beta * y(2) - gamma*(y(3) - A0) )*(y(1)>10);
          delta *1e-8* y(1) - epsilon *0.1* y(2); 
          eta* 1e-7 *ylag1(1) - theta *0.01 * (y(3) - A0) ;]; % For mild and severe group, the coefficient after eta is 1e-7
end

function s = history(t)
global A0
    s(1) =0;
	s(2) = 0;
	s(3) = A0;
end
