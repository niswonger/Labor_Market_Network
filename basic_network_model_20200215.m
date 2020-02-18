clear all 
%% Set key parameters
global N_tot C sigma n m ed lambda w_l h eps
N_tot = 10000; % total population
C = 2; % number of cities
f_s = .3; % initial high skill fraction
w_l = 15; % low skilled wage
w_h1 = 20; % high skilled wage in city 1
w_h2 = 25; % high skilled wage in city 2
sigma = 2; % variance in wage
n = .05; % level of churn in the labor market 
d_l = 2*[3,1;
       1,3]; % mean connections for low skilled from city to city
d_s = 2*[10,3;
       3, 10]; % mean connections for high skilled from city to city
m = 1; % Moving cost
ed = 5; % training cost 
lambda = 10; % objective function penalty for uncertainty 
h = 2e-6; % housing cost
eps = 2; % housing elasticity
%% Set initial population
% Assign cities to the population
c = randsample(repmat(1:C,1,N_tot/C),N_tot)';
% Assign skills to the population
s = rand(N_tot,1)>(1-f_s);
% Set initial mean wages
w_bar = [ones(1,C)*w_l;w_h1,w_h2 ]';
% Define initial wage vector
w = zeros(N_tot,1);
% Instantiate index for city skills
cs_index = zeros(N_tot,C,2);
% Instantiate count of each group 
N = zeros(size(w_bar));
% get wages for the initial population
for i = 1:C 
    for j = 1:2
        cs_index(1:N_tot,i,j) = c == i & s == j-1;  
        N(i,j) = sum(cs_index(1:N_tot,i,j));
        w(logical(cs_index(1:N_tot,i,j))) ...
            = normrnd(w_bar(i,j),sigma,N(i,j),1);
    end
end



%% Play forwar for T periods
T = 100; 
pop = zeros(T,2);
skill = zeros(T,2);
movers = zeros(T,2);
% Prepare to loop for several generations
for i = 1:T
    pop(i,:) = [sum(c==1),sum(c==2)];
    skill(i,:) = [sum(c==1&s==1)/sum(c==1),sum(c==2&s==1)/sum(c==2)];
    % this function returns a new population of locations, skill levels, wages
    % as well as showing how many in the last generation moved or chose to
    % upgrade skills. 
    [c,s,w,move,edu] = ...
        solveNextGen(d_l,d_s,c,s,w,w_bar); 
    movers(i,:) = move;
end



plot(pop)
legend('1','2')
plot(movers)
legend('1','2')


figure
hold on 
plot(skill)
% Add labels
hXLabel = xlabel('Periods');
hYLabel = ylabel('Fraction Skilled');
% Add legend
hLegend = legend('High Return','Low Return');
% Adjust font
set([hXLabel, hYLabel, hLegend], 'FontSize', 15)
% set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w');

figure
plot(pop)
% Add labels
hXLabel = xlabel('Periods');
hYLabel = ylabel('Population');
% Add legend
hLegend = legend('High Return','Low Return');
% Adjust font
set([hXLabel, hYLabel, hLegend], 'FontSize', 15)
% set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w');

figure
plot(movers)
% Add labels
hXLabel = xlabel('Periods');
hYLabel = ylabel('Total Number of Movers');
% Add legend
hLegend = legend('High Return','Low Return');
% Adjust font
set([hXLabel, hYLabel, hLegend], 'FontSize', 15)
% set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'off', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w');
