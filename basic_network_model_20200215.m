clear all 
%% Set key parameters
global N_tot C sigma n m ed lambda w_l
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
m = 1.8; % Moving cost
ed = 5; % training cost 
lambda = 10; % objective function penalty for uncertainty 
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



%% Get next generation
% this function returns a new population of locations, skill levels, wages
% as well as showing how many in the last generation moved or chose to
% upgrade skills.   
[c_next,s_next,w_next,move,edu] = solveNextGen(d_l,d_s,c,s,w,w_bar); 

