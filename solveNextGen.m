% This function takes in the degree distribution means for poisson random
% graph as well as the current skill levels, locations and wages and
% generates a new population where individuals decide where to live and
% whether to upgrade skills based on sempling wages from their network

function [c_next,s_next,w_next,move,edu] = solveNextGen(d_l,d_s,c,s,w,w_bar)
    global N_tot C sigma n m ed lambda w_l h eps
    %% Generate "new" people
    % Get proportional new labor 
    N_new = N_tot*n;
    % Define mean connections for each skill level going from city i to city j
    d = cat(3,d_l,d_s);
    % Store mean observed wages 
    w_obs = zeros(N_new,C);
    d_obs = zeros(N_new,C);
    c_new = zeros(N_new,1);

    skill_new = zeros(N_new,1);
    for j = 1:C 
        index = 1;
        for i = 1:C
            for k = 1:2
                index_new = round(n*sum(c == i & s == k-1)); % New population is proportional
                % get the number of connections from a poisson draw
                d_samp = poissrnd(d(i,j,k),index_new,1);
                % Sample from the relevant wage to get the mean from the sample 
                fun = @(d)mean(randsample(w(c == j & s == 1),min(d,sum(c == j & s == 1))));
                w_obs(index:index+index_new-1,j) = arrayfun(fun,d_samp);
                c_new(index:index+index_new-1) = i;
                d_obs(index:index+index_new-1,j) = d_samp;
                skill_new(index:index+index_new-1,1) = k;
                index = index + index_new;
            end
        end
    end
    w_obs(isnan(w_obs))=0;
    % New people choose skills and location
    out = w_l; %-h*[sum(c==1),sum(c==2)].^(eps); % Outside option of no upgrading
    U_new = w_obs-lambda./(d_obs).^.5-[m,0].*(c_new==2)-[0,m].*(c_new==1)-ed; % -h*[sum(c==1),sum(c==2)].^(eps);
    m_obs_1 = (c_new == 2 & (U_new(:,1) >  U_new(:,2)& U_new(:,1)>max(out-[m,0])));% |...
        %(c_new == 2 & (out(1)-m)>out(2) & max(U_new,[],2)<max(out-[m,0]));
    m_obs_2 = (c_new == 1 & (U_new(:,2) >  U_new(:,1)& U_new(:,2)>max(out-[0,m])));% |...
        %(c_new == 1 & (out(2)-m)>out(1) & max(U_new,[],2)<max(out-[0,m]));
    m_obs = m_obs_1|m_obs_2;
    s_new = (max(U_new,[],2)>max(out-[m,0].*(c_new==2)-[0,m].*(c_new==1),[],2));
    %% Determine the new population's wage, remove "old" people and merge
    % Get new cities
    c_next = c_new;
    c_next(c_new == 1 & m_obs ==1) = 2;
    c_next(c_new == 2 & m_obs ==1) = 1;
    w_new = zeros(N_new,1); % instantiate the wage
    % Get wages for the new population
    for i = 1:C 
        for j = 1:2
            w_new(c_next == i & s_new == j-1) ...
                = normrnd(w_bar(i,j),sigma,sum(c_next == i & s_new == j-1),1);
        end
    end
    % Now we eliminate the same fraction from the population as a whole and 
    keepers = randsample(1:length(c),round(.96*length(c))); 
    c_next = [c_next;c(keepers)];
    s_next = [s_new;s(keepers)];
    w_next = [w_new;w(keepers)];
    move = [sum(m_obs_1),sum(m_obs_2)];
    edu = sum(s_new);
end