function [betaT,betaT_woAR1,beta_bar] = construct_beta(SimPeriod,...
    retro_lb,retro_ub,beta,beta_avg,var_prev,var_share,var_infection)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaT = beta_avg*ones(SimPeriod,1);

var_share_prev = var_prev(end-mean(retro_lb,retro_ub)+1:end,1);

% step 2: relative infection rate in the last 17 weeks
relative_var_infection_prev = 1 + var_share_prev * var_infection;
% step 3:
no_var_beta = beta(end-mean(retro_lb,retro_ub)+1:end) ./ relative_var_infection_prev;
beta_bar = mean(no_var_beta);
betaT = beta_bar*(1+var_infection*var_share);
betaT_woAR1 = betaT;
