function [deltaT,delta_bar] = construct_delta_variant(RetroPeriod,...
    retro_lb,retro_ub,deltaT,delta,delta_average,var_prev,var_share,var_infection_delta,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Delta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var_share_prev = var_prev(end-RetroPeriod+1:end,1);


% step 2: relative infection rate in the last 17 weeks
relative_var_infection_prev = 1 + var_share_prev * var_infection_delta;
% step 4:
no_var_delta = delta(end-RetroPeriod+1:end) ./ relative_var_infection_prev;
delta_bar = sum(no_var_delta.*(I(end-RetroPeriod+1:end)/sum(I(end-RetroPeriod+1:end))));
deltaT = deltaT.*(delta_bar/delta_average).*(1+var_infection_delta*var_share);
