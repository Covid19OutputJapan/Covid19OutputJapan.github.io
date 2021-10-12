function [deltaT] = construct_delta_variant(home,SimPeriod,RetroPeriod,Tdata,dateEN,pref,...
    retro_lb,retro_ub,deltaT,delta,delta_average,var_initial,var_ss,var_infection_delta,var_growth,I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Delta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var_intercept = log(var_initial/(var_ss-var_initial)); % Logit of the variant share, most recently
% step 1: extrapolation in the variant share in the last 17 weeks     -mean(retro_lb,retro_ub)+1 = -17 + 1
var_share = exp((1:SimPeriod)'*var_growth+var_intercept).*var_ss./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
var_share_data = readtable([home 'Variant_rate.xls']);
date_column = table2array(var_share_data(:,1));
pref_column = table2array(var_share_data(:,2));
var_share_mat = zeros(Tdata,1);
ct = 1;
for i = 1:length(date_column)
    if pref_column{i} == convertCharsToStrings(pref)
        var_share_prev(ct,1) = table2array(var_share_data(i,5));
        var_date_prev(ct,1) = table2array(var_share_data(i,1));
        ct = ct + 1;
    end
end
ind_var1 = find(dateEN > datetime(var_date_prev(1,1)),1,'first');
ind_var2 = find(dateEN > datetime(var_date_prev(end,1)),1,'first');
var_share_mat(ind_var1-1:ind_var2-1,1) = var_share_prev(:,1);
logit_initial = log(var_share_mat(ind_var2-1)/(var_ss-var_share_mat(ind_var2-1)));
for i = ind_var2:length(var_share_mat)
    var_share_mat(i) = exp((i-ind_var2+1)*var_growth+logit_initial).*var_ss./(1+exp((i-ind_var2+1)*var_growth+logit_initial));
end
var_share_prev = var_share_mat(end-retro_lb+1:end,1);

% step 2: relative infection rate in the last 17 weeks
relative_var_infection_prev = 1 + var_share_prev * var_infection_delta;
% step 4:
no_var_delta = delta(end-mean(retro_lb,retro_ub)+1:end) ./ relative_var_infection_prev;
delta_bar = sum(no_var_delta.*(I(end-RetroPeriod+1:end)/sum(I(end-RetroPeriod+1:end))));
deltaT = deltaT.*(delta_bar/delta_average).*(1+var_infection_delta*var_share);
