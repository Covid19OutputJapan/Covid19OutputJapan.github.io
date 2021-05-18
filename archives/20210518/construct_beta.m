function [betaT,betaT_woAR1,var_share,var_share_mat] = construct_beta(home,SimPeriod,Tdata,dateEN,pref,...
    retro_lb,retro_ub,beta,beta_avg,var_initial,var_ss,var_infection,var_growth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betaT = beta_avg*ones(SimPeriod,1);

%for iS = 1:length(D1_vec)

var_intercept = log(var_initial/(var_ss-var_initial)); % Logit of the variant share, most recently
% var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept

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
var_share_prev = zeros(retro_lb,1);
var_share_prev = var_share_mat(end-retro_lb+1:end,1);

% step 2: relative infection rate in the last 17 weeks
relative_var_infection_prev = 1 + var_share_prev * var_infection;
% step 3:
no_var_beta = beta(end-mean(retro_lb,retro_ub)+1:end) ./ relative_var_infection_prev;
beta_bar = mean(no_var_beta);
betaT = beta_bar*(1+var_infection*var_share);
betaT_woAR1 = betaT;
