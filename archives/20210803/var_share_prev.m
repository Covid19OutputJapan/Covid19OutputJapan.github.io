function [var_share, var_prev, var_initial] = var_share_prev(Data,SimPeriod,...
    var_ss,var_growth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for iS = 1:length(D1_vec)
% if iPC==1
%     Data = importdata([home '\Covid_weekly_newV.csv']);  % Import weekly Covid data by prefecture
% else
%     Data = importdata([home 'Covid_weekly_newV.csv']);  % Import weekly Covid data by prefecture
% end
% variant_data = Data.data(strcmp(Data.textdata(2:end,1),pref),:); %same length as covid_weekly.csv
var_prev = Data;
var_prev(isnan(var_prev))=0;
ind_var = find(var_prev>0,1,'last');
logit_initial = log(var_prev(ind_var)/(var_ss-var_prev(ind_var)));
if ind_var < length(var_prev)
    for i = ind_var+1:length(var_prev)
        var_prev(i) = exp((i-ind_var)*var_growth+logit_initial).*var_ss./(1+exp((i-ind_var)*var_growth+logit_initial));
    end
end
var_initial = var_prev(end);

var_intercept = log(var_initial/(var_ss-var_initial)); % Logit of the variant share, most recently
% var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept

% Simulated share of variants
var_share = exp((1:SimPeriod)'*var_growth+var_intercept).*var_ss./(1+exp((1:SimPeriod)'*var_growth+var_intercept));


