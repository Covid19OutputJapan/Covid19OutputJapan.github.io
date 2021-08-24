function [delta_wo_alpha,delta_ICU_pref_wo_alpha, delta_Hospital_wo_alpha,...
    delta_ICU_nation_wo_alpha, beta_avg, beta_wo_alpha, beta_se,...
    delta_average, delta_ICU_nation_average, delta_ICU_pref_average, delta_Hospital_average, ...
    var_share, var_prev, var_initial, var_share2, var_prev2, var_initial2] ...
    = eliminate_variant_effects(delta_average, delta_ICU_nation_average, ...
    delta_ICU_pref_average, delta_Hospital_average, ...
    beta, SimPeriod, var_ss, var_growth, var_infection, var_infection_delta, var_growth2, var_infection2, var_infection_delta2,...
    RetroPeriod, RetroPeriodDelta, RetroPeriodICU_nation, RetroPeriodICU_pref, RetroPeriodHospital,...
    Data)

%--- Eliminate the effects of vaccination from delta ---%
% delta_ss = delta_average * (0.1063/1.53); %Share of the death rate among youth to the death rate of all populaiton
% VD_elderly = D1 * V1_elderly + (D2 - D1) * V2_elderly;
% VD_ordinary = (D1 * V1_medical + (D2 - D1) * V2_medical) + (D1 * V1_others + (D2 - D1) * V2_others);
% share = ((1 - sum(VD_elderly(1:end - 2)) / elderly_total) * (delta_average - delta_ss) ...
%     + (1 - sum(VD_ordinary(1:end - 2)) / (ordinary_total)) * delta_ss) / delta_average;
% delta_wo_vaccination = delta_average / share; %Eliminate the effects of vaccination in the past average value, Past 17 without vaccine effects
% 
% %--- Eliminate the effects of vaccination from delta_ICU_nation ---%
% delta_ICU_nation_ss = delta_ICU_nation_average * (0.3916/1.62); %Share of the death rate among youth to the death rate of all populaiton
% share = ((1 - sum(VD_elderly(1:end - 2)) / elderly_total) * (delta_ICU_nation_average - delta_ICU_nation_ss) ...
%     + (1 - sum(VD_ordinary(1:end - 2)) / (ordinary_total)) * delta_ICU_nation_ss) / delta_ICU_nation_average;
% delta_ICU_nation_wo_vaccination = delta_ICU_nation_average / share; %Eliminate the effects of vaccination in the past average value
% 
% %--- Eliminate the effects of vaccination from delta_ICU_pref ---%
% delta_ICU_pref_ss = delta_ICU_pref_average * (0.3916/1.62); %Share of the death rate among youth to the death rate of all populaiton
% share = ((1 - sum(VD_elderly(1:end - 2)) / elderly_total) * (delta_ICU_pref_average - delta_ICU_pref_ss) ...
%     + (1 - sum(VD_ordinary(1:end - 2)) / (ordinary_total)) * delta_ICU_pref_ss) / delta_ICU_pref_average;
% delta_ICU_pref_wo_vaccination = delta_ICU_pref_average / share; %Eliminate the effects of vaccination in the past average value
% 
% %--- Eliminate the effects of vaccination from delta_Hospital ---%
% delta_Hospital_ss = delta_Hospital_average * (0.3916/1.62); %Share of the death rate among youth to the death rate of all populaiton
% share = ((1 - sum(VD_elderly(1:end - 2)) / elderly_total) * (delta_Hospital_average - delta_Hospital_ss) ...
%     + (1 - sum(VD_ordinary(1:end - 2)) / (ordinary_total)) * delta_Hospital_ss) / delta_Hospital_average;
% delta_Hospital_wo_vaccination = delta_Hospital_average / share; %Eliminate the effects of vaccination in the past average value

%--- Eliminate the effects of alpha variant from delta and beta---%
% [var_share, var_prev, var_initial] = var_share_prev(Data(:, 20), SimPeriod, var_ss, var_growth);
%     var_infection_adjustment = (1 - mean(var_prev(end - RetroPeriod + 1:end))) * 1 ...
%         + mean(var_prev(end - RetroPeriod + 1:end)) * (1 + var_infection); %Relative increase of infectiousness (alpha varaint, past 17 weeks)
%     var_infection_delta_adjustment = (1 - mean(var_prev(end - RetroPeriodDelta + 1:end))) * 1 ...
%         + mean(var_prev(end - RetroPeriodDelta + 1:end)) * (1 + var_infection_delta); %Relative increase of death rate (alpha varaint, past 17 weeks)
%     var_infection_delta_ICU_nation_adjustment = (1 - mean(var_prev(end - RetroPeriodICU_nation + 1:end))) * 1 ...
%         + mean(var_prev(end - RetroPeriodICU_nation + 1:end)) * (1 + var_infection_delta); %Relative increase of death rate (alpha varaint, past 17 weeks)
%     var_infection_delta_ICU_pref_adjustment = (1 - mean(var_prev(end - RetroPeriodICU_pref + 1:end))) * 1 ...
%         + mean(var_prev(end - RetroPeriodICU_pref + 1:end)) * (1 + var_infection_delta); %Relative increase of death rate (alpha varaint, past 17 weeks)
%     var_infection_delta_Hospital_adjustment = (1 - mean(var_prev(end - RetroPeriodHospital + 1:end))) * 1 ...
%         + mean(var_prev(end - RetroPeriodHospital + 1:end)) * (1 + var_infection_delta); %Relative increase of death rate (alpha varaint, past 17 weeks)
%
%     delta_wo_alpha = delta_wo_vaccination / var_infection_delta_adjustment; %Eliminate the effects of alpha variant in the past average value
%     delta_ICU_nation_wo_alpha = delta_ICU_nation_wo_vaccination / var_infection_delta_ICU_nation_adjustment;
%     delta_ICU_pref_wo_alpha = delta_ICU_pref_wo_vaccination / var_infection_delta_ICU_pref_adjustment;
%     delta_Hospital_wo_alpha = delta_Hospital_wo_vaccination / var_infection_delta_Hospital_adjustment;

delta_wo_alpha = delta_average /(1 + var_infection_delta); %Eliminate the effects of alpha variant in the past average value
delta_ICU_nation_wo_alpha = delta_ICU_nation_average /(1+var_infection_delta);
delta_ICU_pref_wo_alpha = delta_ICU_pref_average /(1+var_infection_delta);
delta_Hospital_wo_alpha = delta_Hospital_average /(1+var_infection_delta);

%     beta_wo_alpha = beta_avg / var_infection_adjustment; %Eliminate the effects of alpha variatn

beta_wo_alpha_series = beta/(1+var_infection);
beta_wo_alpha = mean(beta_wo_alpha_series(end - RetroPeriod + 1:end));



%--- Eliminate the effects of delta variant from delta and beta---%

%--- Eliminate the effects of delta variant from beta---%
[var_share2, var_prev2, var_initial2] = var_share_prev(Data(:, 24), SimPeriod, var_ss, var_growth2);
[beta_wo_delta_series, beta_avg,beta_se] = variant_adjustment(beta_wo_alpha_series,var_prev2,RetroPeriod,var_infection2);
[delta_wo_delta_series, delta_average,delta_se] = variant_adjustment(delta_wo_alpha,var_prev2,RetroPeriod,var_infection2);
[delta_ICU_nation_wo_delta_series, delta_ICU_nation_average,delta_ICU_nation_se] = variant_adjustment(delta_ICU_nation_wo_alpha,var_prev2,RetroPeriod,var_infection2);
[delta_ICU_pref_wo_delta_series, delta_ICU_pref_average,delta_ICU_pref_se] = variant_adjustment(delta_ICU_pref_wo_alpha,var_prev2,RetroPeriod,var_infection2);
[delta_Hospital_wo_delta_series, delta_Hospital_average,delta_Hospital_se] = variant_adjustment(delta_Hospital_wo_alpha,var_prev2,RetroPeriod,var_infection2);


% var_infection_adjustment2 = (1 - mean(var_prev2(end - RetroPeriod + 1:end))) * 1 ...
%     + mean(var_prev2(end - RetroPeriod + 1:end)) * (1 + var_infection2); %Relative increase of infectiousness (alpha varaint, past 17 weeks)
% var_infection_delta_adjustment2 = (1 - mean(var_prev2(end - RetroPeriodDelta + 1:end))) * 1 ...
%     + mean(var_prev2(end - RetroPeriodDelta + 1:end)) * (1 + var_infection_delta2); %Relative increase of death rate (alpha varaint, past 17 weeks)
% var_infection_delta_ICU_nation_adjustment2 = (1 - mean(var_prev2(end - RetroPeriodICU_nation + 1:end))) * 1 ...
%     + mean(var_prev2(end - RetroPeriodICU_nation + 1:end)) * (1 + var_infection_delta2);
% var_infection_delta_ICU_pref_adjustment2 = (1 - mean(var_prev2(end - RetroPeriodICU_pref + 1:end))) * 1 ...
%     + mean(var_prev2(end - RetroPeriodICU_pref + 1:end)) * (1 + var_infection_delta2);
% var_infection_delta_Hospital_adjustment2 = (1 - mean(var_prev2(end - RetroPeriodHospital + 1:end))) * 1 ...
%     + mean(var_prev2(end - RetroPeriodHospital + 1:end)) * (1 + var_infection_delta2);
% 
% delta_average = delta_wo_alpha / var_infection_delta_adjustment2; %Eliminate the effects of delta variant in the past average value
% delta_ICU_nation_average = delta_ICU_nation_wo_alpha / var_infection_delta_ICU_nation_adjustment2;
% delta_ICU_pref_average = delta_ICU_pref_wo_alpha / var_infection_delta_ICU_pref_adjustment2;
% delta_Hospital_average = delta_Hospital_wo_alpha / var_infection_delta_Hospital_adjustment2;





end