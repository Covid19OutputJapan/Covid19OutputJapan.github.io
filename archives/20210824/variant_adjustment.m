function [adjusted_values,adjusted_average,adjusted_se] = variant_adjustment(values,var_prev,RetroPeriod,var_infection)
%--- Eliminate the effects of variant---%
adjusted_values = values ./ (1+var_infection.*var_prev);
adjusted_average =  mean(adjusted_values(end - RetroPeriod + 1:end));
adjusted_se  = std(adjusted_values(end - RetroPeriod + 1:end))/sqrt(length(adjusted_values(end - RetroPeriod + 1:end)));
