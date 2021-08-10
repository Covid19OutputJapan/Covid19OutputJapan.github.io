function [betaT] = seasonal_adjustment(betaT,retro_lb,retro_ub,dateD,SimPeriod)

retro = mean(retro_lb,retro_ub);
cw = week(datetime(dateD(end),'ConvertFrom','excel')); % current week-of-year number
sine_lag = cw + 52 - 43 + 1;  % reference point in sine function is 43rd week of the year
seasonal_effect = 0.125;   % degree of seasonality
seasonality = 1 + seasonal_effect*sin((sine_lag*2*pi/52)+((2*pi)/52:(2*pi)/52:2*pi*SimPeriod/52));
seasonality_prev = 1 + seasonal_effect*sin(((sine_lag-retro)*2*pi/52)+((2*pi)/52:(2*pi)/52:2*pi*retro/52));
% seasonality = seasonality/seasonality(1);
seasonality = seasonality/mean(seasonality_prev);
betaT = betaT.*transpose(seasonality);

% beta_avg = mean(beta_prev / seasonality_prev);
% beta_avg .* seasonality;