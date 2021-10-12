retro = mean(retro_lb,retro_ub);
cw = week(datetime(dateD(end),'ConvertFrom','excel')); % current week-of-year number
sine_lag = cw + 52 - 43;  % reference point in sine function is 43rd week of the year
seasonal_effect = 0.1;   % degree of seasonality
seasonality = 1 + seasonal_effect*sin((sine_lag*2*pi/52)+((2*pi)/52:(2*pi)/52:2*pi*SimPeriod/52));
seasonality_prev = 1 + seasonal_effect*sin(((sine_lag-retro)*2*pi/52)+((2*pi)/52:(2*pi)/52:2*pi*retro/52));
seasonality = seasonality/seasonality(1);
betaT = betaT.*seasonality;