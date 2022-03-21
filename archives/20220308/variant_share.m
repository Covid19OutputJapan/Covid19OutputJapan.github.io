function share = variant_share(ini, ss, growth, SimPeriod)

% ini             = 0.0001;      % last period value in the data (at Tdata)
% ss              = 1;
% growth          = 0.8;
% SimPeriod       = 52; 
tvec            = 1:SimPeriod;
share           = zeros(SimPeriod, 1);

logit_initial   = log(ini / (ss - ini) ); %logit_initial = log(share_0 / (ss - share_0))



X = exp(logit_initial +  tvec'* growth );
share(1:end)       = (X .* ss) ./  (1 + X );

% figure(1)
% plot(share)
end

% end