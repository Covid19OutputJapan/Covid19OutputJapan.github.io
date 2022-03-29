function  share = variant_share_two_target(T1, T2, target_shareT1, target_shareT2, ss, SimPeriod) 

tvec            = 1:SimPeriod;
share           = zeros(SimPeriod, 1);



% Find growth rate so that s_T = x when ss, share0, and T is giveN:
% T1 = 7;       %SimPeriod + T1
% T2 = 10;
% target_share1    =        0.05;
% target_share2    =        0.8;

target_initial   =      log(target_shareT1 / (ss - target_shareT1) );
growth_target   =       (log(target_shareT2 / (ss - target_shareT2) )  -  target_initial)/(T2-T1);

logit_target_initial_value   =     log(target_shareT1/(1-target_shareT1))  - growth_target*T1;


logit_initial   =   logit_target_initial_value;
growth          =   growth_target;


X = exp(logit_initial +  tvec'* growth );
share(1:end)       = (X .* ss) ./  (1 + X );

% figure(1)
% plot(share)

end