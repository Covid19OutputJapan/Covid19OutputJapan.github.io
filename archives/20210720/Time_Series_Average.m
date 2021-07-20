function [delta,beta_tilde,ERN,beta,ICU_inflow_nation, ICU_inflow_pref, Hospital_inflow,...
    gammaT,delta_average,ICU_inflow_avg_nation,ICU_inflow_avg_pref,Hospital_inflow_avg,delta_sample,beta_avg,beta_se,delta_se]...
    = Time_Series_Average(S,I,D,ICU_nation,ICU_pref,hospital,dD,N,Tdata,SimPeriod,...
        RetroPeriod,POP0,gamma,hconstant,h_all,alpha,k,...
        gamma_ICU,ICU_adjustment,gamma_Hospital,Hospital_adjustment,RetroPeriodDelta)
    %--- Compute the history of time-varying parameters ---%
    delta = (D(2:Tdata+1)-D(1:Tdata))./I(1:Tdata);                              % death rate
    beta_tilde = (POP0.*N(1:Tdata))./((S(1:Tdata).*I(1:Tdata)));   % overall infection rate, p4 of Fujii and Nakata (2020)
    ERN = (S(1:end-1)/POP0).*beta_tilde./(gamma+delta); % effective reproduction number
    if hconstant == 0
        beta = beta_tilde./(1+h_all*alpha).^k; % raw infection rate
    elseif hconstant == 1
        beta = beta_tilde./(1+(h_all(2)/h_all(1))*alpha).^k;
    end

    %--- Construct time series of parameters ---%
    gammaT = gamma*ones(SimPeriod,1);
    delta_sample = delta(end-RetroPeriodDelta+1:end);
    delta_average = sum(delta_sample.*(I(end-RetroPeriodDelta+1:end)/sum(I(end-RetroPeriodDelta+1:end))));
    ICU_inflow_nation = (ICU_nation(2:Tdata+1) - ICU_nation(1:Tdata) + gamma_ICU.*ICU_nation(1:Tdata) + dD(1:Tdata))./(delta(1:Tdata).*I(1:Tdata));
    ICU_inflow_avg_nation = mean(ICU_inflow_nation(end-RetroPeriod+1:end))*ICU_adjustment;
    ICU_inflow_pref = (ICU_pref(2:Tdata+1) - ICU_pref(1:Tdata) + gamma_ICU.*ICU_pref(1:Tdata) + dD(1:Tdata))./(delta(1:Tdata).*I(1:Tdata));
    ICU_inflow_avg_pref = mean(ICU_inflow_pref(end-RetroPeriod+1:end))*ICU_adjustment;
    Hospital_inflow = (hospital(2:Tdata+1) - hospital(1:Tdata) + gamma_Hospital.*hospital(1:Tdata))./(delta(1:Tdata).*I(1:Tdata));
    Hospital_inflow_avg = mean(Hospital_inflow(end-RetroPeriod+1:end))*Hospital_adjustment;
    
    %beta_r = 0;
    %for retrop = retro_lb:retro_ub
    %    beta_r = beta_r + mean(beta(end-retrop+1:end));
    %end
    %beta_avg = beta_r/(retro_ub-retro_lb+1);
    beta_sample = beta(end-RetroPeriodDelta+1:end);
    beta_avg = mean(beta_sample);
    beta_se  = std(beta_sample)/sqrt(length(beta_sample));
    delta_se= std(delta_sample)/sqrt(length(delta_sample));