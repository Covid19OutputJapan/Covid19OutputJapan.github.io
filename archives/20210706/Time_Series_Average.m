function [delta,beta_tilde,ERN,beta,ICU_inflow,...
    gammaT,delta_average,ICU_inflow_avg,delta_sample]...
    = Time_Series_Average(S,I,D,ICU,dD,N,Tdata,SimPeriod,...
        RetroPeriod,POP0,gamma,hconstant,h_all,alpha,k,...
        gamma_ICU,ICU_adjustment,RetroPeriodDelta)
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
    ICU_inflow = (ICU(2:Tdata+1) - ICU(1:Tdata) + gamma_ICU.*ICU(1:Tdata) + dD(1:Tdata))./(delta(1:Tdata).*I(1:Tdata));
    ICU_inflow_avg = mean(ICU_inflow(end-RetroPeriod+1:end))*ICU_adjustment;
