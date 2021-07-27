function [delta,beta_tilde,ERN,beta,ICU_nation_inflow, ICU_pref_inflow, Hospital_inflow,...
    gammaT,delta_average,delta_ICU_nation_average,delta_ICU_pref_average,delta_Hospital_average,...
    ICU_nation_inflow_avg,ICU_pref_inflow_avg,Hospital_inflow_avg,beta_avg,beta_se,...
    delta_se, delta_ICU_nation_se, delta_ICU_pref_se, delta_Hospital_se]...
    = Time_Series_Average(S,I,D,ICU_nation,ICU_pref,hospital,dD,N,Tdata,SimPeriod,...
        RetroPeriod,POP0,hconstant,h_all,alpha,k,...
        gamma,gamma_ICU_nation,gamma_ICU_pref,gamma_Hospital,...
        ICU_nation_adjustment,ICU_pref_adjustment,Hospital_adjustment,...
        RetroPeriodDelta,RetroPeriodICU_nation,RetroPeriodICU_pref,RetroPeriodHospital)
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
    weight = I(end-RetroPeriodDelta+1:end)/sum(I(end-RetroPeriodDelta+1:end));
    weight_delta_sample = delta(end-RetroPeriodDelta+1:end).*weight;
    delta_average = sum(weight_delta_sample);
    delta_se= sqrt(var(delta_sample)*sum((weight/sum(weight)).^2));
    
    % ICU (National standard)
    % delta_ICU_nation_sample = delta(end-RetroPeriodICU_nation+1:end);
    % weight = I(end-RetroPeriodICU_nation+1:end)/sum(I(end-RetroPeriodICU_nation+1:end));
    % delta_ICU_nation_average = sum(delta_ICU_nation_sample.*(I(end-RetroPeriodICU_nation+1:end)/sum(I(end-RetroPeriodICU_nation+1:end))));
    ICU_nation_inflow = (ICU_nation(2:Tdata+1) - ICU_nation(1:Tdata) + gamma_ICU_nation.*ICU_nation(1:Tdata) + dD(1:Tdata))./(delta(1:Tdata).*I(1:Tdata));
    ICU_nation_inflow_avg = mean(ICU_nation_inflow(end-RetroPeriodICU_nation+1:end))*ICU_nation_adjustment;
    delta_ICU_nation = delta;
    delta_ICU_nation_sample = delta_ICU_nation(end-RetroPeriodICU_nation+1:end);
    weight = I(end-RetroPeriodICU_nation+1:end)/sum(I(end-RetroPeriodICU_nation+1:end));
    weight_delta_ICU_nation_sample = delta_ICU_nation(end-RetroPeriodICU_nation+1:end).*weight;
    delta_ICU_nation_average = sum(weight_delta_ICU_nation_sample);
    delta_ICU_nation_se= sqrt(var(delta_ICU_nation_sample)*sum((weight/sum(weight)).^2));

    % ICU (Tokyo standard)
    % delta_ICU_pref_sample = delta(end-RetroPeriodICU_pref+1:end);
    % delta_ICU_pref_average = sum(delta_ICU_pref_sample.*(I(end-RetroPeriodICU_pref+1:end)/sum(I(end-RetroPeriodICU_pref+1:end))));
    % weight_ICU
    ICU_pref_inflow = (ICU_pref(2:Tdata+1) - ICU_pref(1:Tdata) + gamma_ICU_pref.*ICU_pref(1:Tdata) + dD(1:Tdata))./(delta(1:Tdata).*I(1:Tdata));
    ICU_pref_inflow_avg = mean(ICU_pref_inflow(end-RetroPeriodICU_pref+1:end))*ICU_pref_adjustment;
    
    delta_ICU_pref = delta;
    delta_ICU_pref_sample = delta_ICU_pref(end-RetroPeriodICU_pref+1:end);
    weight = I(end-RetroPeriodICU_pref+1:end)/sum(I(end-RetroPeriodICU_pref+1:end));
    weight_delta_ICU_pref_sample = delta_ICU_pref(end-RetroPeriodICU_pref+1:end).*weight;
    delta_ICU_pref_average = sum(weight_delta_ICU_pref_sample);
    delta_ICU_pref_se= sqrt(var(delta_ICU_pref_sample)*sum((weight/sum(weight)).^2));

    % Hospital
    % delta_Hospital_sample = delta(end-RetroPeriodHospital+1:end);
    % delta_Hospital_average = sum(delta_Hospital_sample.*(I(end-RetroPeriodHospital+1:end)/sum(I(end-RetroPeriodHospital+1:end))));
    delta_Hospital = delta;
    delta_Hospital_sample = delta_Hospital(end-RetroPeriodHospital+1:end);
    weight = I(end-RetroPeriodHospital+1:end)/sum(I(end-RetroPeriodHospital+1:end));
    weight_delta_Hospital_sample = delta_Hospital(end-RetroPeriodHospital+1:end).*weight;
    delta_Hospital_average = sum(weight_delta_Hospital_sample);
    delta_Hospital_se= sqrt(var(delta_Hospital_sample)*sum((weight/sum(weight)).^2));

    Hospital_inflow = (hospital(2:Tdata+1) - hospital(1:Tdata) + gamma_Hospital.*hospital(1:Tdata))./(delta(1:Tdata).*I(1:Tdata));
    Hospital_inflow_avg = mean(Hospital_inflow(end-RetroPeriodHospital+1:end))*Hospital_adjustment;
    
    %beta_r = 0;
    %for retrop = retro_lb:retro_ub
    %    beta_r = beta_r + mean(beta(end-retrop+1:end));
    %end
    %beta_avg = beta_r/(retro_ub-retro_lb+1);
    beta_sample = beta(end-RetroPeriod+1:end);
    beta_avg = mean(beta_sample);
    beta_se  = std(beta_sample)/sqrt(length(beta_sample));
%     delta_se= std(delta_sample)/sqrt(length(delta_sample));
    