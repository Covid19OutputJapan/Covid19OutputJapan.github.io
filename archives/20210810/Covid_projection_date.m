function [CumD,GDPLoss,alphapath,SimData,SimN,SimERN,SimICU_nation,SimICU_pref,SimHospital,betaAR1] = ...
    Covid_projection_date(InitialValues,alpha_on,alpha_off,...
    betaT,gammaT,deltaT,delta_ICU_nation,delta_ICU_pref,delta_Hospital,V,h,k,POP0,...
    hconstant,alpha_duration,ICU_nation_inflow_avg, ICU_pref_inflow_avg, Hospital_inflow_avg, ...
    gamma_ICU_nation, gamma_ICU_pref, gamma_Hospital,beta_shock_after_emergency,rho_after_emergency,alpha_jump,th_off_date)
% Simulation for cumulative deaths and output loss given parameters
% InitialValues = [S(end),I(end),R(end),D(end),ICU_nation(end),ICU_pref(end)];
% - InitialValues = 1-by-4 vector of initial values: S, I, R, D, ICU(national definition), ICU(prefecture-specific definitoin)
% - alpha_on is the target alpha with status of emergency
% - beta_sim, gammaT, deltaT and V are T-by-1 vectors of time-varying
% parameters
% - h, k are scalar parameters
% - T is a simulation period
%beta_shock = 1.0;

beta_shock = beta_shock_after_emergency;
beta_sim = betaT;
T = length(beta_sim);
SimData = zeros(T+1,length(InitialValues));
SimData(1,:) = InitialValues;
SimN = zeros(T,1);
SimICU_nation = zeros(T+1,1);
SimICU_pref = zeros(T+1,1);
SimHospital = zeros(T+1,1);
SimICU_nation(1) = InitialValues(5);
SimICU_pref(1) = InitialValues(6);
SimHospital(1) = InitialValues(7);

alpha_off_path = alpha_on+(alpha_off-alpha_on)*alpha_jump+(1-alpha_jump)*((alpha_off-alpha_on)/(alpha_duration+1)):(1-alpha_jump)*((alpha_off-alpha_on)/(alpha_duration+1)):alpha_off;
delta_init = deltaT(1);
% coef_th_on_path = (th_on2 - th_on1)/(1-deltaT(T)/delta_init);
% THonPath = th_on1 + coef_th_on_path.*(1-deltaT./delta_init);

alphapath = alpha_off .* ones(T,1);
alphapath(1:th_off_date-1) = alpha_on;
alphapath(th_off_date:th_off_date+alpha_duration) = alpha_off_path;

for t = th_off_date:T
    beta_sim(t) = beta_sim(t)*(1+beta_shock);
    beta_shock = beta_shock * rho_after_emergency;
end

for t = 1:T
    alpha_sim = alphapath(t);
    
    if hconstant == 0
        SimN(t) = ((1 + h*alpha_sim)^k)*beta_sim(t)*SimData(t,1)*SimData(t,2)*(1/POP0);
    elseif hconstant == 1
        SimN(t) = ((1+(h(2)/h(1))*alpha_sim)^k)*beta_sim(t)*SimData(t,1)*SimData(t,2)*(1/POP0);
    end
    SimData(t+1,1) = SimData(t,1) - SimN(t) - V(t);
    SimData(t+1,2) = SimData(t,2) + SimN(t) - gammaT(t)*SimData(t,2) - deltaT(t)*SimData(t,2);
    SimData(t+1,3) = SimData(t,3) + gammaT(t)*SimData(t,2) + V(t);
    SimData(t+1,4) = SimData(t,4) + deltaT(t)*SimData(t,2);
    SimICU_nation(t+1) = max(0, SimICU_nation(t) + ICU_nation_inflow_avg*delta_ICU_nation(t)*SimData(t,2) ...
        - gamma_ICU_nation*SimICU_nation(t) - deltaT(t)*SimData(t,2));
    SimICU_pref(t+1) = max(0, SimICU_pref(t) ...
        + ICU_pref_inflow_avg*delta_ICU_pref(t)*SimData(t,2) ...
        - gamma_ICU_pref*SimICU_pref(t) - deltaT(t)*SimData(t,2));
    SimHospital(t+1) = max(0, SimHospital(t) + Hospital_inflow_avg*delta_Hospital(t)*SimData(t,2) ...
        - gamma_Hospital*SimHospital(t));
    
end
SimERN = (SimData(1:end-1,1)./POP0).*(((1+(h(2)/h(1))*alphapath).^k).*beta_sim)./(gammaT+deltaT);
CumD = SimData(end,4);              % Cumulative deaths during the simulation period
GDPLoss = mean(alphapath);              % Average output loss during the simulation period
betaAR1 = beta_sim;
