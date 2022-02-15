function [CumD,GDPLoss,alphaPath,SimData,SimN,SimERN,SimICU_nation,SimICU_pref,SimHospital,betaPath] = ...
    Covid_projection_withSOE(InitialValues,alpha_on,alpha_off,th_on,...
    th_off,betaT,gamma,deltaT,delta_ICU_nation,delta_ICU_pref,delta_Hospital,V,h,k,POP0,...
    hconstant,DRi,ICU_nation_inflow_avg, ICU_pref_inflow_avg, Hospital_inflow_avg, ...
    gamma_ICU_nation, gamma_ICU_pref, gamma_Hospital,...
    relative_beta_SOE, beta_jump, beta_goal, seasonality, alphaBox, state, simple_beta_avg)
% Simulation for cumulative deaths and output loss given parameters
% InitialValues = [S(end),I(end),R(end),D(end),ICU_nation(end),ICU_pref(end)];
% - InitialValues = 1-by-4 vector of initial values: S, I, R, D, ICU(national definition), ICU(prefecture-specific definitoin)
% - alpha_on is the target alpha with status of emergency
% - beta, gamma, delta and V are T-by-1 vectors of time-varying
% parameters
% - h, k are scalar parameters
% - T is a simulation period
%beta_shock = 1.0;

% betaT = beta;

betaT_Path = betaT;
alpha_Path = alphaBox;
alpha_duration = DRi;
T = length(betaT);
betaPath = zeros(T,1);
alphaPath = zeros(T,1);

SimData = zeros(T+1,length(InitialValues));
SimData(1,:) = InitialValues;
SimN = zeros(T,1);
SimICU_nation = zeros(T+1,1);
SimICU_pref = zeros(T+1,1);
SimHospital = zeros(T+1,1);
SimICU_nation(1) = InitialValues(5);
SimICU_pref(1) = InitialValues(6);
SimHospital(1) = InitialValues(7);

for t = 1:T
    beta = betaT_Path(t);
    alpha = alpha_Path(t);

    if hconstant == 0
        SimN(t) = ((1 + h*alpha)^k)*beta*SimData(t,1)*SimData(t,2)*(1/POP0);
    elseif hconstant == 1
        SimN(t) = ((1+(h(2)/h(1))*alpha)^k)*beta*SimData(t,1)*SimData(t,2)*(1/POP0);
    end
    
    SimData(t+1,1) = SimData(t,1) - SimN(t) - V(t);
    SimData(t+1,2) = SimData(t,2) + SimN(t) - gamma(t)*SimData(t,2) - deltaT(t)*SimData(t,2);
    SimData(t+1,3) = SimData(t,3) + gamma(t)*SimData(t,2) + V(t);
    SimData(t+1,4) = SimData(t,4) + deltaT(t)*SimData(t,2);
    SimICU_nation(t+1) = max(0, SimICU_nation(t) + ICU_nation_inflow_avg*delta_ICU_nation(t)*SimData(t,2) ...
        - gamma_ICU_nation*SimICU_nation(t) - deltaT(t)*SimData(t,2));
    SimICU_pref(t+1) = max(0, SimICU_pref(t) ...
        + ICU_pref_inflow_avg*delta_ICU_pref(t)*SimData(t,2) ...
        - gamma_ICU_pref*SimICU_pref(t) - deltaT(t)*SimData(t,2));
    SimHospital(t+1) = max(0, SimHospital(t) + Hospital_inflow_avg*delta_Hospital(t)*SimData(t,2) ...
        - gamma_Hospital*SimHospital(t));
    
        
    if (SimN(t) >= th_on) 
        state = 1;
        alpha_Path = alpha_on*ones(T,1);
        betaT_Path = relative_beta_SOE * simple_beta_avg*ones(T,1);
        
    elseif (SimN(t) < th_off) 
        if state == 1
            betaT_Path(t + 1:t + alpha_duration + 1) = ...
                linspace(beta_jump*beta_goal,beta_goal,alpha_duration + 1).*seasonality(t + 1:t + alpha_duration + 1);
            alpha_Path(t + 1:t + alpha_duration + 1) = ...
                   [alpha_on:(alpha_off - alpha_on)/(alpha_duration):alpha_off];
            betaT_Path(t + alpha_duration + 2:end) = betaT(t + alpha_duration + 2:end);
            alpha_Path(t + alpha_duration + 2:end) = alpha_off;
        end
        
        state = 0;
    end
        
    betaPath(t) = beta;
    alphaPath(t) = alpha;
end

SimERN = (SimData(1:end-1,1)./POP0).*(((1+(h(2)/h(1))*alphaPath).^k).*betaPath)./(gamma+deltaT);
CumD = SimData(end,4);              % Cumulative deaths during the simulation period
GDPLoss = mean(alphaPath);              % Average output loss during the simulation period
