function [CumD,GDPLoss,alphaPath,SimData,SimN,SimERN,SimICU_nation,SimICU_pref,SimHospital,betaPath,beta_tilde_path,SimBRN] = ...
    Covid_projection_omicron(InitialValues,alpha_on,alpha_off,th_on,...
    th_off,betaT,gamma,deltaT,SimICU_nation_rate, SimICU_pref_rate, SimHospital_rate,h,k,POP0,...
    lagged_cumsumVT1, lagged_cumsumVT2, lagged_cumsumVT3, E1, E2, E3, cum_in_R, ...
    hconstant,DRi, ...
    gamma_ICU_nation, gamma_ICU_pref, gamma_Hospital,...
    scale_beta, beta_jump, beta_goal, seasonality, alphaBox, state, simple_beta_avg)
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
cum_in_simI    = 0;
for t = 1:T    
    beta = betaT_Path(t);
    alpha = alpha_Path(t);
    if E2(t) < E1(t)
        error(['E2 is lower than E1', ', t = ' num2str(t), ', E2 = ' num2str(E2(t)), ' E1 = ',  num2str(E1)])
    end
    SimData(t,3) = min(E1(t) * lagged_cumsumVT1(t) ...
                 + (E2(t) - E1(t)) * lagged_cumsumVT2(t) ...
                 + (E3(t) - E2(t)) * lagged_cumsumVT3(t) ...
                 + cum_in_R + cum_in_simI,POP0); %R = laggedV1*E1 + laggedV2 * E2 + laggedV3 * E3 + sum(\gamma * I_past)
    SimData(t,1) = max(0,POP0 - sum(SimData(t,2:4)));  %S = Pop0 - (I+R+D)
    cum_in_simI  = cum_in_simI + gamma(t) * SimData(t,2);
    if hconstant == 0
        SimN(t) = ((1 + h*alpha)^k)*beta*SimData(t,1)*SimData(t,2)*(1/POP0);
    elseif hconstant == 1
        SimN(t) = ((1+(h(2)/h(1))*alpha)^k)*beta*SimData(t,1)*SimData(t,2)*(1/POP0);
    end
    
    SimData(t+1,2) = SimData(t,2) + SimN(t) - gamma(t)*SimData(t,2) - deltaT(t)*SimData(t,2); %I
    SimData(t+1,4) = SimData(t,4) + deltaT(t)*SimData(t,2); %D
    
    SimICU_nation(t+1) = max(0, SimICU_nation(t) + SimICU_nation_rate(t)*SimData(t,2) ...
        - gamma_ICU_nation*SimICU_nation(t) - deltaT(t)*SimData(t,2));
    SimICU_pref(t+1) = max(0, SimICU_pref(t) ...
        + SimICU_pref_rate(t)*SimData(t,2) ...
        - gamma_ICU_pref*SimICU_pref(t) - deltaT(t)*SimData(t,2));
    SimHospital(t+1) = max(0, SimHospital(t) + SimHospital_rate(t)*SimData(t,2) ...
        - gamma_Hospital*SimHospital(t));
    
        
    if (SimICU_nation(t+1) >= th_on)  %(SimN(t) >= th_on)  %
        state = 1;
        alpha_Path = alpha_on*ones(T,1);
        betaT_Path = scale_beta * simple_beta_avg*ones(T,1);
        
    elseif (SimICU_nation(t+1) < th_off)%(SimN(t) < th_off)
        if state == 1
            if t + alpha_duration + 1 > T
                remainT = T-t;
                nextBetaPath(1:alpha_duration+1) = linspace(beta_jump*beta_goal(t),beta_goal(t),alpha_duration + 1).*seasonality(t + 1:t + alpha_duration + 1);
                betaT_Path(t + 1:end) = nextBetaPath(1:remainT);
                nextAlphaPath(1:alpha_duration+1) = [alpha_on:(alpha_off - alpha_on)/(alpha_duration):alpha_off];
                alpha_Path(t + 1:end) = nextAlphaPath(1:remainT);
            else
                betaT_Path(t + 1:t + alpha_duration + 1) = ...
                    linspace(beta_jump*beta_goal(t+1),beta_goal(t + alpha_duration + 1),alpha_duration + 1).*seasonality(t + 1:t + alpha_duration + 1);
                alpha_Path(t + 1:t + alpha_duration + 1) = ...
                   [alpha_on:(alpha_off - alpha_on)/(alpha_duration):alpha_off];
                betaT_Path(t + alpha_duration + 2:end) = betaT(t + alpha_duration + 2:end);
                alpha_Path(t + alpha_duration + 2:end) = alpha_off;
            end
        end
        
        state = 0;
    end
        
    betaPath(t) = beta;
    alphaPath(t) = alpha;
end

beta_tilde_path = (((1+(h(2)/h(1))*alphaPath).^k).*betaPath);
SimBRN = beta_tilde_path./(gamma+deltaT);
SimERN = (SimData(1:end-1,1)./POP0).* SimBRN;

CumD = SimData(end,4);              % Cumulative deaths during the simulation period
GDPLoss = mean(alphaPath);              % Average output loss during the simulation period
