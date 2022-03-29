function [CumD,GDPLoss,alphaPath,SimData,SimN,SimERN,SimICU_nation,SimICU_pref,SimHospital,SimNewICU_pref,betaPath,beta_tilde_path,SimBRN] = ...
    Covid_projection(InitialValues,betaT,gamma,deltaT,...
                     SimICU_nation_rate, SimICU_pref_rate, ...
                     SimHospital_rate,SimNewICU_pref_rate, ...
                     h,k,POP0, V, E1, E2, E3, ...
                     hconstant,DRi, gamma_ICU_nation, gamma_ICU_pref, ...
                     gamma_Hospital,gamma_newICU_pref, alphaBox)
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
SimNewICU_pref = zeros(T+1,1);
SimHospital = zeros(T+1,1);
SimICU_nation(1) = InitialValues(5);
SimICU_pref(1) = InitialValues(6);
SimHospital(1) = InitialValues(7);
SimNewICU_pref(1) = InitialValues(8);

for t = 1:T    
    beta = betaT_Path(t);
    alpha = alpha_Path(t);
    if hconstant == 0
        SimN(t) = ((1 + h*alpha)^k)*beta*SimData(t,1)*SimData(t,2)*(1/POP0);
    elseif hconstant == 1
        SimN(t) = ((1+(h(2)/h(1))*alpha)^k)*beta*SimData(t,1)*SimData(t,2)*(1/POP0);
    end
    SimData(t+1,1) = SimData(t,1) - SimN(t) - V(t); %S
    SimData(t+1,2) = SimData(t,2) + SimN(t) - gamma(t)*SimData(t,2) - deltaT(t)*SimData(t,2); %I
    SimData(t+1,3) = SimData(t,3) + gamma(t)*SimData(t,2) + V(t); %R
    SimData(t+1,4) = SimData(t,4) + deltaT(t)*SimData(t,2); %D
    
    SimICU_nation(t+1) = max(0, SimICU_nation(t) + SimICU_nation_rate(t)*SimData(t,2) ...
        - gamma_ICU_nation(t)*SimICU_nation(t) - deltaT(t)*SimData(t,2));
    SimICU_pref(t+1) = max(0, SimICU_pref(t) ...
        + SimICU_pref_rate(t)*SimData(t,2) ...
        - gamma_ICU_pref(t)*SimICU_pref(t) - deltaT(t)*SimData(t,2));
    SimNewICU_pref(t+1) = max(0, SimNewICU_pref(t) ...
        + SimNewICU_pref_rate(t)*SimData(t,2) ...
        - gamma_newICU_pref(t)*SimNewICU_pref(t) - deltaT(t)*SimData(t,2));
    SimHospital(t+1) = max(0, SimHospital(t) + SimHospital_rate(t)*SimData(t,2) ...
        - gamma_Hospital(t)*SimHospital(t));
    
%         
%     if t == SoE_date %(SimNewICU_pref(t+1) >= th_on && t >=2) &&  (SimN(t) >= th_on_N)  %
%         state = 1;
%         alpha_Path = alpha_on*ones(T,1);
%         betaT_Path = scale_beta * simple_beta_avg*ones(T,1);
%         
%     elseif (t >= SoE_off_date)  %(SimNewICU_pref(t+1) < th_off) || (SimN(t) < th_off_N)   %(SimN(t) < th_off)
%         if state == 1
%             if t + alpha_duration + 1 > T
%                 remainT = T-t;
%                 nextBetaPath(1:alpha_duration+1) = linspace(beta_jump*beta_goal,beta_goal,alpha_duration + 1).*seasonality(t + 1:t + alpha_duration + 1);
%                 betaT_Path(t + 1:end) = nextBetaPath(1:remainT);
% %                 nextAlphaPath(1:alpha_duration+1) = alpha_on:(alpha_off - alpha_on)/(alpha_duration):alpha_off;
% %                 alpha_Path(t + 1:end) = nextAlphaPath(1:remainT);
%             else
%                 betaT_Path(t + 1:t + alpha_duration + 1) = ...
%                     linspace(beta_jump*beta_goal,beta_goal,alpha_duration + 1).*seasonality(t + 1:t + alpha_duration + 1);
% %                 alpha_Path(t + 1:t + alpha_duration + 1) = ...
% %                    alpha_on:(alpha_off - alpha_on)/(alpha_duration):alpha_off;
%                 betaT_Path(t + alpha_duration + 2:end) = betaT(t + alpha_duration + 2:end);
% %                 alpha_Path(t + alpha_duration + 2:end) = alpha_off;
%             end
%         end
%     elseif (t <= SoE_off_date)
%         if t- t_SoE_lag <= 0
%             beta = (sum(past_beta(end +t -t_SoE_lag : end)) + sum(betaT_Path(1 : t-1))) / t_SoE_lag;
%         else
%             beta = mean(betaT_Path(t-t_SoE_lag:t-1));
%         end
%     
% 
%         state = 0;
%     end
        
    betaPath(t) = beta;
    alphaPath(t) = alpha;
end

beta_tilde_path = (((1+(h(2)/h(1))*alphaPath).^k).*betaPath);
SimBRN = beta_tilde_path./(gamma+deltaT);
SimERN = (SimData(1:end-1,1)./POP0).* SimBRN;

CumD = SimData(end,4);              % Cumulative deaths during the simulation period
GDPLoss = mean(alphaPath);              % Average output loss during the simulation period
