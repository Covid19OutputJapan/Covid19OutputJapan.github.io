function [CumD,GDPLoss,alphaPath,SimData,SimN,ERN_path_mat,SimICU_nation,SimICU_pref,SimHospital,beta_path_mat,beta_tilde_path_mat,BRN_path_mat,omicron_share] = ...
    Covid_projection_endogenous_omicron(SimData,alpha_on,alpha_off,th_on,...
    th_off,betaT,gamma,deltaT,delta_ICU_nation,delta_ICU_pref,delta_Hospital,V,V_omicron,h,k,POP0,...
    omicron_relative_infectivity, omicron_realtive_severity, omicron_immunity,  ...
    hconstant,DRi,ICU_nation_inflow_avg, ICU_pref_inflow_avg, Hospital_inflow_avg, ...
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
nV = length(SimData(1,1,:));
betaT_Path = betaT;
alpha_Path = alphaBox;
alpha_duration = DRi;
T = length(betaT);
betaPath = zeros(T,1);
alphaPath = betaPath;
omicron_share = betaPath;

SimN = zeros(T,1);
SimN_mat = zeros(T,nV);
beta_path_mat = SimN_mat;
delta_path_mat = SimN_mat;
beta_tilde_path_mat = SimN_mat;
BRN_path_mat = SimN_mat;
ERN_path_mat = SimN_mat;
SimICU_nation = zeros(T+1,1);
SimICU_pref = zeros(T+1,1);
SimHospital = zeros(T+1,1);
POP0_mat = [POP0, POP0-SimData(1,4,1)];


for t = 1:T    
    beta = betaT_Path(t);
    alpha = alpha_Path(t);
    
    if hconstant == 0
        SimN_mat(t,1) = min(((1 + h*alpha)^k)*beta*SimData(t,1,1)*SimData(t,2,1)*(1/POP0), SimData(t,1,1));
        SimN_mat(t,2) = min(omicron_relative_infectivity * ((1 + h*alpha)^k)*beta*SimData(t,1,2)*SimData(t,2,2)*(1/POP0), SimData(t,1,2));
    elseif hconstant == 1
        SimN_mat(t,1) = min(((1+(h(2)/h(1))*alpha)^k)*beta*SimData(t,1,1)*SimData(t,2,1)*(1/POP0), SimData(t,1,1));
        SimN_mat(t,2) = min(omicron_relative_infectivity *((1+(h(2)/h(1))*alpha)^k)*beta*SimData(t,1,2)*SimData(t,2,2)*(1/POP0), SimData(t,1,2)); 
    end
    SimN(t,1) = SimN_mat(t,1) + SimN_mat(t,2);
    weightedSimN = omicron_immunity * SimN_mat(t,1) + SimN_mat(t,2);
%     inR          = gamma(t) * (SimData(t,2,1) + SimData(t,2,2));
%     weightedInR  = omicron_immunity * gamma(t)*SimData(t,2,1) + gamma(t)*SimData(t,2,2);
    omicron_share(t) = SimN_mat(t, 2) / SimN(t,1);
    
    delta_path_mat(t,1) = deltaT(t);
    delta_path_mat(t,2) = omicron_realtive_severity*deltaT(t);
    dD_delta    =  delta_path_mat(t,1)*SimData(t,2,1);
    dD_omicron  = delta_path_mat(t,2)*SimData(t,2,2);

%     SimData(t+1,1,1) = max(0, SimData(t,1,1) - SimN(t,1) - V(t)); %S
    SimData(t+1,1,1) = max(0, SimData(t,1,1) - SimN(t,1) - V(t) - dD_omicron) ; %S
    SimData(t+1,2,1) = SimData(t,2,1) + SimN_mat(t,1) - gamma(t)*SimData(t,2,1) - deltaT(t)*SimData(t,2,1); %I
    SimData(t+1,3,1) = SimData(t,3,1) + gamma(t) * SimData(t,2,1) + V(t) + SimN_mat(t,2); % R
    SimData(t+1,4,1) = SimData(t,4,1) + dD_delta; % death
    SimData(t+1,5,1) = max(0, SimData(t,5,1) ...
                     + ICU_nation_inflow_avg*delta_ICU_nation(t)*SimData(t,2,1) ...
                     - gamma_ICU_nation*SimData(t,5,1) - dD_delta); %ICU nation
    SimData(t+1,6,1) = max(0, SimData(t,6,1) ...
                     + ICU_pref_inflow_avg*delta_ICU_pref(t)*SimData(t,2,1) ...
                     - gamma_ICU_pref*SimData(t,6,1) - dD_delta); % ICU pref
    SimData(t+1,7,1) = max(0, SimData(t,7,1) ...
                     + Hospital_inflow_avg*delta_Hospital(t)*SimData(t,2,1) ...
                     - gamma_Hospital*SimData(t,7,1)); % hospital
                 
    
%     SimData(t+1,1,2) = SimData(t,1,2) - weightedSimN - V_omicron(t); %S
    SimData(t+1,1,2) = SimData(t,1,2) - weightedSimN - V_omicron(t) - dD_delta; %S
    if SimData(t,2,2) > 0
        SimData(t+1,2,2) = SimData(t,2,2) + SimN_mat(t,2) - gamma(t)*SimData(t,2,2) - dD_omicron; %I
    end
    SimData(t+1,3,2) = SimData(t,3,2) + gamma(t)*SimData(t,2,2) + V_omicron(t) + omicron_immunity * SimN_mat(t,1); % R
    SimData(t+1,4,2) = SimData(t,4,2) + dD_omicron; % death
    SimData(t+1,5,2) = max(0, SimData(t,5,2) ...
                     + omicron_realtive_severity*ICU_nation_inflow_avg*delta_ICU_nation(t)*SimData(t,2,2) ...
                     - gamma_ICU_nation*SimData(t,5,2) - dD_omicron); %ICU nation
    SimData(t+1,6,2) = max(0, SimData(t,6,2) ...
                     + omicron_realtive_severity*ICU_pref_inflow_avg*delta_ICU_pref(t)*SimData(t,2,2) ...
                     - gamma_ICU_pref*SimData(t,6,2) - dD_omicron); % ICU pref
    SimData(t+1,7,2) = max(0, SimData(t,7,2) ...
                     + omicron_realtive_severity*Hospital_inflow_avg*delta_Hospital(t)*SimData(t,2,2) ...
                     - gamma_Hospital*SimData(t,7,2)); % hospital
    
    SimICU_nation(t+1,1)   = sum(SimData(t+1,5,:));
    SimICU_pref(t+1,1)     = sum(SimData(t+1,6,:));
    SimHospital(t+1,1)     = sum(SimData(t+1,7,:));
    SimD(t+1,1)            = sum(SimData(t+1,4,:));
        
%     if (SimICU_pref(t+1) >= th_on) %(SimN(t) >= th_on) 
%         state = 1;
%         alpha_Path = alpha_on*ones(T,1);
%         betaT_Path = scale_beta * simple_beta_avg*ones(T,1);
%         
%     elseif (SimICU_pref(t+1) < th_off) %(SimN(t) < th_off)
%         if state == 1
% %             betaT_Path(t + 1:t + alpha_duration + 1) = ...
% %                 linspace(beta_jump*beta_goal(t+1),beta_goal(t + alpha_duration + 1),alpha_duration + 1).*seasonality(t + 1:t + alpha_duration + 1);
%             betaT_Path(t + 1:t + alpha_duration + 1) = ...
%                  linspace(beta_jump*beta_goal(t+1),beta_goal(t + alpha_duration + 1),alpha_duration + 1).*seasonality(t + 1:t + alpha_duration + 1);
%             alpha_Path(t + 1:t + alpha_duration + 1) = ...
%                    [alpha_on:(alpha_off - alpha_on)/(alpha_duration):alpha_off];
%             betaT_Path(t + alpha_duration + 2:end) = betaT(t + alpha_duration + 2:end);
%             alpha_Path(t + alpha_duration + 2:end) = alpha_off;
%         end
%         
%         state = 0;
%     end
    
    
    betaPath(t) = beta * (1 * (1-omicron_share(t)) + omicron_relative_infectivity * omicron_share(t));
    alphaPath(t) = alpha;
    
    beta_path_mat(t,1) = beta;
    beta_path_mat(t,2) = omicron_relative_infectivity * beta;
end
% weightedS = SimData(1:end-1,1,1).*(1-omicron_share) + SimData(1:end-1,1,2).*(omicron_share);
% beta_tilde_path = (((1+(h(2)/h(1))*alphaPath).^k).*betaPath);
% SimBRN = beta_tilde_path./(gamma+deltaT);
% SimERN = (weightedS./POP0).* SimBRN;
for iV = 1:nV
    beta_tilde_path_mat(:,iV) = (((1+(h(2)/h(1))*alphaPath).^k).*beta_path_mat(:,iV));
    BRN_path_mat(:,iV) = beta_tilde_path_mat(:,iV)./(gamma+delta_path_mat(:,iV));
    ERN_path_mat(:,iV) = BRN_path_mat(:,iV) .* (SimData(1:end-1,1,iV) / POP0_mat(iV));
end



CumD = SimData(end,4,1)+SimData(end,4,2);              % Cumulative deaths during the simulation period
GDPLoss = mean(alphaPath);              % Average output loss during the simulation period
