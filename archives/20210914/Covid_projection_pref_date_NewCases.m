function [CumD,GDPLoss,alphapath,SimData,SimN,SimERN,THonPath,SimICU_nation,SimICU_pref,SimHospital,beta] = ...
    Covid_projection_pref_date_NewCases(InitialValues,alpha_on,alpha_off,th_on,th_on2,...
    th_off,th_off2,th_off3,beta,gamma,deltaT,delta_ICU_nation,delta_ICU_pref,delta_Hospital,V,h,k,POP0,...
    hconstant,alpha_duration,state,ICU_nation_inflow_avg, ICU_pref_inflow_avg, Hospital_inflow_avg, ...
    gamma_ICU_nation, gamma_ICU_pref, gamma_Hospital,beta_shock_after_emergency,rho_after_emergency,alpha_jump,th_off_date, ...
    ICU_limit_pref_vec, BED_vec,...
    alpha_lockdown, lockdown_state,...
    lockdown_th_on, lockdown_th_off,...
    beta_jump, beta_goal,...
    SOE_end,...
    betaT_woAR,simple_beta_avg)
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

betaT_ref = beta;

betaT_Path = beta;

beta_additive_path = transpose(0 + (beta_goal) * beta_jump:(1 - beta_jump) * ((beta_goal) / (alpha_duration)):beta_goal);

lockdown_counter = 1;

lockdown_wave = 0;



T = length(beta);
SimData = zeros(T+1,length(InitialValues));
SimData(1,:) = InitialValues;
SimN = zeros(T,1);
SimICU_nation = zeros(T+1,1);
SimICU_pref = zeros(T+1,1);
SimHospital = zeros(T+1,1);
SimICU_nation(1) = InitialValues(5);
SimICU_pref(1) = InitialValues(6);
SimHospital(1) = InitialValues(7);

alpha_soe_path_ori = [alpha_on:(alpha_off - alpha_on)/(alpha_duration ):alpha_off];
alpha_soe_path(1:alpha_duration) = alpha_soe_path_ori(2:alpha_duration+1);
alpha_off_lockdown_path = [alpha_jump:(alpha_off - alpha_jump)/(alpha_duration - 1):alpha_off];

delta_init = deltaT(1);
coef_th_on_path = (th_on2 - th_on)/(1-deltaT(T)/delta_init);
THonPath = th_on + coef_th_on_path.*(1-deltaT./delta_init);

alphapath = alpha_off .* ones(T,1);


for i = 1:T
    
    
    
    if i <= SOE_end
        alpha = alpha_on;
        alphapath(i) = alpha;
    end
    
    alphapath(i) = alpha;
    
    beta(i) = betaT_Path(i);
    
    
    
    if hconstant == 0
        SimN(i) = ((1 + h*alpha)^k)*beta(i)*SimData(i,1)*SimData(i,2)*(1/POP0);
    elseif hconstant == 1
        SimN(i) = ((1+(h(2)/h(1))*alpha)^k)*beta(i)*SimData(i,1)*SimData(i,2)*(1/POP0);
    end
    
    SimData(i+1,1) = SimData(i,1) - SimN(i) - V(i);
    SimData(i+1,2) = SimData(i,2) + SimN(i) - gamma(i)*SimData(i,2) - deltaT(i)*SimData(i,2);
    SimData(i+1,3) = SimData(i,3) + gamma(i)*SimData(i,2) + V(i);
    SimData(i+1,4) = SimData(i,4) + deltaT(i)*SimData(i,2);
    SimICU_nation(i+1) = max(0, SimICU_nation(i) + ICU_nation_inflow_avg*delta_ICU_nation(i)*SimData(i,2) ...
        - gamma_ICU_nation*SimICU_nation(i) - deltaT(i)*SimData(i,2));
    SimICU_pref(i+1) = max(0, SimICU_pref(i) ...
        + ICU_pref_inflow_avg*delta_ICU_pref(i)*SimData(i,2) ...
        - gamma_ICU_pref*SimICU_pref(i) - deltaT(i)*SimData(i,2));
    SimHospital(i+1) = max(0, SimHospital(i) + Hospital_inflow_avg*delta_Hospital(i)*SimData(i,2) ...
        - gamma_Hospital*SimHospital(i));
    

    if i == SOE_end
            betaT_ref(i:i + alpha_duration) = beta(i:i + alpha_duration) + beta_additive_path;
            betaT_ref(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
            beta_goal_ref = betaT_ref(i + alpha_duration) - betaT_ref(i-1);
            betaT_ref_test = betaT_ref(i - 1) + beta_goal_ref * beta_jump;
            beta_jump = (betaT_ref_test - beta(i));
            betaT_Path(i:i + alpha_duration) = beta(i:i + alpha_duration) + beta_additive_path;
            betaT_Path(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
    end
        
    if (SimN(i) >= th_on) && i >= SOE_end
        lockdown_counter = 1;
        lockdown_state = 1;
        alpha = alpha_lockdown;
        betaT_Path(i:end) = simple_beta_avg;%betaT_woAR(i:end);
        
        
    elseif (SimN(i) < th_off) && i >= SOE_end
        
        
        if lockdown_state == 1
            lockdown_wave = lockdown_wave + 1;
            
            betaT_ref(i:i + alpha_duration) = beta(i:i + alpha_duration) + beta_additive_path;
            betaT_ref(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
            beta_goal_ref = betaT_ref(i + alpha_duration) - betaT_ref(i-1);
            betaT_ref_test = betaT_ref(i - 1) + beta_goal_ref * beta_jump;
            beta_jump = (betaT_ref_test - beta(i));
            betaT_Path(i:i + alpha_duration) = beta(i:i + alpha_duration) + beta_additive_path;
            betaT_Path(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
            
            alpha = alpha_lockdown;
        end
        
        lockdown_state = 0;
             
        if (lockdown_counter <= alpha_duration) && lockdown_wave == 0
            alpha = alpha_soe_path(lockdown_counter);
            %                 alphapath(i) = alpha;
            lockdown_counter = lockdown_counter + 1;
            %             beta(i) = betaT(i);
        elseif (lockdown_counter <= alpha_duration) && lockdown_wave >= 1
            alpha = alpha_off_lockdown_path(lockdown_counter);
            %                 alphapath(i) = alpha;
            lockdown_counter = lockdown_counter + 1;
        else
            alpha = alpha_off;
        end
        
    elseif th_off <= SimN(i) && SimN(i) < th_on && i >= SOE_end
        
        
        if lockdown_state == 1
            alpha = alpha_lockdown;

        elseif lockdown_state == 0 && i >= SOE_end
            if (lockdown_counter <= alpha_duration) && lockdown_wave == 0 && i>= SOE_end
                alpha = alpha_soe_path(lockdown_counter);
                lockdown_counter = lockdown_counter + 1;
                
            elseif (lockdown_counter <= alpha_duration) && lockdown_wave >= 1 && i>= SOE_end
                alpha = alpha_off_lockdown_path(lockdown_counter);
                lockdown_counter = lockdown_counter + 1;
                
            else
                
                alpha = alpha_off;

            end
        end
         
    end
    
    
end

SimERN = (SimData(1:end-1,1)./POP0).*(((1+(h(2)/h(1))*alphapath).^k).*beta)./(gamma+deltaT);
CumD = SimData(end,4);              % Cumulative deaths during the simulation period
GDPLoss = mean(alphapath);              % Average output loss during the simulation period
