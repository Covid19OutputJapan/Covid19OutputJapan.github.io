function [CumD,GDPLoss,alphapath,SimData,SimN,SimERN,THonPath,SimICU_nation,SimICU_pref,SimHospital,beta] = ...
    Covid_projection_pref_date_lockdown(InitialValues,alpha_on,alpha_off,th_on1,th_on2,...
    th_off1,th_off2,th_off3,beta,gamma,deltaT,delta_ICU_nation,delta_ICU_pref,delta_Hospital,V,h,k,POP0,...
    hconstant,alpha_duration,state,ICU_nation_inflow_avg, ICU_pref_inflow_avg, Hospital_inflow_avg, ...
    gamma_ICU_nation, gamma_ICU_pref, gamma_Hospital,beta_shock_after_emergency,rho_after_emergency,alpha_jump,th_off_date, ...
    ICU_limit_pref_vec, BED_vec,...
    alpha_lockdown, lockdown_state,...
    lockdown_th_on, lockdown_th_off,...
    beta_jump, beta_goal,...
    beta_shock_at_emergency)
% Simulation for cumulative deaths and output loss given parameters
% InitialValues = [S(end),I(end),R(end),D(end),ICU_nation(end),ICU_pref(end)];
% - InitialValues = 1-by-4 vector of initial values: S, I, R, D, ICU(national definition), ICU(prefecture-specific definitoin)
% - alpha_on is the target alpha with status of emergency
% - beta, gamma, delta and V are T-by-1 vectors of time-varying
% parameters
% - h, k are scalar parameters
% - T is a simulation period
%beta_shock = 1.0;

% beta_shock = beta_shock_after_emergency;
beta_counter = 0;

betaT = beta;

betaT_ref = beta;
beta_additive_path = transpose(0 + (beta_goal) * beta_jump:(1 - beta_jump) * ((beta_goal) / (alpha_duration)):beta_goal);

lockdown_counter = 1;

lockdown_wave = 0;

end_counter = 0;




if lockdown_state == 1 % 1 = state of emergency is on, 0 = it's lifted
    alpha = alpha_lockdown;
else
    alpha = alpha_on+(alpha_off-alpha_on)*alpha_jump;
end

% beta_shock = beta_shock_after_emergency;

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

% alpha_off_path = alpha_lockdown+(alpha_off-alpha_on)*alpha_jump+(1-alpha_jump)*((alpha_off-alpha_on)/(alpha_duration+1)):(1-alpha_jump)*((alpha_off-alpha_on)/(alpha_duration+1)):alpha_off;

% alpha_off_path = alpha_lockdown+(alpha_off-alpha_on)*alpha_jump+(1-alpha_jump)*((alpha_off-alpha_on)/(alpha_duration+1)):(1-alpha_jump)*((alpha_off-alpha_on)/(alpha_duration+1)):alpha_off;
alpha_off_lockdown_path = [alpha_jump:(alpha_off - alpha_jump)/alpha_duration:alpha_off];

delta_init = deltaT(1);
coef_th_on_path = (th_on2 - th_on1)/(1-deltaT(T)/delta_init);
THonPath = th_on1 + coef_th_on_path.*(1-deltaT./delta_init);

% alphapath =  alpha_lockdown .* ones(T,1);
% alphapath =  [alpha_on;alpha_lockdown .* ones(T-1,1)];

alphapath = alpha_off .* ones(T,1);
% alphapath(2:alpha_duration) = alpha_off_lockdown_path;
% alphapath(th_off_date:th_off_date+alpha_duration) = alpha_off_path;

% for i = th_off_date:T
%     beta(i) = beta(i)*(1+beta_shock);
%     beta_shock = beta_shock * rho_after_emergency;
% end

for i = 1:T
    %     if i + alpha_duration > 156
    %         end_counter = 1;
    %     end
    alpha = alphapath(i);
    
    %     if i == 1
    %         alpha = alpha_on;
    %     else
    
    %     alphapath(i) = alpha_lockdown;
    %     alpha = alphapath(i);
    
    
    
    if SimICU_nation(i) >= lockdown_th_on *BED_vec(i) || SimICU_pref(i) >= lockdown_th_on * ICU_limit_pref_vec(i)
        lockdown_counter = 1;
        lockdown_state = 1;
        alpha = alpha_lockdown;
        alphapath(i) = alpha;
        
        
        %         beta_shock = 0.1 * i * beta_shock_at_emergency;
        %         if i < T
        %             beta(i) = beta(i)*(1+beta_shock);
        % %             beta(i) = beta(i)*(1+beta_shock);
        %         end
        %         alphapath(i) = alpha;
        
        %         if lockdown_wave == 1 && lockdown_counter > (alpha_duration + 1)
        %             alpha = alpha_off;
        %             alphapath(i) = alpha;
        %         elseif lockdown_counter <= (alpha_duration + 1) && lockdown_wave == 1
        %             alpha = alpha_off_lockdown_path(lockdown_counter);
        %             alphapath(i) = alpha;
        %             lockdown_counter = lockdown_counter + 1;
        %             %         else
        %             %              alpha = alpha_off;
        %             %              alphapath(i) = alpha;
        %             %         elseif lockdown_wave > 1
        %             %             alpha = alpha_off;
        %             %             alphapath(i) = alpha;
        %         end
        
    elseif SimICU_nation(i) < lockdown_th_off * BED_vec(i) && SimICU_pref(i) < lockdown_th_off * ICU_limit_pref_vec(i)
        
        
        if lockdown_state == 1
            lockdown_wave = lockdown_wave + 1;
        end
        
        
        lockdown_state = 0;
        
        if lockdown_state == 0 %&& beta_counter == 0
            if i + alpha_duration <= 208
                betaT_ref(i:i + alpha_duration) = beta(i:i + alpha_duration) + beta_additive_path;
                betaT_ref(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
                beta_goal_ref = betaT_ref(i + alpha_duration) - betaT_ref(i-1);
                betaT_ref_test = betaT_ref(i - 1) + beta_goal_ref * beta_jump;
                beta_jump = (betaT_ref_test - beta(i));
                betaT(i:i + alpha_duration) = beta(i:i + alpha_duration) + beta_additive_path;
                betaT(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
                %                 if i == 1
                %                     betaT_ref(i:i + alpha_duration) = beta(i:i + alpha_duration) + beta_additive_path;
                %                     betaT_ref(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
                %                     beta_goal_ref = betaT_ref(i + alpha_duration) - betaT_ref(1);
                %                     betaT_ref_test = betaT_ref(1) + beta_goal_ref * beta_jump;
                %                     beta_jump = (betaT_ref_test - beta(i));
                %                     betaT(i:i + alpha_duration) = beta(i:i + alpha_duration) + beta_additive_path;
                %                     betaT(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
                %                 else
                %                     betaT_ref(i:i + alpha_duration) = beta(i:i + alpha_duration) + beta_additive_path;
                %                     betaT_ref(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
                %                     beta_goal_ref = betaT_ref(i + alpha_duration) - betaT_ref(i-1);
                %                     betaT_ref_test = betaT_ref(i - 1) + beta_goal_ref * beta_jump;
                %                     beta_jump = (betaT_ref_test - beta(i));
                %                     betaT(i:i + alpha_duration) = beta(i:i + alpha_duration) + beta_additive_path;
                %                     betaT(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
                %                 end
                %             elseif i + alpha_duration > 156 && end_counter == 1
                %
                %
                %
                %                 betaT_ref(i:end) = beta(i:end) + beta_additive_path(1:T - i);
                % %                 betaT_ref(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
                %                 beta_goal_ref = betaT_ref(end) - betaT_ref(i-1);
                %                 betaT_ref_test = betaT_ref(i - 1) + beta_goal_ref * beta_jump;
                %                 beta_jump = (betaT_ref_test - beta(i));
                %                 betaT(i:end) = beta(i:i + alpha_duration) + beta_additive_path;
                %
                %                 end_counter = end_counter + 1;
                % %                 betaT(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
                %             end
            elseif end_counter == 0 && i + alpha_duration > 208
                
                
                add_end = T-i;
                betaT_ref(i:end) = beta(i:end) + beta_additive_path(1:add_end);
                %                 betaT_ref(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
                beta_goal_ref = betaT_ref(end) - betaT_ref(i-1);
                betaT_ref_test = betaT_ref(i - 1) + beta_goal_ref * beta_jump;
                beta_jump = (betaT_ref_test - beta(i));
                betaT(i:end) = beta(i:i + alpha_duration) + beta_additive_path;
                
                
                %                 betaT(i + alpha_duration + 1:end) = beta(i + alpha_duration + 1:end) + beta_additive_path(end);
            end
            
            end_counter = end_counter + 1;
            
            beta_counter = beta_counter + 1;
        end
        %         if i < T
        %             beta(i) = beta(i)*(1+beta_shock);
        %             %beta(i) = beta(i)*(1+beta_shock);
        %         end
        %         beta_shock = beta_shock * rho_after_emergency;
        %         alpha = alpha_off_lockdown_path(1);
        %         alphapath(i) = alpha;
        %         lockdown_counter = 1;
        %         lockdown_state = 0;
        if lockdown_counter <= alpha_duration + 1
            alpha = alpha_off_lockdown_path(lockdown_counter);
            alphapath(i) = alpha;
            lockdown_counter = lockdown_counter + 1;
            %             beta(i) = betaT(i);
        elseif lockdown_counter > alpha_duration + 1
            alpha = alpha_off;
            alphapath(i) = alpha;
            %             beta(i) = betaT(i);
        end
        
    elseif ((lockdown_th_off * BED_vec(i) <= SimICU_nation(i)) && (SimICU_nation(i) < lockdown_th_on * BED_vec(i))) || ((lockdown_th_off * ICU_limit_pref_vec(i) <= SimICU_pref(i)) && (SimICU_pref(i) < lockdown_th_on * ICU_limit_pref_vec(i)))
        
        %         if lockdown_state == 1
        %             alpha = alpha_lockdown;
        %         end
        %
        %         if lockdown_counter <= alpha_duration + 1
        %             alpha = alpha_off_lockdown_path(lockdown_counter);
        %             alphapath(i) = alpha;
        %             lockdown_counter = lockdown_counter + 1;
        %             %             beta(i) = betaT(i);
        %         else
        %             alpha = alpha_off;
        %             alphapath(i) = alpha;
        %             %             beta(i) = betaT(i);
        %         end
        if lockdown_state == 1
            alpha = alpha_lockdown;
            alphapath(i) = alpha;
        elseif lockdown_state == 0
            if lockdown_counter <= alpha_duration + 1
                alpha = alpha_off_lockdown_path(lockdown_counter);
                alphapath(i) = alpha;
                lockdown_counter = lockdown_counter + 1;
                %             beta(i) = betaT(i);
            else
                alpha = alpha_off;
                alphapath(i) = alpha;
                %             beta(i) = betaT(i);
            end
        end
        
    end
    
    
    
    
    
    beta(i) = betaT(i);
    
    if i <= 3
        alpha = alpha_on;
        alphapath(i) = alpha;
    end
    
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
    
    
    
end
SimERN = (SimData(1:end-1,1)./POP0).*(((1+(h(2)/h(1))*alphapath).^k).*beta)./(gamma+deltaT);
CumD = SimData(end,4);              % Cumulative deaths during the simulation period
GDPLoss = mean(alphapath);              % Average output loss during the simulation period
