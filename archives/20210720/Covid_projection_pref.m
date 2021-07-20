function [CumD,GDPLoss,alphapath,SimData,SimN,SimERN,THonPath,SimICU_nation,SimICU_pref,SimHospital,beta] = ...
    Covid_projection_pref(InitialValues,alpha_on,alpha_off,th_on1,th_on2,...
    th_off1,th_off2,th_off3,beta,gamma,delta,delta_init,delta_ICU,V,h,k,POP0,...
    hconstant,alpha_duration,state,ICU_inflow_avg_nation,ICU_inflow_avg_pref,gamma_ICU,...
    Hospital_inflow_avg,gamma_Hospital,beta_shock_after_emergency,rho_after_emergency,alpha_jump)
% Simulation for cumulative deaths and output loss given parameters
% InitialValues = [S(end),I(end),R(end),D(end),ICU_nation(end),ICU_pref(end)];
% - InitialValues = 1-by-4 vector of initial values: S, I, R, D, ICU(national definition), ICU(prefecture-specific definitoin)
% - alpha_on is the target alpha with status of emergency
% - beta, gamma, delta and V are T-by-1 vectors of time-varying
% parameters
% - h, k are scalar parameters
% - T is a simulation period
%beta_shock = 1.0;

beta_shock = beta_shock_after_emergency;

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
alphapath = zeros(T,1);
if state == 1 % 1 = state of emergency is on, 0 = it's lifted
    alpha = alpha_on;
else
    alpha = alpha_on+(alpha_off-alpha_on)*alpha_jump;
end
% alpha_off_path = alpha_on:(alpha_off-alpha_on)/(alpha_duration+1):alpha_off;
alpha_off_path = alpha_on+(alpha_off-alpha_on)*alpha_jump:(1-alpha_jump)*((alpha_off-alpha_on)/(alpha_duration+1)):alpha_off;
counter = 1;
wave = 0;
wave_off = 0;
th_on1_ori = th_on1;
th_on2_ori = th_on2;
% th1 = find(delta < delta_ss,1);
% th_on1i = th_on1:(th_on2-th_on1)/(th1-1):th_on2;
% THonPath = [transpose(th_on1i);th_on2.*ones(T-th1,1)];
coef_th_on_path = (th_on2 - th_on1)/(1-delta(T)/delta_init);
THonPath = th_on1 + coef_th_on_path.*(1-delta./delta_init);
for i = 1:T
%     th_on = th_on1 + constant * deltapath(i);
    th_on = THonPath(i);
%     if i < th1
%         th_on = th_on1i(i);
%     else
%         th_on = th_on2;
%     end
    if wave == 0
        th_off = th_off1;
    elseif wave == 1
        th_off = th_off2;
    elseif wave >= 2
        th_off = th_off3;
    end
    alphapath(i) = alpha;
    if hconstant == 0
        SimN(i) = ((1 + h*alpha)^k)*beta(i)*SimData(i,1)*SimData(i,2)*(1/POP0);
    elseif hconstant == 1
        SimN(i) = ((1+(h(2)/h(1))*alpha)^k)*beta(i)*SimData(i,1)*SimData(i,2)*(1/POP0);
    end
    SimData(i+1,1) = SimData(i,1) - SimN(i) - V(i);
    SimData(i+1,2) = SimData(i,2) + SimN(i) - gamma(i)*SimData(i,2) - delta(i)*SimData(i,2);
    SimData(i+1,3) = SimData(i,3) + gamma(i)*SimData(i,2) + V(i);
    SimData(i+1,4) = SimData(i,4) + delta(i)*SimData(i,2);
    SimICU_nation(i+1) = max(0, SimICU_nation(i) + ICU_inflow_avg_nation*delta_ICU(i)*SimData(i,2) ...
        - gamma_ICU*SimICU_nation(i) - delta(i)*SimData(i,2));
    SimICU_pref(i+1) = max(0, SimICU_pref(i) + ICU_inflow_avg_pref*delta_ICU(i)*SimData(i,2) ...
        - gamma_ICU*SimICU_pref(i) - delta(i)*SimData(i,2));
    SimHospital(i+1) = max(0, SimHospital(i) + Hospital_inflow_avg*delta_ICU(i)*SimData(i,2) ...
        - gamma_Hospital*SimHospital(i));
    if th_on <= SimN(i)
        alpha = alpha_on;
        beta_shock = beta_shock_after_emergency;
        counter = 1;
        state = 1;
    elseif th_off >= SimN(i)
        if state == 1
            wave = wave + 1;
            %beta(i+1) = beta(i+1)*(1+beta_shock);
        end
        state = 0;
        if i < T
            beta(i+1) = beta(i+1)*(1+beta_shock);
            %beta(i) = beta(i)*(1+beta_shock);
        end
        beta_shock = beta_shock * rho_after_emergency;
        if counter <= alpha_duration + 1
            alpha = alpha_off_path(counter+1);
            counter = counter + 1;
        else
            alpha = alpha_off;
        end
    elseif th_on > SimN(i) && th_off < SimN(i)
        if state == 1
            alpha = alpha_on;
        elseif state == 0
            if i < T
                beta(i+1) = beta(i+1)*(1+beta_shock);
                %beta(i) = beta(i)*(1+beta_shock);
            end
            beta_shock = beta_shock * rho_after_emergency;
            if counter <= alpha_duration + 1
                alpha = alpha_off_path(counter+1);
                counter = counter + 1;
            else
                alpha = alpha_off;
            end
        end
    end

end
SimERN = (SimData(1:end-1,1)./POP0).*(((1+(h(2)/h(1))*alphapath).^k).*beta)./(gamma+delta);
CumD = SimData(end,4);              % Cumulative deaths during the simulation period
GDPLoss = mean(alphapath);              % Average output loss during the simulation period
