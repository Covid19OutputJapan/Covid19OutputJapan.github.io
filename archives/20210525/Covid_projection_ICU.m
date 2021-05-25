function [CumD,GDPLoss,alphapath,SimData,SimN,SimERN,ICU] = Covid_projection_ICU(InitialValues,alpha_on,alpha_off,th_on,th_off1,th_off2,th_off3,beta,gamma,delta,V,h,k,POP0,hconstant,alpha_duration,current_alpha,state,ICU_death,gamma_ICU)
% Simulation for cumulative deaths and output loss given parameters

% - InitialValues = 1-by-5 vector of initial values: S, I, R, D,ICU
% - alpha_on is the target alpha with status of emergency
% - beta, gamma, delta and V are T-by-1 vectors of time-varying
% parameters
% - h, k are scalar parameters
% - T is a simulation period
% - state = 1 --> state of emergency is on, 0 otherwise

T = length(beta);
SimData = zeros(T+1,length(InitialValues));
SimData(1,:) = InitialValues;
SimN = zeros(T,1);
ICU = zeros(T+1,1);
ICU(1) = InitialValues(5);
alphapath = zeros(T,1);
%state = SOE;   % 1 = state of emergency is on, 0 = it's lifted
alpha_off_path = alpha_on:(alpha_off-alpha_on)/(alpha_duration+1):alpha_off;
alpha_off_path_now = current_alpha:(alpha_off-current_alpha)/(alpha_duration+1):alpha_off;
counter = 1;
wave = 0;

if state == 0
    alpha = current_alpha;
elseif state == 1
    alpha = alpha_on;
end


for i = 1:T
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
    ICU(i+1) = ICU(i) + ICU_death*delta(i)*SimData(i,2) - gamma_ICU*ICU(i);
    if th_on <= ICU(i)
        alpha = alpha_on;
        counter = 1;
        state = 1;
    elseif th_off >= ICU(i)
        if state == 1
            wave = wave + 1;
        end
        state = 0;
        if counter <= alpha_duration + 1
            if wave == 0 % when your current N is smaller than th_off1
                alpha = alpha_off_path_now(counter+1);
            else % N path after lifting the state of emergency
                alpha = alpha_off_path(counter+1);
            end
            counter = counter + 1;
        else
            alpha = alpha_off;
        end
    elseif th_on > ICU(i) && th_off < ICU(i)
        if state == 1
            alpha = alpha_on;
        elseif state == 0
            if counter <= alpha_duration + 1
                if wave == 0 % when your current N is smaller than th_off1
                    alpha = alpha_off_path_now(counter+1);
                else % N path after lifting the state of emergency
                    alpha = alpha_off_path(counter+1);
                end
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
