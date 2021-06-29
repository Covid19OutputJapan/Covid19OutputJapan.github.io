function [V,deltaT,VT,delta_ICU] = ...
    vaccine_distribution(V1_medical,V2_medical,...
    V1_elderly,V2_elderly,V1_others,V2_others,...
    ind_date,ind_date2,date_slowdown,lagged_VT2share,VT3share,...
    elderly,elderly_total,medical,ordinary,delta_average,lag,...
    medical_duration,paces_ori,paces2,paces3,sw_vacpath,gradual_paces,gradual_paces2,...
    E1,E2,D1,D2,ps,POP0,SimPeriod,Tdata)

    % ================================================== %
    % ========= Vaccine Distribution Simulation ======== %
    % == Assuming vaccine distributed accoditng to ps == %
    % ================================================== %
    %Algorithm
    % (1) Calculate the elderly 2nd shots and medical for the first "lag" week
    % (2) Elderly 1st shots are the difference between (1) and avaialbe vaccine
    % (3) After the "lag" periods, distribute Elderly 1st shots so that
    % all elderly will get vaccinated by the given date (ind_date)
    %... should be 3 weeks before the target terminal week fo the second shots
    % (4) Elderly second shots are assinged so that total # of available
    % vaccine shots are fixed.
    % (5) After every elderly get vaccinated, ordinary first shots are assigned
    % so that the total # of available vaccine is constant for the "lag" periods.
    % (6) After "lag" period, Ordinary first shots are distributed with
    % a given share you chose.
    % However, acutal share in the code is set so that the total weeks for
    % ordinary first shots are the same with a given share.
    % (7) Ordinary seconds shots are assinged so that total vaccine shots are fixed.

    
    lag_moderna = lag + 1;
    
    paces_now = V1_elderly(end)+V2_elderly(end)+V1_medical(end)+V2_medical(end)+V1_others(end)+V2_others(end);
    paces = ps*paces_ori; 
    vacpath_ori = zeros(SimPeriod,1);
    if paces_now <= paces
        % vacpath_ori(1+sw_vacpath:gradual_paces) = (paces/(gradual_paces-sw_vacpath)):(paces/(gradual_paces-sw_vacpath)):paces;
        vacpath_ori(1+sw_vacpath:gradual_paces) = paces_now+((paces-paces_now)/(gradual_paces-sw_vacpath)):((paces-paces_now)/(gradual_paces-sw_vacpath)):paces;
        vacpath_ori(gradual_paces+1:end) = paces*ones(SimPeriod-gradual_paces,1);
        % elderly = elderly - (sum(V1_elderly));
    elseif paces_now > paces
        vacpath_ori = paces*ones(SimPeriod,1);
    end

    T = length(vacpath_ori);
    VT = zeros(T,6);

    % (1) Calculate the elderly 2nd shots and medical for the first "lag" week
    VT(1:lag,6) ...
    = VT(1:lag,6)+(sum(V1_medical)-sum(V2_medical) )/lag; % MEDICAL second shots in the first lag periods
    VT(1:medical_duration,5) ...
    =max(0,VT(1:medical_duration,5) ...
    + (medical - sum(V1_medical))/medical_duration); % MEDICAL first shots in the first lag periods
    VT(lag+1:medical_duration+lag,6) ...
    = max(0,(medical-sum(V1_medical))/medical_duration); %distribute the remaining of MEDICAL-second shots
    VT(1:lag,2) = lagged_VT2share.*V1_elderly(end-lag+1:end); % second ELDERLY shots in the first few weeks of simulation
    VT(1:lag_moderna,4) = V1_others(end-lag_moderna+1:end); % second OTHERS shots in the first few weeks of simulation
    sumV = sum(VT,2); % Current total shots at each week after calculating # of vaccination dependent on before-simulation values

    for i = 1:medical_duration + lag
        if sumV(i) > vacpath_ori(i)
            error(['Vaccinated in Data > Simulated Vaccine Path at ' num2str(i) 'th period of simulation'])
        end
    end

    % (2) Elderly 1st shots are the difference between (1) and avaialbe vaccine
    vacpath = max(vacpath_ori- sumV,0); %Compute available slots for vaccination
    VT(1:lag,1) = vacpath(1:lag,1); %First shot for elderly = total shots - (medical shots + second elderly shots)

    % (3) After the "lag" periods over, distribute Elderly 1st shots so that
    % all elderly will get vaccinated by the given date (ind_date)
    remainVT1 = elderly - (sum(VT(:,1))+sum(V1_elderly)); %# of elderly who haven't vaccinated after lag week
    finish_e1 = ind_date-Tdata+1; %the date at which elderly-first shots should end
    if finish_e1 <= lag
        error("lag is larger than finish_e1.")
    end
    VT(lag+1:finish_e1,1) = remainVT1/(finish_e1 - lag); %distribute remaining elderly-first shots
    VT(lag+1:finish_e1,2) = vacpath(lag+1:finish_e1) - VT(lag+1:finish_e1,1); %distribute elderly-second shots until the date of finish_e1
    % (4) Elderly second shots are assinged so that total # of available
    % vaccine shots are fixed.
    remainVT2 = elderly - (sum(VT(:,2))+sum(V2_elderly)); % Number of elderly who does not take second shots at the date of finish_e1
    remainVT2_week = floor(remainVT2/paces)+1; % Number of weeks required for the elderly-second shots
    finish_e2 = finish_e1 + remainVT2_week; %the date at which elderly-second shots will end
    VT(finish_e1+1:finish_e2-1,2) ...
    = vacpath(finish_e1+1:finish_e2-1); %only elderly-second shots during finish_e1 and finish_e2
    VT(finish_e2,2) = remainVT2 - paces*(remainVT2_week-1); %# of second shots left at the date of finish_e2

    % (5) After every elderly get vaccinated, ordinary first shots are assigned
    % so that the total # of available vaccine is constant for the "lag" periods.
    start_o1 = max(1,ind_date2-Tdata);
    enddate = finish_e2+1;
         
    if paces2 > 0
        VT(start_o1:start_o1+gradual_paces2-1,3) = ((paces2)/(gradual_paces2-sw_vacpath)):((paces2)/(gradual_paces2-sw_vacpath)):paces2;
        VT(start_o1+gradual_paces2:enddate,3) = paces2;
    end
    
%     VT(finish_e2,3) = vacpath(finish_e2) - VT(finish_e2,2); %# of ordinary-first shots at the date of finish_e2
    
    
%     VT(finish_e2+1:finish_e2+lag-1,3) = vacpath(finish_e2+1:finish_e2+lag-1); 

    %# of ordinary-first shots after lag weeks since finish_e2
    % VT(finish_e2,3) = vacpath(finish_e2) - VT(finish_e2,2); %# of ordinary-first shots at the date of finish_e2

    % (6) After "lag" period, Ordinary first shots are distributed with
    % a given share you chose.
    % However, acutal share in the code is set so that the total weeks for
    % ordinary first shots are the same with a given share.

    % (7) Ordinary seconds shots are assinged so that total vaccine shots are fixed.
    start_o2 = start_o1 + lag_moderna; % (lag of Moderna) = (lag of Pfizer) + 1   %finish_e2+lag;  %date at which ordinary-second shots start
    vacpath(finish_e2+1:date_slowdown) = vacpath(finish_e2+1:date_slowdown) + paces2;
    vacpath(date_slowdown+1:end) = paces3;
    
    
    %paces_w_shokuiki = paces_ori +  paces2;

    for t = start_o2:finish_e2
        remainVT_others = vacpath(t) - VT(t, 1) - VT(t, 2);
        % remainVT_others = paces2 - VT(t, 3) - VT(t, 4);
        VT(t, 3) = VT(t, 3) + remainVT_others *VT3share;
        VT(t, 4) = VT(t, 4) + remainVT_others *(1-VT3share);
    end
    
    remainVT3 = ordinary - sum(VT(:,3)); %# of finish ordinary-first shots remained
    remainVT3_week  = floor(remainVT3/(vacpath(finish_e2+1)*VT3share))+ 1;
    %remainVT3_week  = floor(remainVT3/(paces*VT3share))+ 1;
    %# of weeks required to finish ordinary-first shots

%     vacpath = max(vacpath_ori- sumV,0); %Compute available slots for vaccination
    
    for t = finish_e2+1:T
        if sum(VT(:,3)) < ordinary % Compute ordinary-first
            VT(t, 3) =  min(vacpath(t)*VT3share, ordinary - sum(VT(:,3)));  %min(remainVT3/remainVT3_week, ordinary - sum(VT(:,3)));
        else
            VT(t, 3) = 0;
        end

        if sum(VT(:,4)) < ordinary % Compute ordinary-second
            VT(t, 4) =  min(vacpath(t) - VT(t,3),ordinary - sum(VT(:,4)));
        else
            break
        end
    end

    % disp('Elderly First Shot Ends at the week of ')
    % disp(date(Tdata + finish_e1))

    % disp('Elderly Second Shot Ends at the week of ')
    % disp(date(Tdata + finish_e2))

    % For debugging purpose
    % disp('Total Vaccinated Elderly for the First Shot')
    % disp(sum(VT(:,1))+sum(V1_elderly))
    %
    % disp('Total Vaccinated Elderly for the Second Shot')
    % disp(sum(VT(:,2))+sum(V2_elderly))
    %
    % disp('# of elderly * accept_share')
    % disp(elderly)
    %
    % disp('Total Vaccinated Medical for the First Shot')
    % disp(sum(VT(:,5))+sum(V1_medical))
    %
    % disp('Total Vaccinated Medical for the Second Shot')
    % disp(sum(VT(:,6))+sum(V2_medical))
    %
    % disp('# of medical')
    % disp(medical)
    %

    V1 = E1*(VT(:,1)+VT(:,3)+VT(:,5))+(E2-E1)*(VT(:,2)+VT(:,4)+VT(:,6));
    V2 = D1*VT(:,1)+(D2-D1)*VT(:,2);
    V_ord = D1*(VT(:,3)+VT(:,5))+(D2-D1)*(VT(:,4)+VT(:,6));
    V2 = [0;0;V2(1:end-2)];
    V = [0;0;V1(1:end-2)];
    V_ord = [0;0;V_ord(1:end-2)];
    
    V2(1) = sum(D1*V1_elderly(1:end-1)+(D2-D1)*V2_elderly(1:end-1));
    V_ord(1) = sum(D1*V1_medical(1:end-1)+(D2-D1)*V2_medical(1:end-1) + D1*V1_others(1:end-1)+(D2-D1)*V2_others(1:end-1));
    V2(2) = D1*V1_elderly(end)+(D2-D1)*V2_elderly(end);
    V_ord(2) = (D1*V1_medical(end)+(D2-D1)*V2_medical(end)) + (D1*V1_others(end)+(D2-D1)*V2_others(end));
    
    V1_prev = E1*(V1_elderly+V1_medical+V1_others)+(E2-E1)*(V2_elderly+V2_medical+V2_others);
    V(1) = V1_prev(end-1);
%     V(1) = V(1) + sum(V1_prev(1:end-1));
    V(2) = V1_prev(end);
    
    delta_ss = delta_average*(0.1063/1.53); 
    % Sheet1 (age_heterogeneity/severe_risk_2021JUN06.xlsx)
    deltaT = (1-(cumsum(V2)/elderly_total))*(delta_average-delta_ss)+(1-(cumsum(V_ord)/(POP0-elderly_total)))*delta_ss;
    
    ICU_ss = delta_average*(0.3916/1.62); 
    % Sheet2 (age_heterogeneity/severe_risk_2021JUN06.xlsx)
%         ICU_ss = delta_ss;
%     ICU_ss = delt a_average*0.1; % Sheet2 (age_heterogeneity/severe_risk_2021JUN06.xlsx)
    delta_ICU = (1-(cumsum(V2)/elderly_total))*(delta_average-ICU_ss)+(1-(cumsum(V_ord)/(POP0-elderly_total)))*ICU_ss;

    if any(VT(:) < 0) % the number of shots exceeds  elderly_total or ordinary_total or medical_total.
        error('Negative Vaccination Value Exists')
    end
