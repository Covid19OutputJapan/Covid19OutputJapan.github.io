function [V,deltaT,VT] = ...
    vaccine_distribution_medical(V1_medical,V2_medical,...
    V1_elderly,V2_elderly,elderly_jp,medical_jp,ordinary_jp,accept_share,...
    delta_average,lag,medical_duration,paces_ori,sw_vacpath,gradual_paces,...
    E1,E2,D1,D2,ps,POP0,SimPeriod)
% Compute the path of vaccinations for each group
paces = ps*paces_ori; %3600000;
vacpath = zeros(SimPeriod,1);
vacpath(1+sw_vacpath:gradual_paces) = (paces/(gradual_paces-sw_vacpath)):(paces/(gradual_paces-sw_vacpath)):paces;
vacpath(gradual_paces+1:end) = paces*ones(SimPeriod-gradual_paces,1);
elderly_total = ps*elderly_jp;
medical_total = ps*medical_jp;
ordinary_total = ps*ordinary_jp;
%medical = medical_total*accept_share;
medical = medical_total;
elderly = elderly_total*accept_share;
ordinary = ordinary_total*accept_share;
elderly = elderly - (sum(V1_elderly));


T = length(vacpath);
VT = zeros(T,6);
for t = 1:T
    CV = cumsum(VT);
    if t <= lag
        if CV(t,1) + vacpath(t) <= elderly
            VT(t,1) = vacpath(t);
        elseif CV(t,1) + vacpath(t) > elderly && CV(t,1) < elderly
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1);
        elseif CV(t,1) >= elderly
            if CV(t,3) + vacpath(t) <= ordinary
                VT(t,3) = vacpath(t);
            elseif CV(t,3) + vacpath(t) > ordinary && CV(t,3) < ordinary
                VT(t,3) = ordinary - CV(t,3);
            end
        end
    elseif CV(t,1) + vacpath(t) <= elderly
        if VT(t-lag,1) == 0
            VT(t,1) = vacpath(t);
        elseif VT(t-lag,1) > 0
            VT(t,2) = VT(t-lag,1);
            VT(t,1) = vacpath(t) - VT(t,2);
        end
    elseif CV(t,1) + vacpath(t) > elderly && CV(t,1) + vacpath(t) - VT(t-lag,1) <= elderly
        if VT(t-lag,1) == 0
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1);
        elseif VT(t-lag,1) > 0
            VT(t,2) = VT(t-lag,1);
            VT(t,1) = vacpath(t) - VT(t,2);
        end
    elseif CV(t,1) + vacpath(t) - VT(t-lag,1) >= elderly && CV(t,1) < elderly
        if VT(t-lag,1) == 0
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1);
        elseif VT(t-lag,1) > 0
            VT(t,2) = VT(t-lag,1);
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1) - VT(t,2);
        end
    elseif CV(t,1) >= elderly
        if VT(t-lag,1) > 0
            if VT(t-lag,3) == 0
                VT(t,2) = VT(t-lag,1);
                VT(t,3) = vacpath(t) - VT(t,2);
            else
                VT(t,2) = VT(t-lag,1);
                VT(t,4) = VT(t-lag,3);
                VT(t,3) = vacpath(t) - VT(t,2) - VT(t,4);
            end
        elseif CV(t,3) < ordinary
            if VT(t-lag,3) > 0
                VT(t,4) = VT(t-lag,3);
                VT(t,3) = vacpath(t) - VT(t,4);
            else
                VT(t,3) = vacpath(t);
            end
        elseif CV(t,3) >= ordinary
            if VT(t-lag,3) > 0
                VT(t,4) = VT(t-lag,3);
            end
        end
    end
end

VT(1:lag,6) = VT(1:lag,6)+(sum(V1_medical)-sum(V2_medical) )/lag;
VT(1:medical_duration,5) =VT(1:medical_duration,5) + (medical - sum(V1_medical))/medical_duration;
VT(1+lag:medical_duration+lag,6) = VT(1+lag:medical_duration+lag,6)+  (medical-sum(V1_medical))/medical_duration;
VT(1:lag,2) = VT(1:lag,2)+(sum(V1_elderly)-sum(V2_elderly) )/lag;
V1 = E1*(VT(:,1)+VT(:,3)+VT(:,5))+(E2-E1)*(VT(:,2)+VT(:,4)+VT(:,6));
V2 = D1*VT(:,1)+(D2-D1)*VT(:,2);
V2 = [0;0;V2(1:end-2)];
V = [0;0;V1(1:end-2)];
% V_ord = D1*(VT(:,3)+VT(:,5))+D2*(VT(:,4)+VT(:,6));
V_ord = D1*(VT(:,3)+VT(:,5))+(D2-D1)*(VT(:,4)+VT(:,6));
V_ord = [0;0;V_ord(1:end-2)];

delta_ss = delta_average*(0.1063/1.53);
% delta_ss = delta_average*(0.09/1.28);

deltaT = (1-(cumsum(V2)/elderly_total))*(delta_average-delta_ss)+(1-(cumsum(V_ord)/(POP0-elderly_total)))*delta_ss;
