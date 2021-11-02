
S_elderly_prev = elderly_total - sum(E1*V1_elderly(1:end-2) + (E2 - E1)*V2_elderly(1:end-2));
V1_elderly_Sim = [V1_elderly(end-1);V1_elderly(end);VT(1:end-2,5)];
V2_elderly_Sim = [V2_elderly(end-1);V2_elderly(end);VT(1:end-2,6)];
S_elderly_path = S_elderly_prev - cumsum(E1 * V1_elderly_Sim+ (E2 - E1) * V2_elderly_Sim);

S_young_prev =  ordinary_total + medical_total...
    - sum(E1*(V1_others(1:end-2) + V1_medical(1:end-2))) ...
    - sum((E2 - E1)*(V2_others(1:end-2) + V2_medical(1:end-2)));
V1_others_Sim = [V1_others(end-1);V1_others(end);VT(1:end-2,3)];
V2_others_Sim = [V2_others(end-1);V2_others(end);VT(1:end-2,4)];
V1_medical_Sim = [V1_medical(end-1);V1_medical(end);VT(1:end-2,1)];
V2_medical_Sim = [V2_medical(end-1);V2_medical(end);VT(1:end-2,2)];
V1_young_Sim = V1_others_Sim + V1_medical_Sim;
V2_young_Sim = V2_others_Sim + V2_medical_Sim;
S_young_path = S_young_prev - cumsum(E1 * V1_young_Sim+ (E2 - E1) * V2_young_Sim);
elderly_share_Sim = S_elderly_path./(S_young_path+S_elderly_path);
elderly_share_Tdata = S_elderly_prev(end)/(S_young_prev(end)+S_elderly_prev(end));

deltaT = elderly_share_Sim./elderly_share_Tdata * delta_average .* lambda_delta;

%     delta_ss = delta_average*(0.1063/1.53); % Sheet1 (age_heterogeneity/severe_risk_2021JUN06.xlsx)
%     deltaT = (1-(cumsum(V2)/elderly_total))*(delta_average-delta_ss)+(1-(cumsum(V_ord)/(POP0-elderly_total)))*delta_ss;

%     delta_ICU_nation_ss = delta_ICU_nation_average*(0.3916/1.62); % Sheet2 (age_heterogeneity/severe_risk_2021JUN06.xlsx)
%     delta_ICU_nation = (1-(cumsum(V2)/elderly_total))*(delta_ICU_nation_average-delta_ICU_nation_ss)+(1-(cumsum(V_ord)/(POP0-elderly_total)))*delta_ICU_nation_ss;
delta_ICU_nation = elderly_share_Sim./elderly_share_Tdata * delta_ICU_nation_average .* lambda_delta_ICU_nation;

%     delta_ICU_pref_ss = delta_ICU_pref_average*(0.3916/1.62); % Sheet2 (age_heterogeneity/severe_risk_2021JUN06.xlsx)
%     delta_ICU_pref = (1-(cumsum(V2)/elderly_total))*(delta_ICU_pref_average-delta_ICU_pref_ss)+(1-(cumsum(V_ord)/(POP0-elderly_total)))*delta_ICU_pref_ss;
delta_ICU_pref = elderly_share_Sim./elderly_share_Tdata * delta_ICU_pref_average .* lambda_delta_ICU_pref;

%     delta_Hospital_ss = delta_Hospital_average*(0.3916/1.62); % Sheet2 (age_heterogeneity/severe_risk_2021JUN06.xlsx)
%     delta_Hospital = (1-(cumsum(V2)/elderly_total))*(delta_Hospital_average-delta_Hospital_ss)+(1-(cumsum(V_ord)/(POP0-elderly_total)))*delta_Hospital_ss;
delta_Hospital = elderly_share_Sim./elderly_share_Tdata * delta_Hospital_average .* lambda_delta_Hospital;
