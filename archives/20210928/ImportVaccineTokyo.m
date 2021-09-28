function [V1_medical_w,V2_medical_w,V1_elderly_w,V2_elderly_w,V1_others_w,V2_others_w]...
    = ImportVaccineTokyo(home,dateEN,Data,Tdata,ps,iPC)
%%

vaccine_non_medical = readmatrix([home 'vaccination_Tokyo_cumulative.xlsx']);
dateV_ori = vaccine_non_medical(:,1) + 21916;
dateV = datetime(dateV_ori,'ConvertFrom','excel');
V1_elderly_d = cum_to_new(vaccine_non_medical(:,3));
V2_elderly_d = cum_to_new(vaccine_non_medical(:,4));
V1_others_d = cum_to_new(vaccine_non_medical(:,6));
V2_others_d = cum_to_new(vaccine_non_medical(:,7));

%Non-medical personel (age 65 or older) ... elderly
V1_elderly_d(isnan(V1_elderly_d)) = 0;
V2_elderly_d(isnan(V2_elderly_d)) = 0;
Vdata = size(V1_elderly_d,1);
V1_elderly_w = zeros(length(dateEN),1);
V2_elderly_w = zeros(length(dateEN),1);
for i = 1:1:length(dateEN)
    index_date = (dateV >= dateEN(i)-3 & dateV <= dateEN(i)+3);
    V1_elderly_w(i) = sum(V1_elderly_d(index_date));
    V2_elderly_w(i) = sum(V2_elderly_d(index_date));
end

%Non-medical personel (age 15-64)
V1_others_d(isnan(V1_others_d)) = 0;
V2_others_d(isnan(V2_others_d)) = 0;
Vdata = size(V1_others_d,1);
V1_others_w = zeros(length(dateEN),1);
V2_others_w = zeros(length(dateEN),1);
for i = 1:1:length(dateEN)
    index_date = (dateV >= dateEN(i)-3 & dateV <= dateEN(i)+3);
    V1_others_w(i) = sum(V1_others_d(index_date));
    V2_others_w(i) = sum(V2_others_d(index_date));
end

% Medical personels
vaccine_medical = readmatrix([home 'vaccine_daily_medical.xls']);
[V1_medical_past, V2_medical_past] = vaccine_daily_to_weekly_table(vaccine_medical, ps, dateEN,iPC);

M_first = cum_to_new(Data(:,12));
M_second = cum_to_new(Data(:,13));
M_ps_first = Data(:,15);
M_ps_second = Data(:,16);
indM1 = find(M_ps_first>0,1,'first');
indM21 = find(M_ps_second>0,1,'first');

V1_medical_w = zeros(Tdata,1);
V1_medical_w(1:length(V1_medical_w)) = (V1_medical_past/ps)*M_ps_first(indM1);
V1_medical_w(indM1+1:end) = M_first(indM1+1:end);

V2_medical_w = zeros(Tdata,1);
V2_medical_w(1:length(V2_medical_w)) = (V2_medical_past/ps)*M_ps_second(indM1);
V2_medical_w(indM21+1:end) = M_second(indM21+1:end);
if V1_medical_w(Tdata) == 0
    V1_medical_w(Tdata) = V1_medical_w(Tdata-1);
    V2_medical_w(Tdata) = V2_medical_w(Tdata-1);
end
if V1_medical_w(Tdata) == 0
    V1_medical_w(Tdata) = V1_medical_w(Tdata-1);
    V2_medical_w(Tdata) = V2_medical_w(Tdata-1);
end
