function [V1_medical, V1_medical2, V2_medical, V2_medical2, ...
    V1_elderly, V1_elderly2, V2_elderly, V2_elderly2,...
    V1_others, V1_others2, V2_others, V2_others2,...
    vs_MT,vs_ET,vs_M1,vs_E1,vs_M2,vs_E2] ...
    = ImportVaccineData(home,iPC,Data,pref,dateEN,ps,vaccine_disp_switch)
%%
% if iPC==1
%     covid = importdata([home '\Covid_weekly_newV.csv']);  % Import weekly Covid data by prefecture
% else
%     covid = importdata([home 'Covid_weekly_newV.csv']);  % Import weekly Covid data by prefecture
% end
% Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
Tdata= size(Data,1);
E_total = cum_to_new(Data(:,8));
E_first = cum_to_new(Data(:,9));
E_second = cum_to_new(Data(:,10));
M_total = cum_to_new(Data(:,11));
M_first = cum_to_new(Data(:,12));
M_second = cum_to_new(Data(:,13));
% % ====================================
% O_total = cum_to_new(Data(:,11));
% O_first = cum_to_new(Data(:,12));
% O_second = cum_to_new(Data(:,13));
% % ====================================
M_ps_total = Data(:,14);
M_ps_first = Data(:,15);
M_ps_second = Data(:,16);
E_ps_total = Data(:,17);
E_ps_first = Data(:,18);
E_ps_second = Data(:,19);
% % ====================================
% O_ps_total = Data(:,17);
% O_ps_first = Data(:,18);
% O_ps_second = Data(:,19);
% % ====================================

vs_MT = M_ps_total(end);
vs_M1 = M_ps_first(end);
vs_M2 = M_ps_second(end);
vs_ET = E_ps_total(end);
vs_E1 = E_ps_first(end);
vs_E2 = E_ps_second(end);
% % ====================================
% vs_OT = O_ps_total(end);
% vs_O1 = O_ps_first(end);
% vs_O2 = O_ps_second(end);
% % ====================================


vaccine_medical = readmatrix([home 'vaccine_daily_medical.xls']);
vaccine_elderly = readmatrix([home 'vaccine_daily_elderly.xls']);
% ====================================
vaccine_others = readmatrix([home 'vaccine_daily_others.xls']);
% ====================================

[V1_medical, V2_medical] = vaccine_daily_to_weekly_table(vaccine_medical, ps, dateEN,iPC);
[V1_elderly, V2_elderly] = vaccine_daily_to_weekly_table(vaccine_elderly, ps, dateEN,iPC);
[V1_others, V2_others] = vaccine_daily_to_weekly_table(vaccine_others, ps, dateEN,iPC);

% Calculate the ratio of the number of vaccines for elderly and others
V1_elderly_share = V1_elderly./(V1_elderly + V1_others);
V1_elderly_share(isnan(V1_elderly_share)) = 0;
V2_elderly_share = V2_elderly./(V2_elderly + V2_others);
V2_elderly_share(isnan(V2_elderly_share)) = 0;

indM1 = find(M_ps_first>0,1,'first');
indM2 = find(M_ps_first>0,1,'last');
indM21 = find(M_ps_second>0,1,'first');
indM22 = find(M_ps_second>0,1,'last');

V1_medical2 = zeros(Tdata,1);
V1_medical2(1:length(V1_medical)) = (V1_medical/ps)*M_ps_first(indM1);
V1_medical2(indM1+1:end) = M_first(indM1+1:end);

V2_medical2 = zeros(Tdata,1);
V2_medical2(1:length(V2_medical)) = (V2_medical/ps)*M_ps_second(indM1);
V2_medical2(indM21+1:end) = M_second(indM21+1:end);
if V1_medical(Tdata) == 0
    V1_medical(Tdata) = V1_medical(Tdata-1);
    V2_medical(Tdata) = V2_medical(Tdata-1);
end
if V1_medical2(Tdata) == 0
    V1_medical2(Tdata) = V1_medical2(Tdata-1);
    V2_medical2(Tdata) = V2_medical2(Tdata-1);
end


% V1_elderly2 = E_first;
% V2_elderly2 = E_second;
V1_elderly2 = E_first.*V1_elderly_share;
V2_elderly2 = E_second.*V2_elderly_share;
V1_others2 = E_first - V1_elderly2;
V2_others2 = E_second - V2_elderly2;

if V1_elderly(Tdata) == 0
    V1_elderly(Tdata) = V1_elderly(Tdata-1);
    V2_elderly(Tdata) = V2_elderly(Tdata-1);
end
if V1_elderly2(Tdata) == 0
    V1_elderly2(Tdata) = V1_elderly2(Tdata-1);
    V2_elderly2(Tdata) = V2_elderly2(Tdata-1);
end



indM = find(V1_medical>0,1,'first');
VacMat_M = [V1_medical, V1_medical2, V2_medical, V2_medical2];
indE = find(V1_elderly>0,1,'first');
VacMat_E = [V1_elderly, V1_elderly2, V2_elderly, V2_elderly2];
title = ["First shots (pop share)", "First shot", "Second shots (pop share)", "Second shot"];
cum_E_first = sum(V1_elderly2);
cum_E_second = sum(V2_elderly2);
cum_M_first = sum(V1_medical2);
cum_M_second = sum(V2_medical2);
if vaccine_disp_switch == 1
    disp("Name of Prefecture:")
    disp(pref)
    disp("Since the week of ")
    disp(dateEN(indM))
    disp("Vaccinated Madical Personnels")
    disp(title)
    disp(VacMat_M(indM:end,:))
    
    
    disp("Since the week of ")
    disp(dateEN(indE))
    disp("Vaccinated Elderly")
    disp(title)
    disp(VacMat_E(indE:end,:))
    
    
    
    disp("Latest # of Cumulative First Shots, Medical ")
    disp(cum_M_first)
    disp("Latest # of Cumulative Second Shots, Medical ")
    disp(cum_M_second)
    disp("Latest # of Cumulative Total Shots, Medical ")
    disp(cum_M_first+cum_M_second)
    disp("Latest # of Cumulative First Shots, Elderly ")
    disp(cum_E_first)
    disp("Latest # of Cumulative Second Shots, Elderly ")
    disp(cum_E_second)
    disp("Latest # of Cumulative Total Shots, Elderly ")
    disp(cum_E_first+cum_E_second)
end


% Total in the week of 4/23 for Medical: 726496 (177672+150489+152884+126419+119032)
% First in the week of 4/23 for Medical: 726496 (147989+128634+126188+98622+64725)
% sumt = sum(M_total(indM1:indM2),1);
% sumf = sum(M_first(indM1:indM2),1);
% sums = sum(M_second(indM1:indM2),1);

%    disp("Prefecture:")
%    disp(pref)
%    disp(title)
%    disp(compare(pindex,:))
