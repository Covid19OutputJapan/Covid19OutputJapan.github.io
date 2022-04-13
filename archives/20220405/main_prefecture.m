% This m-file executes simulation and generates figures for the  analysis of 
% the spread of Covid-19 in Tokyo shown at https://covid19outputjapan.github.io
clear all;  close all;
tic;
if ispc == 1
    HOME        = 'C:\Users\畝矢\Dropbox\PC\Documents\GitHub\github_Codes\';
    DATA_PATH   = 'C:\Users\畝矢\Dropbox\fujii_nakata\Website\Codes\Data\';
    FIGURE_PATH = 'C:\Users\畝矢\Dropbox\fujii_nakata\Website\Codes\Figures\';
    FONT_NAME          = 'Yu Gothic'; % Font style for xaxis, yaxis, title
elseif ismac == 1
    HOME        = '/Users/ymaeda/Desktop/GitHub/github_Codes/'; %include \ or / at the end
    DATA_PATH   = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/Data/';
    FIGURE_PATH = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/Figures/';
    FONT_NAME   = 'YuGothic';
end
cd(HOME);
%====================== Program parameter values ======================%
PREF = 'Tokyo';
PREF_GDP = 106; %Cho-yen; Trillion yen
FIGURE_SAVE = 1;% 0 = figures won't be saved, 1 = they will be saved
DATA_SAVE = 1; % save back data
DATA_SWITCH = 1; % Use I_data and gamma = mean(gamma_data(end-17+1:end))
FONT_SIZE = 12; % common font size for many figures
LEGEND_FONT_SIZE = 6;
LEGEND_FONT_SIZE_MAIN = 12;
AXIS_FONT_SIZE = 8;
ft = '%.1f';
Y_FORMAT = '%.0f';
LANGUAGE = {'EN', 'JP'};

%===================== Figure Names =================%
figname_main = string(['MainResults']);
figname_beta_tilde = 'Beta Tilde Path';
figname_beta = 'Beta Path';
figname_var = 'Variant_Share';
figname_delta = 'Death Rate transition';
figname_ICU_nation = 'Severity Rate transition National Standard';
figname_ICU_local = 'Severity Rate transition Local Standard';
figname_new_ICU_local = 'Severity Rate transition New Local Standard';
figname_ERN = 'ERN transition';
figname_Hospital = 'Hospital Rate transition';
figname_dD = 'dD transition';
figname_BRN = 'BRN';
figname_omicron_share = 'Omicron share';

%================== Model Fixed Parameter Values ============================%
parameter
SimPeriod = 52;

%=============== Import data ============%
import_prefecture


newICU_pref     = zeros(Tdata+1, 1); %Read # of severe cases under the new definition in Tokyo
NewSevere       = readmatrix([DATA_PATH 'tokyo_new_severe_cases.csv']);
initialICUdata  = datetime(2020,12,17);
indNewICU       = find(date == initialICUdata,1,'first');
newICU_pref(indNewICU+1:end,1) = NewSevere(:,6);

% sample_period = find(date == datetime(2021, 6, 17)):find(date == datetime(2021, 11, 4)); %The periods of the 5th wave
if find(date <= datetime(2022, 3, 31),1,'last') < Tdata
    sample_period = find(date >= datetime(2022, 2, 1) & date <= datetime(2022, 3, 31));
else
    sample_period = find(date >= datetime(2022, 2, 1), 1):Tdata;
end

% Interpolate missing values
%ICU_nation(end) = 12; %Update Saturday values from Tokyo website : https://www.fukushihoken.metro.tokyo.lg.jp/iryo/kansen/corona_portal/info/kunishihyou.html
ICU_nation(109) = 618;
ICU_nation(103) = round((ICU_nation(104) + ICU_nation(102))/2, 0);
ICU_nation(115) = 407;
ICU_nation(116) = 350;
BED(102) = BED(101);

% Retrieve gamma from the time-series of I and R
if DATA_SWITCH == 1
    I = I_data;
    gamma_data = zeros(Tdata, 1);
    %     delta_data = zeros(Tdata,1);
    for i = 1:Tdata
        if I(i) > 0
            gamma_data(i) = dR_data(i) / I_data(i);
            %             delta_data(i) = dD_data(i)/I_data(i);
        end
    end
    gamma = mean(gamma_data(end - RetroPeriod + 1:end));
end

%============ Import vaccine data ============%
Vtable                          = readtable([DATA_PATH 'vaccination_Tokyo_newly.xlsx']);
vac1stweek                      = Vtable.week(1); % index of the first week of the COVID-19 vaccination in Tokyo
V1_elderly                      = zeros(Tdata, 1);
V2_elderly                      = V1_elderly;
V1_others                       = V1_elderly;
V2_others                       = V1_elderly;
V1_elderly(vac1stweek:110)    = Vtable.elderly_first;
V2_elderly(vac1stweek:110)    = Vtable.elderly_second;
V1_others(vac1stweek:110)     = Vtable.others_first;
V2_others(vac1stweek:110)     = Vtable.others_second;

% Medical personels
vaccine_medical                    = readmatrix([DATA_PATH 'vaccine_daily_medical.xls']);
[V1_medical_past, V2_medical_past] = vaccine_daily_to_weekly_table(vaccine_medical, ps, dateEN);

M_first                     = cum_to_new(tablepref.Medical_cum_first); 
M_second                    = cum_to_new(tablepref.Medical_cum_second);
M_ps_first                  = cum_to_new(tablepref.Medical_ps_first); 
M_ps_second                 = cum_to_new(tablepref.Medical_ps_second); 
M_first(find(date == datetime(2021, 8, 5)):end)     = zeros(length(find(date == datetime(2021, 8, 5)):length(M_first)),1);
M_second(find(date == datetime(2021, 8, 5)):end)    = zeros(length(find(date == datetime(2021, 8, 5)):length(M_second)),1);
indM1                       = find(M_ps_first > 0, 1, 'first');   
indM21                      = find(M_ps_second > 0, 1, 'first');

V1_medical                  = (V1_medical_past / ps) * M_ps_first(indM1);
V1_medical(indM1 + 1:end)   = M_first(indM1 + 1:end);
V2_medical                  = (V2_medical_past / ps) * M_ps_second(indM1);
V2_medical(indM21 + 1:end)  = M_second(indM21 + 1:end);

cumsumPastV1                = cumsum(V1_elderly + V1_others + V1_medical);
cumsumPastV2                = cumsum(V2_elderly + V2_others + V2_medical);

V3_d = vaccine_medical(:,5);
Vdata = size(vaccine_medical,1);
V3_w = zeros(length(dateEN),1);
dateV = datetime(vaccine_medical(:,1),'ConvertFrom','excel');
for i = find(dateEN == datetime(2021,12,02)):Tdata
    V3_w(i) = max(V3_d((dateEN(i)-3<dateV) & (dateV<dateEN(i)+3)));
end
V3_ps = tablepref.ps_third(end);
V3_w = round(V3_w * V3_ps);
V3_total = round(cum_to_new(V3_w));
V3_total(end) = round(V3_total(end)*7/4);
%============ Simulated Vaccine Path ============%
VT = zeros(SimPeriod,9);
VT(1,2) = V1_elderly(end-2);
VT(2,2) = V1_elderly(end-1);
VT(3,2) = V1_elderly(end);
VT(1,5) = V1_others(end-2);
VT(2,5) = V1_others(end-1);
VT(3,5) = V1_others(end);
V3_elderly = zeros(Tdata,1);
V3_medical = zeros(Tdata,1);
V3_others = zeros(Tdata,1);

thirdVtable = readtable([DATA_PATH 'tokyo_3rd_vaccination.csv']);
e_share     = thirdVtable.elderly_rate_of_third_dose(end);
V3_elderly  = V3_total * e_share;
V3_medical(V3_total>0)  = (V2_medical(find(V2_medical>0,1):find(V2_medical>0,1) + find(V3_total>0,1,'last')-find(V3_total>0,1))*accept_share);
V3_medical = min(V3_medical, V3_total - V3_elderly);
V3_others   = V3_total - V3_elderly - V3_medical;

indJan2022 = findDateIndex(SimDateEN, "2022,1,1", "2022,1,31");
indFeb2022 = findDateIndex(SimDateEN, "2022,2,1", "2022,2,28");
indMar2022 = findDateIndex(SimDateEN, "2022,3,1", "2022,3,31");
indApr2022 = findDateIndex(SimDateEN, "2022,4,1", "2022,4,30");
indMay2022 = findDateIndex(SimDateEN, "2022,5,1", "2022,5,31");
indJun2022 = findDateIndex(SimDateEN, "2022,6,1", "2022,6,30");

totalVT                 = zeros(SimPeriod,1);
if ~isempty(indJan2022); totalVT(indJan2022) = 200000; end
if ~isempty(indFeb2022); totalVT(indFeb2022) = 500000; end
if ~isempty(indMar2022); totalVT(indMar2022) = 1050000; end
if ~isempty(indApr2022); totalVT(indApr2022) = 500000; end
if ~isempty(indMay2022); totalVT(indMay2022) = 350000; end
if ~isempty(indJun2022); totalVT(indJun2022) = 80000; end


VT(:,9)             = totalVT * 0.5; % Medical
VT9_ind             = find(sum(V2_medical)*accept_share < cumsum(VT(:,9)+sum(V3_medical)),1,'first');
VT(VT9_ind:end,9)   = zeros(SimPeriod-(VT9_ind-1),1);
VT(VT9_ind,9)       = sum(V2_medical)*accept_share - sum(VT(1:VT9_ind-1,9)) - sum(V3_medical);

VT(:,3)             = totalVT - VT(:,9); %Elderly
VT3_ind             = find(sum(V2_elderly)*accept_share<cumsum(VT(:,3))+sum(V3_elderly),1,'first');
VT(VT3_ind:end,3)   = zeros(SimPeriod-(VT3_ind-1),1);
VT(VT3_ind,3)       = sum(V2_elderly)*accept_share -    (sum(VT(1:VT3_ind-1,3)) + sum(V3_elderly));

VT(:,6)             = totalVT -VT(:,9) - VT(:,3); %others
VT6_ind             = find(sum(V2_others)*accept_share_ordinary < cumsum(VT(:,6))+sum(V3_others),1,'first');
VT(VT6_ind:end,6)   = zeros(SimPeriod-(VT6_ind-1),1);
VT(VT6_ind,6)       = sum(V2_others)*accept_share_ordinary - (sum(VT(1:VT6_ind-1,6)) +sum(V3_others));

cumVT3              = (cumsum(VT(:,3))+sum(V3_elderly))/(sum(V2_elderly));
cumVT6              = (cumsum(VT(:,6))+sum(V3_others))/(sum(V2_others));
cumVT9              = (cumsum(VT(:,9))+sum(V3_medical))/(sum(V2_medical));

% figure
% plot(SimDateEN, cumVT3,'r')
% hold on
% plot(SimDateEN, cumVT6,'b')
% plot(SimDateEN, cumVT9, 'k')
% ylim([0 1])

cumsumVT1           = cumsum(VT(:,1) + VT(:,4) + VT(:,7))+ cumsumPastV1(end);
lagged_cumsumVT1    = [cumsumPastV1(end-1);cumsumPastV1(end);cumsumVT1(1:end-2)];
cumsumVT2           = cumsum(VT(:,2) + VT(:,5) + VT(:,8))+ cumsumPastV2(end);
lagged_cumsumVT2    = [cumsumPastV2(end-1);cumsumPastV2(end);cumsumVT2(1:end-2)];

cumsumPastV3        = cumsum(V3_elderly + V3_others + V3_medical);
cumsumVT3           = cumsum(VT(:,3) + VT(:,6) + VT(:,9))+ cumsumPastV3(end);
lagged_cumsumVT3    = [cumsumPastV3(end-1);cumsumPastV3(end);cumsumVT3(1:end-2)];

%============ Constructing the reference level of output ===========%
[potentialGDP, referenceGDP, alpha] = construct_GDP(GDP, TdataGDP);

[Malt, h_all, h_all_se, h, h_se] = estimate_h(M, alpha, TdataGDP, RetroH, hconstant);
% plot GDP and mobilitya
figname = string(['Mobility_GDP_' char(PREF)]);
f       = figure('Name', figname);
plot_mobility(Malt, alpha, Tdata, TdataGDP, YearMonthWeekJP, xtick1, FONT_SIZE, 16)
if FIGURE_SAVE == 1
    saveas(f, [FIGURE_PATH char(PREF) '/MobilityGDPLine.png']);
end

%===== Compute the history of S, I, R, D in the data period ====%
if DATA_SWITCH == 0
    [S, I, R, D,cum_in_R] ...
        = SIRD(Tdata, POP0, N, E1, E2, E3,...
        V1_elderly, V1_medical, V1_others, V2_elderly, V2_medical, V2_others, ...
        V3_elderly,V3_medical,V3_others,...
        gamma, dD, TdataGDP, referenceGDP, alpha);
else
    pastV = zeros(Tdata, 1);
    pastV(3:end) = E1 * (V1_elderly(1:end - 2) + V1_medical(1:end - 2) + V1_others(1:end - 2)) ...
        + (E2 - E1) * (V2_elderly(1:end - 2) + V2_medical(1:end - 2) + V2_others(1:end - 2)) ...
        + (E3 - E2) * (V3_elderly(1:end - 2) + V3_medical(1:end - 2) + V3_others(1:end - 2)) ;
    S = S_data;
    S(2:end) = S(2:end) - cumsum(pastV);
    R = R_data;
    R(2:end) = R(2:end) + cumsum(pastV);
    D = D_data;
end

%--- Compute the history of time-varying parameters ---%
delta           = (D(2:Tdata+1)-D(1:Tdata))./I(1:Tdata);            % death rate
beta_tilde      = (POP0.* N(1:Tdata))./((S(1:Tdata).*I(1:Tdata)));  % overall infection rate, p4 of Fujii and Nakata (2020)
BRN             = beta_tilde ./ (gamma + delta);                    % Basic reproduction number
ERN             = (S(1:end-1)/POP0).* BRN; % effective reproduction number
if hconstant == 0
    beta        = beta_tilde./(1+h_all*alpha).^k; % raw infection rate
elseif hconstant == 1
    beta        = beta_tilde./(1+(h_all(2)/h_all(1))*alpha).^k;
end
%--- Time series average ---%
equal_weight_retro  = ones(RetroPeriod,1)/RetroPeriod;
equal_weight_sample = ones(length(sample_period), 1)/length(sample_period);
I_weight     = I(sample_period)/sum(I(sample_period));

ICU_nation_rate     = (ICU_nation(2:Tdata+1) - ICU_nation(1:Tdata) ...
                    + gamma_ICU_nation.*ICU_nation(1:Tdata) ...
                    + dD(1:Tdata))./I(1:Tdata);
ICU_pref_rate       = (ICU_pref(2:Tdata+1) - ICU_pref(1:Tdata) ...
                    + gamma_ICU_pref.*ICU_pref(1:Tdata) ...
                    + dD(1:Tdata))./I(1:Tdata);
Hospital_rate       = (hospital(2:Tdata+1) - hospital(1:Tdata) ...
                    + gamma_Hospital.*hospital(1:Tdata))./I(1:Tdata);
newICU_pref_rate    = (newICU_pref(2:Tdata+1) - newICU_pref(1:Tdata) ...
                    + gamma_newICU_pref.*newICU_pref(1:Tdata) ...
                    + dD(1:Tdata))./I(1:Tdata);
newICU_pref_rate(1:find(newICU_pref > 0, 1, 'first')-1) = nan;

[delta_average, delta_se] = ts_avg(delta, sample_period, I_weight);
[ICU_nation_rate_average] = ts_avg(ICU_nation_rate, sample_period, I_weight);
[ICU_pref_rate_average]   = ts_avg(ICU_pref_rate, sample_period, I_weight);
[Hospital_rate_average]   = ts_avg(Hospital_rate, sample_period, I_weight);
[newICU_pref_rate_average] = ts_avg(newICU_pref_rate, sample_period, I_weight);
[simple_beta_avg, beta_se] = ts_avg(beta, Tdata-RetroPeriod+1:Tdata, equal_weight_retro);

ICU_nation_rate_average     = ICU_nation_rate_average  *ICU_nation_adjustment;
ICU_pref_rate_average       = ICU_pref_rate_average    *ICU_pref_adjustment;
Hospital_rate_average       = Hospital_rate_average    *Hospital_adjustment;
newICU_pref_rate_average    = newICU_pref_rate_average *newICU_pref_adjustment;

past_omicron_share  = zeros(Tdata,1);
past_omicron_share(find(dateEN == datetime(2021, 12, 23)))= 0.079;
past_omicron_share(find(dateEN == datetime(2021, 12, 30)))= 0.446;
past_omicron_share(find(dateEN == datetime(2022, 1, 6))) = 0.81;
past_omicron_share(find(dateEN == datetime(2022, 1, 13)))= 0.894; %https://www3.nhk.or.jp/news/html/20220118/k10013437211000.html
past_omicron_share(find(dateEN == datetime(2022, 1, 20)))= 0.95;
past_omicron_share(find(dateEN == datetime(2022, 1, 27)):Tdata) = 1.00;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Simulaiton starts here %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=========== simulation setting =============%
% Sim_BA2_share = variant_share_two_target(indApr2022(1), indApr2022(end), 0.05, 0.8, 1, SimPeriod);
% BA2_share       = [zeros(Tdata, 1); Sim_BA2_share];
T1_ind = find(date == datetime(2022, 2, 17));
T2_ind = find(date == datetime(2022, 3, 17));
share_T1 = 0.074;
share_T2 = 0.523;
BA2_share = variant_share_two_target(T1_ind, T2_ind, share_T1, share_T2, 1, Tdata+SimPeriod);
Sim_BA2_share = BA2_share(Tdata+1:end);
BA2_relative_infectivity = 1.3;
BA2_infectivity_path = (1-Sim_BA2_share)  * 1 + Sim_BA2_share * BA2_relative_infectivity;
past_BA2_infectivity_path = (1 - BA2_share(1:Tdata)) + BA2_share(1:Tdata) * BA2_relative_infectivity;


SoE_off_date        = find(date == datetime(2022, 3, 24)) + 1; % まん防 終了...
beta_jump_duration  = find(SimDateEN == datetime(2022, 5, 19));
% 現実世界とのLagを考慮して+1
SoE_date            = SimPeriod+1;  %End of Feb., First wekk of Mar., Second week of Mar
alphaAug2021        = mean( alpha(findDateIndex(date, "2021,8,1", "2021,8,31"))  );

state               = 1;    %=1 if currently under SOE
th_on_percent       = 0.35;
th_on               = max(newICU_limit_pref_vec)*th_on_percent; %8000*7; %ICU_limit_pref_vec(end)*0.5; %
th_off              = th_on*0.5; %1000 * 7;  %ICU_limit_pref_vec(end)*0.25; %
th_on_N             = 25000*7;
th_off_N            = th_on_N * 0.1;
relative_beta_SOE   = 0.3;
alpha_on            = alphaAug2021;
alpha_off           = 0.00;
alpha_jump          = 0;

originalE1          = E1;
originalE2          = E2;
originalE3          = E3;

InitialValues       = [S(end), I(end), R(end), D(end), ...
                       ICU_nation(end), ICU_pref(end), hospital(end), newICU_pref(end)];

gammaT                  = gamma             * ones(SimPeriod,1);
gamma_ICU_nation_path   = gamma_ICU_nation  * ones(SimPeriod,1);
gamma_ICU_pref_path     = gamma_ICU_pref    * ones(SimPeriod,1);
gamma_newICU_pref_path  = gamma_newICU_pref * ones(SimPeriod,1);
gamma_Hospital_path     = gamma_Hospital    * ones(SimPeriod,1);

seasonality = seasonal_adjustment(retro_lb, retro_ub, dateD, SimPeriod+DR+1, seasonal_effect);

% omicron_realtive_severity           = 0.2;%[0.4,0.2,0.05];%[1, 0.5, 0.25];
% omicron_realtive_severity_nation    = omicron_realtive_severity; %[0.6,0.35,0.10];%[1, 0.5, 0.25];
% omicron_realtive_hospitalized_rate  = omicron_realtive_severity; 

% Simulation parameters for loop
%iX
ori_BRN_goal_vec = [4.5, 3.5, 2.5];
ori_beta_goal_vec = ori_BRN_goal_vec .* (gamma + delta_average);
omicron_relative_infectivity = 1.0;%1.2 

betaT_temp_ini_vec = [0, 0, 0];   
betaT_rho_vec      = [0.99, 0.99, 0.99];


%iY
target_relative_severity = [1.25 ,1, 0.75];
transitionLength         = 4;
beta_jump_vec        = [0, 0, 0, 0]; 

severity_new_pref_standard_vector = target_relative_severity;
severity_pref_standard_vector     = target_relative_severity;
hospitalization_standard_vector   = target_relative_severity;
severity_nation_standard_vector   = target_relative_severity;
fatality_rate_standard_vector     = target_relative_severity;

%iX
omicron_E2_vector = 1;%[1, 0.6, 0.2]; %Relative

xvec = ori_BRN_goal_vec;
yvec = target_relative_severity;
zvec = omicron_E2_vector;
nX = length(xvec);
nY = length(yvec);
nZ = length(zvec);

% Scenario Name and Line
figfolder = 'Baseline';



titlevec            = {'Scenario', 'シナリオ'};
figname_xvar        = '_relative_severity_';

Scenario = cell(nX,1);
ScenarioEN = cell(nX,1);
for iX = 1:nX
    Scenario{iX}            = ['基本再生産数 = ',num2str(xvec(iX))];
    ScenarioEN{iX}          = ['Basic Reproduction Number = ',num2str(xvec(iX))];
end

Scenario_vec        = [ScenarioEN,Scenario];
fig_Scenario_vec    = ["Scenario_A", "Scenario_B", "Scenario_C"];
linecolor           = {"r", "k", "b"};
LineStyles          = {"-", "-", "-"};
lineWidth           = [1.5, 1.5, 1.5];
markertype          = {'o','o','o'};

lineNameJP = cell(nY,1);
lineNameEN = cell(nY,1);
for iY = 1:nY
    lineNameJP{iY}          = ['2月~3月の重症化率の',num2str(yvec(iY)*100),'%'];
    lineNameEN{iY}          = ['Severity Rate ',num2str(yvec(iY)*100),'% (Relatieve to Feburary and March) '];
end

xmax                = indJun2022(end) + Tdata;
xmin                = find(date == datetime(2021, 7, 1));
column_num_main     = nZ;

% initialize_matirix
DMat                = nan(nX, nY, nZ);
AlphaMat            = DMat;

SimData                 = nan(SimPeriod + 1, length(InitialValues), nX, nY, nZ);
AlphaPath               = nan(SimPeriod, nX, nY, nZ);
NPath                   = AlphaPath;
Sim_dD                  = AlphaPath;
SimERN                  = AlphaPath;
SimBRN                  = AlphaPath;
betaPath                = AlphaPath;
betaTildePath           = AlphaPath;
deltaPath               = AlphaPath;
delta_ICU_nationPath    = AlphaPath;
delta_ICU_prefPath      = AlphaPath;
delta_newICU_prefPath   = AlphaPath;
delta_HospitalPath      = AlphaPath;

SimICU_nation           = nan(SimPeriod+1, nX, nY, nZ);
SimICU_pref             = SimICU_nation;
SimNewICU_pref          = SimICU_nation;
SimHospital             = SimICU_nation;


%==== Simulation Starts ======%

for iX = 1:nX %different figures
    ori_beta_goal = ori_beta_goal_vec(iX);
    
    betaT_temp_ini = betaT_temp_ini_vec(iX);
    beta_rho = betaT_rho_vec(iX);

    
    for iY = 1:nY %different curves within a figure
        severity_nation_standard    = severity_nation_standard_vector(iY);
        severity_new_pref_standard  = severity_new_pref_standard_vector(iY);
        severity_pref_standard      = severity_pref_standard_vector(iY);
        hospitalization_standard    = hospitalization_standard_vector(iY);
        fatality_rate_standard      = fatality_rate_standard_vector(iY);
        beta_jump_path = zeros(SimPeriod,1);
%         beta_jump_path(SoE_off_date:SoE_off_date+beta_jump_duration)   = linspace(beta_jump_vec(iY),1,beta_jump_duration+1);
%         beta_jump_path(SoE_off_date+beta_jump_duration+1:end)              = ones(size(beta_jump_path(SoE_off_date+beta_jump_duration+1:end)));
        
        
        for iZ = 1:nZ %different curves within a figure

            omicronE3 = originalE3;
            omicronE2 = originalE2 * omicron_E2_vector(iZ);
            omicronE1 = originalE1 * omicron_E2_vector(iZ);
            
            % calculate vaccine effectiveness
            E1 = omicronE1; %vector
            E2 = omicronE2; %vector
            E3 = omicronE3; %vector
            
            D1 = (D1 - E1)/(1-E1);      % Find reduction of death conditoinal on infection after first does
            D2 = (D2 - E2)/(1-E2);      % Find reduction of death conditoinal on infection after second does
            D3 = (D3 - E3)/(1-E3);      % Find reduction of death conditoinal on infection after second does
            
            %calculate_omicron_share
            sim_omicron_share   = ones(SimPeriod,1);
            omicron_share       = [past_omicron_share; sim_omicron_share];
            
            VE = E1.*(VT(:,1)+VT(:,4)+VT(:,7))...
                +(E2-E1).*(VT(:,2)+VT(:,5)+VT(:,8)) ...
                +(E3-E2).*(VT(:,3)+VT(:,6)+VT(:,9));
            VE_prev = originalE1.*(V1_elderly+V1_medical+V1_others)...
                +(originalE2-originalE1).*(V2_elderly+V2_medical+V2_others)...
                +(originalE3-originalE2).*(V3_elderly+V3_medical+V3_others);
            V = [VE_prev(end-1);VE_prev(end);VE(1:end-2)];
            
            %Construct betapath
            %Add Influence of BA.1, BA.2, and Seasonality 
            beta_goal   = ori_beta_goal;     

         
            betaT      = beta_goal* omicron_relative_infectivity...
                .* BA2_infectivity_path...
                .* transpose(seasonality(1:SimPeriod)) ...
                .* construct_onetime_AR1(beta_rho, betaT_temp_ini, SimPeriod);
            
            initial = beta(end)/betaT(1);
            transitionPath = linspace(initial,1,transitionLength+1);
            betaT(1:transitionLength)  = betaT(1:transitionLength)  .* transitionPath(2:end)';
%             betaT      = adjustment_data_to_simulation(beta, betaT, transitionLength);
                                       % ...
%                     .* beta_jump_path(SoE_off_date+1:end) ...
%                     + (1 - beta_jump_path(SoE_off_date+1:end)) .* betaT(SoE_off_date);
            
           
            betaBox     = [beta; betaT]; 
            
            %Construct alphapath
            alpha_recovery_initial = alpha(Tdata) + (alpha_off - alpha(Tdata)) * alpha_jump;
            durationAfterSOE = DR - (Tdata + 1 - SoE_off_date)*(Tdata>=SoE_off_date);
            alphaAfterSOE   = [alpha_recovery_initial:(alpha_off - alpha_recovery_initial) / (durationAfterSOE):alpha_off];
            %alphaT          = [alphaAfterSOE'; alpha_off * ones(SimPeriod, 1)];
            alpha_date_index = max(1,SoE_off_date-Tdata-1);
            currentDateIndex = Tdata - (SoE_off_date - 1) + 1;
            alphaT          = [alphaAfterSOE(currentDateIndex:end)';alpha_off * ones(SimPeriod-alpha_date_index+2,1)];
            alphaBox        = [alpha; alphaT];
            
            % death rate, severity rate, hospoital rate paths
            deltaT                = delta_average               *ones(SimPeriod,1);
            SimICU_nation_rate    = ICU_nation_rate_average     *ones(SimPeriod,1);
            SimICU_pref_rate      = ICU_pref_rate_average       *ones(SimPeriod,1);
            SimNewICU_pref_rate   = newICU_pref_rate_average    *ones(SimPeriod,1);
            SimHospital_rate      = Hospital_rate_average       *ones(SimPeriod,1);
            
            % Adjust for different sevirity rates based on scneario settings 
            deltaT               = deltaT * fatality_rate_standard;
            SimICU_nation_rate   = SimICU_nation_rate   *  severity_nation_standard;
            SimNewICU_pref_rate  = SimNewICU_pref_rate  *  severity_new_pref_standard;
            SimHospital_rate     = SimHospital_rate     *  hospitalization_standard;
            SimICU_pref_rate     = SimICU_pref_rate     *  severity_pref_standard;
            
            deltaT   = ...
                adjustment_data_to_simulation(delta, deltaT, transitionLength);
            SimICU_nation_rate   = ...
                adjustment_data_to_simulation(ICU_nation_rate, SimICU_nation_rate, transitionLength);
            SimICU_pref_rate   = ...
                adjustment_data_to_simulation(ICU_pref_rate, SimICU_pref_rate, transitionLength);
            SimNewICU_pref_rate   = ...
                adjustment_data_to_simulation(newICU_pref_rate, SimNewICU_pref_rate, transitionLength);
            SimHospital_rate   = ...
                adjustment_data_to_simulation(Hospital_rate, SimHospital_rate, transitionLength);
            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation SIRD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [DMat(iX,iY,iZ), AlphaMat(iX,iY,iZ), AlphaPath(:, iX,iY,iZ), ...
                SimData(:, :, iX,iY,iZ), NPath(:, iX,iY,iZ), SimERN(:, iX,iY,iZ),...
                SimICU_nation(:, iX,iY,iZ), SimICU_pref(:, iX,iY,iZ), SimHospital(:, iX,iY,iZ), SimNewICU_pref(:, iX,iY,iZ),...
                betaPath(:,iX,iY,iZ),betaTildePath(:,iX,iY,iZ),SimBRN(:,iX,iY,iZ)] = ...
                Covid_projection(InitialValues,betaT,gammaT,deltaT,...
                     SimICU_nation_rate, SimICU_pref_rate, ...
                     SimHospital_rate,SimNewICU_pref_rate, ...
                     h,k,POP0, V, E1, E2, E3, ...
                     hconstant,DR, gamma_ICU_nation_path, gamma_ICU_pref_path, ...
                     gamma_Hospital_path,gamma_newICU_pref_path,  alphaBox(Tdata + 1:end));
            Sim_dD(:, iX, iY, iZ) = squeeze(SimData(2:end,4,iX,iY,iZ)-SimData(1:end-1,4,iX,iY,iZ));
            
            deltaPath(:, iX,iY,iZ)  = deltaT;
            SimICU_nation_rate_Path(:, iX,iY,iZ)    = SimICU_nation_rate;
            SimICU_pref_rate_Path(:, iX,iY,iZ)      = SimICU_pref_rate;
            SimNewICU_pref_rate_Path(:, iX,iY,iZ)   = SimNewICU_pref_rate;
            SimHospital_rate_Path(:, iX,iY,iZ)      = SimHospital_rate;
            
        end
        
    end
    
end

% Save variable names and values for backdata files
minAlpha = alpha_off; % 経済損失0の基準
AlphaM = AlphaMat;
% AlphaM = AlphaMat(~isnan(AlphaMat));
AlphaM = (AlphaM - minAlpha) * PREF_GDP * 10000 * SimPeriod / 52;
DM = DMat;
% DM = DMat(~isnan(DMat));
% BackDataDA(1:nX*nY*nZ) = [round(AlphaM') / 10000, round(DM'), Scenario'];
%--- Record how many times on and off are triggered ---%
waves = zeros(nX,nY,nZ);
for ixx = 1:nX
    for iyy = 1:nY
        for izz = 1:nZ
            svec = zeros(SimPeriod - 1, 1);
            
            for t = 1:SimPeriod - 1
                svec(t) = AlphaPath(t + 1, ixx,iyy,izz) - AlphaPath(t, ixx,iyy,izz);
            end
            
            waves(ixx,iyy,izz) = sum(svec > 0);
        end
    end
end

% disp('Pessimistic')
% disp(round(squeeze(NPath(3,1,1,1)/7)'))
% disp('Baseline')
% disp(round(squeeze(NPath(3,2,1,1)/7)'))
% disp('Optimistic')
% disp(round(squeeze(NPath(3,3,1,1)/7)'))

toc;


% %==================== Plot and Backdata ====================%

lgdLocation = 'northwest';
column_num = 3;

for l = 1:2 %1:2 when english version needed
    
    for iX = 1:nX
        if l == 1
            lineName = lineNameEN;
        elseif l == 2
            lineName = lineNameJP;
        end
        
        % Generate graphs for the website
        lng = LANGUAGE{l};
        figname = [char(figname_main) '_' char(fig_Scenario_vec(iX)) '_' char(lng)];
        f = figure('Name', figname);
        t = tiledlayout(3,3, 'TileSpacing', 'compact');
        title(t,Scenario_vec(iX,l),'FontSize',20)
        f.Position = [100, 100, 1200, 800];
        
        
        %--- Number of people who are in ICU (local standard)---%
        %         subplot(3, 3, 1)
        nexttile
        title_vec = ["ICU (Old Tokyo Standard)", "重症患者数（旧都基準）"];
        plot_3Dfunction(ICU_pref(2:end), SimICU_pref(2:end,:,:), iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            FONT_NAME, FONT_SIZE, LEGEND_FONT_SIZE, AXIS_FONT_SIZE,Y_FORMAT,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        hold on
        plot([ICU_limit_pref_vec; ones(SimPeriod, 1) * max(ICU_limit_pref_vec)], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        text(xmin, ICU_limit_pref_vec(xmin) * 0.85, '100%', 'FontSize', FONT_SIZE)
        hold on
        plot([ICU_limit_pref_vec * 0.5; ones(SimPeriod, 1) * max(ICU_limit_pref_vec)*0.5], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        hold on
        text(xmin, ICU_limit_pref_vec(xmin) * 0.5 * 0.8, '50%', 'FontSize', FONT_SIZE)
        ylim([0 max(ICU_limit_pref_vec)*1.2])
       
        %         subplot(3, 3, 2)
        nexttile
        
        title_vec = ["ICU (New Tokyo Standard)", "重症患者数（新都基準）"];

        plot_3Dfunction(newICU_pref(2:end), SimNewICU_pref(2:end,:,:), iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            FONT_NAME, FONT_SIZE, LEGEND_FONT_SIZE, AXIS_FONT_SIZE,Y_FORMAT,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        hold on
        plot([newICU_limit_pref_vec; ones(SimPeriod, 1) * max(newICU_limit_pref_vec)], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        text(xmin, newICU_limit_pref_vec(xmin) * 0.85, '100%', 'FontSize', FONT_SIZE)
        hold on
        plot([newICU_limit_pref_vec * th_on_percent; ones(SimPeriod, 1) * th_on], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        hold on
        %         plot([nan(Tdata,1); ones(SimPeriod, 1) * th_on], '-r', 'HandleVisibility', 'off', 'LineWidth', 1.0)
        %         hold on
        %         plot([nan(Tdata,1); ones(SimPeriod, 1) * th_off], '-b', 'HandleVisibility', 'off', 'LineWidth', 1.0)
        text(xmin, newICU_limit_pref_vec(xmin) * th_on_percent * 0.8, [num2str(th_on_percent*100) '%'], 'FontSize', FONT_SIZE)
        ylim([0 max(newICU_limit_pref_vec)*1.2])
        
        %--- Number of people who are in ICU (Naitonal standard) ---%
        %         subplot(3, 3, 3)
%         nexttile
%         title_vec = ["ICU (National Standard)", "重症患者数（国基準）"];
%         plot_3Dfunction(ICU_nation(2:end), SimICU_nation(2:end,:,:), iX, ...
%             WeekNumber, YearMonth, xmin, xmax, ...
%             FONT_NAME, FONT_SIZE, LEGEND_FONT_SIZE, AXIS_FONT_SIZE,Y_FORMAT,...
%             lgdLocation, column_num, l, title_vec, ...
%             lineWidth,linecolor, LineStyles, lineName)
%         Handle = legend;
%         set(Handle, 'Visible', 'off');
%         hold on
%         plot([BED; ones(SimPeriod, 1) * BED(end)], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
%         text(xmin, BED(xmin) * 0.85, '100%', 'FontSize', FONT_SIZE)
%         hold on
%         plot([BED * 0.5; ones(SimPeriod, 1) * BED(end) * 0.5], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
%         text(xmin, BED(xmin) * 0.4, '50%', 'FontSize', FONT_SIZE)
%         ylim([0 max(BED)*1.5])
%         
        %--- Number of newly hospitalized ---%
        %         subplot(3, 3, 3)
        nexttile
        title_vec = ["Hospitalized Patients", "入院患者数"];
        plot_3Dfunction(hospital(2:end), SimHospital(2:end,:,:), iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            FONT_NAME, FONT_SIZE, LEGEND_FONT_SIZE, AXIS_FONT_SIZE,Y_FORMAT,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        hold on
        plot([Hospital_limit_vec; ones(SimPeriod, 1) * max(Hospital_limit_vec)], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        text(xmin, Hospital_limit_vec(xmin) * 0.85, '100%', 'FontSize', FONT_SIZE)
        hold on
        plot([Hospital_limit_vec * 0.5; ones(SimPeriod, 1) * max(Hospital_limit_vec) * 0.5], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        text(xmin, Hospital_limit_vec(xmin) * 0.4, '50%', 'FontSize', FONT_SIZE)
        ylim([0 max(Hospital_limit_vec)*1.2])
%         
        %         subplot(3, 3, 4)
        nexttile
        title_vec = ["New Deaths (Daily Average)", "新規死亡者数（1日平均）"];
        plot_3Dfunction(dD / 7, Sim_dD/ 7, iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            FONT_NAME, FONT_SIZE, LEGEND_FONT_SIZE, AXIS_FONT_SIZE,Y_FORMAT,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        %         ylim([0 ymax_D])
        %
        %--- Number of new cases ---%
        %         subplot(3, 3, 5)
        nexttile
        title_vec = ["New Cases (Daily Average)", "新規感染者数（1日平均）"];
        plot_3Dfunction(N/7, NPath/7, iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            FONT_NAME, FONT_SIZE, LEGEND_FONT_SIZE, AXIS_FONT_SIZE,Y_FORMAT,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        %         ylim([0,10000])
        %--- GDP Path ---%
        %         Y_FORMAT = '%.2f';
        %         subplot(3, 3, 6)
        nexttile
        title_vec = ["GDP", "GDP"];
        plot_3Dfunction(100*(1-alpha), 100*(1-AlphaPath), iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            FONT_NAME, FONT_SIZE, LEGEND_FONT_SIZE, AXIS_FONT_SIZE,Y_FORMAT,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        ylim([90 100])
        Handle = legend;
        set(Handle, 'Visible', 'off');
        
        %         subplot(3, 3, 7)
        nexttile
        title_vec = ["Transitions of ERN", "実効再生産数"];
        plot_3Dfunction(ERN, SimERN, iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            FONT_NAME, FONT_SIZE, LEGEND_FONT_SIZE, AXIS_FONT_SIZE,Y_FORMAT,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        
        
        %         subplot(3, 3, 8) %Cumulative vaccinated
        nexttile
        
        VT_23 = VT;
        VT_23(:,[1,4,7]) = [];
        plot_vaccinepath_separate_percentage(2, VT_23, V2_medical, V3_medical, V2_elderly, V3_elderly, V2_others, V3_others, ps, YearMonthWeek(:, l), WeekNumber, Tdata, FONT_SIZE, 8, FONT_NAME, xmin, xmax, l, POP0);
        yticks(0:10:100)
        ylim([0 100])
        grid on
        xlim([xmin, xmax])
        %     xticks(find(WeekNumber == 1 & abs(1 - mod(MonthNumber, 4)) < 0.01))
        xticks(find(WeekNumber == 1))
        %     xticklabels(YearMonthWeekJP(xticks))
        xticklabels(YearMonthJP(xticks))
        lgd = legend;
        lgd.FontSize = 12;
        lgd.Location = 'Northwest';
        ax = gca;
        ax.YAxis.FontSize = AXIS_FONT_SIZE;
        ax.XAxis.FontSize = AXIS_FONT_SIZE;
        %     xline(xline_ind, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        %--- Trade-off Curve ---%
        %         subplot(3, 3, 9)
        
        nexttile
        title_vec = ["Transition of the Share of Omicron Variant(BA.2 sub-lineage)", "オミクロン株BA.2系統割合の推移"];
        
        plot(BA2_share)
        hold on
        xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
        xlim([xmin xmax])
        ylim([0 1])
        xticks(find(WeekNumber == 1))
        xticklabels(YearMonth(xticks,l))
        title(title_vec(l),'FontSize',FONT_SIZE,'FontWeight','normal','FontName',FONT_NAME)
        ax = gca;
        ax.YAxis.FontSize = AXIS_FONT_SIZE;   ax.XAxis.FontSize = AXIS_FONT_SIZE;   ax.YAxis.Exponent = 0;
        xtickangle(45)
        
        Handle = legend;
        set(Handle, 'Visible', 'off');
        ylim([0 1])
        grid on
        ax = gca;
        box on
        xtickangle(45)
        
        
        %         for iyy = 1:nY
        %             for izz = 1:nZ
        %                 scatter(AlphaM(iX,iyy,izz),DM(iX,iyy,izz),250,...
        %                     markertype{iyy,izz}, linecolor{iyy,izz},'filled');
        %                 hold on
        %             end
        %         end
        %         if l == 1
        %             xlabel('Output Loss (hundred million yen)','FontSize',FONT_SIZE)
        %             ylabel('Cumulative Deaths','FontSize',FONT_SIZE)
        %             title('Relationship between Covid-19 and Output(within 10 years)','FontSize',FONT_SIZE,'FontWeight','normal')
        %         elseif l == 2
        %             xlabel('経済損失 (兆円)','FontSize',FONT_SIZE,'FontName',FONT_NAME)
        %             ylabel('累計死亡者数','FontSize',FONT_SIZE,'FontName',FONT_NAME)
        %             title('コロナ感染と経済の関係(今後10年)','FontSize',FONT_SIZE,'FontWeight','normal','FontName',FONT_NAME)
        %         end
        %         xlim([0,inf])
        %         xtickangle(45)
        %         grid on
        %         ax = gca;
        %         ax.YAxis.FontSize = 12;
        %         ax.XAxis.FontSize = 12;
        %         ax.YAxis.Exponent = 0;
        %         ax.XAxis.Exponent = 0;
        %         ytickformat('%,6.0f')
        %         box on
        %         grid on
        %
        lgd = legend(nexttile(1));
        lgd.Layout.Tile  = 'south';
        lgd.NumColumns = 3;
        lgd.FontSize = LEGEND_FONT_SIZE_MAIN;
        %     xticks(0:2.5:30)
        %     ylim([max(0, mean(DM(iX,:,:),[2,3]) - 7500) mean(DM(iX,:,:),[2,3]) + 7500])
        
        if FIGURE_SAVE == 1
            %         saveas(f, [home 'Figures/' char(PREF) '/' char(figname_main) char(lng) '.png']);
            saveas(f, [FIGURE_PATH char(PREF) '/' char(figfolder) '/' char(figname) '.png']);
        end
        
        
    end
end %End of LANGUAGE loop = figure loop

%% Other Plot
LEGEND_FONT_SIZE = 12;
AXIS_FONT_SIZE = 12;
lineWidth = 2.0*ones(nY);
plot_parameter

% figure('Name','SR_transition')
% title_vec = ["Transition of S and R", "SとRの推移"];
% for iX = 1:nX
%     for iY = 1:nY
%         plot([S(2:end); squeeze(SimData(2:end,1,iX,iY))])
%         hold on
%         plot([R(2:end); squeeze(SimData(2:end,3,iX,iY))])
%         hold on
%     end
% end

%% Backdata
BRNpast = beta_tilde ./ (gamma + delta);
BackDataN           = zeros(8 + SimPeriod, nX, nY, nZ);
BackDataAlpha       = BackDataN;
BackDataERN         = BackDataN;
BackDataBRN         = BackDataN;
BackDatadD          = BackDataN;
BackDataICU_pref    = BackDataN;
BackDataNewICU_pref = BackDataN;
BackDataICU_nation  = BackDataN;
BackDataHospital    = BackDataN;
for iX = 1:nX
    for iY = 1:nY
        for iZ = 1:nZ
            BackDataN(:, iX, iY, iZ)            = [N(end - 7:end);          NPath(:, iX, iY, iZ)];
            BackDataAlpha(:, iX, iY, iZ)        = [alpha(end - 7:end);      AlphaPath(:, iX, iY, iZ)];
            BackDataERN(:, iX, iY, iZ)          = [ERN(end - 7:end);        SimERN(:, iX, iY, iZ)];
            BackDataBRN(:, iX, iY, iZ)          = [BRNpast(end - 7:end);    SimBRN(:, iX, iY, iZ)];
            BackDatadD(:, iX, iY, iZ)           = [dD(end - 7:end);         Sim_dD(:, iX, iY, iZ)];
            BackDataICU_pref(:, iX, iY, iZ)     = [ICU_pref(end - 7:end);   SimICU_pref(2:end, iX, iY, iZ)];
            BackDataNewICU_pref(:, iX, iY, iZ)     = [newICU_pref(end - 7:end);   SimNewICU_pref(2:end, iX, iY, iZ)];
            BackDataICU_nation(:, iX, iY, iZ)   = [ICU_nation(end - 7:end); SimICU_nation(2:end, iX, iY, iZ)];
            BackDataHospital(:, iX, iY, iZ)     = [hospital(end - 7:end);   SimHospital(2:end, iX, iY, iZ)];
        end
    end
end
%%
if DATA_SAVE == 1
%     for iX = 1:nX
%         titleN = strings(1, 1 + nZ * 9);
%         titleN(1) = "週";
%         for ti = 1:nZ
%             titleN(1, 1 + ti) = append("新規感染者数（", Scenario(iX), "）");
%             titleN(1, 1 + nZ + ti) = append("経済活動（", Scenario(iX), "）");
%             titleN(1, 1 + nZ * 2 + ti) = append("実効再生産数（", Scenario(iX), "）");
%             titleN(1, 1 + nZ * 3 + ti) = append("基本再生産数（", Scenario(iX), "）");
%             titleN(1, 1 + nZ * 4 + ti) = append("入院患者数（", Scenario(iX), "）");
%             titleN(1, 1 + nZ * 5 + ti) = append("重症者数_国基準（", Scenario(iX), "）");
%             titleN(1, 1 + nZ * 6 + ti) = append("重症者数_旧都基準（", Scenario(iX), "）");
%             titleN(1, 1 + nZ * 7 + ti) = append("重症者数_新都基準（", Scenario(iX), "）");
%             titleN(1, 1 + nZ * 7 + ti) = append("重症者数_新都基準（", Scenario(iX), "）");
%             titleN(1, 1 + nZ * 7 + ti) = append("重症者数_新都基準（", Scenario(iX), "）");
%             titleN(1, 1 + nZ * 8 + ti) = append("新規死亡者数（", Scenario(iX), "）");
%         end
%         TN = table([
%             titleN;
%             YearMonthWeekJP(Tdata - 7:end - 1), ...
%             squeeze(round(BackDataN(:, iX, iY, :) / 7)), ...
%             squeeze(round(100 * (1 - BackDataAlpha(:, iX, iY, :)), 1)), ...
%             squeeze(round(BackDataERN(:, iX, iY, :), 2)), ...
%             squeeze(round(BackDataBRN(:, iX, iY, :), 2)), ...
%             squeeze(round(BackDataHospital(:, iX, iY, :))), ...
%             squeeze(round(BackDataICU_nation(:, iX, iY, :))), ...
%             squeeze(round(BackDataICU_pref(:, iX, iY, :))), ...
%             squeeze(round(BackDataNewICU_pref(:, iX, iY, :))), ...
%             squeeze(round(BackDatadD(:, iX, iY, :) / 7))
%             ]);
%         writetable(TN, [FIGURE_PATH char(PREF) '\' char(figfolder) '\BackData_' char(figname_main) '_' char(ScenarioEN(iX)) '.xls'], 'Sheet', '新規感染者数（1日平均）', 'WriteVariableNames', false);
%     end
    for iX = 1:nX
        titleN = strings(1, 15);
        titleN(1) = "週";
        titleN(1, 2) = "新規感染者数";
        titleN(1, 3) = "経済活動";
        titleN(1, 4) = "実効再生産数";
        titleN(1, 5) = "基本再生産数";
        titleN(1, 6) = append("入院患者数（", lineNameJP{1}, "）");
        titleN(1, 7) = append("入院患者数（", lineNameJP{2}, "）");
        titleN(1, 8) = append("入院患者数（", lineNameJP{3}, "）");
        titleN(1, 9) = append("重症者数_旧都基準（", lineNameJP{1}, "）");
        titleN(1, 10) = append("重症者数_旧都基準（", lineNameJP{2}, "）");
        titleN(1, 11) = append("重症者数_旧都基準（", lineNameJP{3}, "）");
        titleN(1, 12) = append("重症者数_新都基準（", lineNameJP{1}, "）");
        titleN(1, 13) = append("重症者数_新都基準（", lineNameJP{2}, "）");
        titleN(1, 14) = append("重症者数_新都基準（", lineNameJP{3}, "）");
        titleN(1, 15) = "新規死亡者数";

        
        TN = table([
            titleN;
            YearMonthWeekJP(Tdata - 7:end - 1), ...
            squeeze(round(BackDataN(:, iX, 2, 1) / 7)), ...
            squeeze(round(100 * (1 - BackDataAlpha(:, iX, 2, 1)), 1)), ...
            squeeze(round(BackDataERN(:, iX, 2, 1), 2)), ...
            squeeze(round(BackDataBRN(:, iX, 2, 1), 2)), ...
            squeeze(round(BackDataHospital(:, iX, :, 1))), ...
            squeeze(round(BackDataICU_pref(:, iX, :, 1))), ...
            squeeze(round(BackDataNewICU_pref(:, iX, :, 1))), ...
            squeeze(round(BackDatadD(:, iX, 2, 1) / 7))
            ]);
        writetable(TN, [FIGURE_PATH, char(PREF) '/' char(figfolder) '/BackData_' char(figname_main) '_' char(ScenarioEN(iX)) '.xls'], 'Sheet', '新規感染者数（1日平均）', 'WriteVariableNames', false);
    end
end

% %%
% if DATA_SAVE == 1
%     titleN = strings(1, 1 + length(TH_index) * 8);
%     titleN(1) = "週";
%
%     TN = table([titleN; YearMonthWeekJP(Tdata - 7:end - 1), ...
%                 round(BackDataN(:, 1:length(TH_index)) / 7), ...
%                 round(100 * (1 - BackDataAlpha(:, 1:length(TH_index))), 1), ...
%                 round(BackDataERN(:, 1:length(TH_index)), 2), ...
%                 round(BackDataICU_nation(:, 1:length(TH_index))), ...
%                 round(BackDataICU_pref(:, 1:length(TH_index))), ...
%                 round(BackDatadD(:, 1:length(TH_index)) / 7), ...
%                 round(BackDataHospital(:, 1:length(TH_index))),...
%                 round(BackDataBRN(:, 1:length(TH_index)), 2)]);
%
%     titleAD = ["経済損失（兆円）", "死亡者数", "ケース"];
%     TAD = table([titleAD; BackDataDA(1:length(TH), :)]);
%
%     %     TVAC = table([titleVAC; ...
%     %         MonthWeekJP(x_left-3:x_right),...
%     %         round(BackData_Area_nv_all(x_left-3:x_right,:))*10000]);
%     TVAC_new = table([titleVAC; ...
%                 YearMonthWeekJP(x_left - 3:x_right), ...
%                     round(BackData_Area_nv_all(x_left - 3:x_right, :) * 10000)]);
%
%     TVAC_cum = table([titleVAC; ...
%                     YearMonthWeekJP(x_left - 3:x_right), ...
%                     round(BackData_Area_cv_all(x_left - 3:x_right, :) * 100000000)]);
%
%     writetable(TN, [home 'Figures/' char(PREF) '/' char(figfolder) '/BackData_' char(figname_main) char(PREF) '.xls'], 'Sheet', '新規感染者数（1日平均）', 'WriteVariableNames', false);
%     writetable(TVAC_new, [home 'Figures/' char(PREF) '/' char(figfolder) '/BackData_' char(figname_main) char(PREF) '.xls'], 'Sheet', 'ワクチン新規接種パス', 'WriteVariableNames', false);
%     writetable(TVAC_cum, [home 'Figures/' char(PREF) '/' char(figfolder) '/BackData_' char(figname_main) char(PREF) '.xls'], 'Sheet', 'ワクチン累計接種パス', 'WriteVariableNames', false);
%     writetable(TAD, [home 'Figures/' char(PREF) '/' char(figfolder) '/BackData_' char(figname_main) char(PREF) '.xls'], 'Sheet', '経済損失と死亡者数', 'WriteVariableNames', false);
% end
%

%% Local Functions
function indexDate = findDateIndex(date, firstDate, lastDate)
% find index in date lies between "firstDate" and "secondDate"
    indexDate = find( (datenum(date)  >= datenum(firstDate)) ...
                     & (datenum(date) <= datenum(lastDate))    );
    indexDate = indexDate(indexDate>=1);
end

function simSeries = adjustment_data_to_simulation(pastSeries, simSeries, transitionLength)
    Y = pastSeries(end-1:end);
    X = [1, 1; 1, 2];
    b = (X'*X)^(-1)*X'*Y;
    
    initial_trend =  b(1) + b(2) * 3;
    
%     initial = (weight * pastSeries(end)/simSeries(1) + initial_trend/simSeries(1) * (1-weight));
    initial = initial_trend/simSeries(1);
    
    transitionPath = linspace(initial,1,transitionLength);
    simSeries(1:transitionLength)  = simSeries(1:transitionLength)  .* transitionPath(1:end)';

end