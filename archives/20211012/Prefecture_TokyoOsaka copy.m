% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata
%

clear variables
close all
iPC = 0; % 0 for Mac, 1 for Windows
if iPC == 1
    home = '\Users\masam\Dropbox\fujii_nakata\Policy_request\2021年7月\BSFuji_2021JUL23\';
else
%     home = '/Users/sohtakawawaki/Dropbox/fujii_nakata (1)/Website/Codes/';
    home = '/Users/ymaeda/Dropbox/fujii_nakata/Policy_request/2021年7月/BSFuji_2021JUL23/';
    %     home = '/Users/koheimachi/Dropbox/fujii_nakata/Website/Codes/';
end
cd(home);
%====================== Program parameter values ======================%
pindex = 1;
figure_save = 0; % 0 = figures won't be saved, 1 = they will be saved
data_save = 0; % save back data
vaccine_figure_loop = 0; % =0 appear only once; =1 appear every loop;
beta_figure_loop = 0; % =0 appear only once; =1 appear every loop;
vaccine_disp_switch = 0; % =0 not display the summary of # of the vaccinated ; =1 display
% ICU_nation = 1; % = 1 use national definition (NHK data), = 0 use data from Tokyo Keizai
% in the "Figure" folder
fs = 12; % common font size for many figures
ldfs = 12; % legend font size for vaccine path
ft = '%.1f';
if iPC == 1
    fn = 'Yu Gothic'; % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
linecolor = {'blue', 'black', 'blue','red', 'black', 'red'};
language = {'EN','JP'};
%======================================================================%

%================== Model Fixed Parameter Values ============================%
parameter

%================ Parameter Values (Prefecture Specific) ====================%
prefecture_parameter

%--- Import data ---%
import_prefecture
% Covid data are recorded at weekly frequencies (Mon-Sun)
% The first week start on January 20 (Mon), 2020

%--- Construct weekly vaccine data ---%

[V1_medical_ori, V1_medical, V2_medical_ori, V2_medical,...
    V1_elderly_ori, V1_elderly, V2_elderly_ori, V2_elderly, ...
    V1_others_ori, V1_others, V2_others_ori, V2_others,...
    vs_MT,vs_ET,vs_M1,vs_E1,vs_M2,vs_E2] ...
    = ImportVaccineData(home,iPC,Data,pref,dateEN,ps,vaccine_disp_switch);
%
V1_medical(end) = 0;
V2_medical(end) = 22708;
V1_elderly(end) = 250000;
V2_elderly(end) = 3.4784e+05;
V2_others(end) = 2.8962e+05;
V1_others(end) = 8400000*ps - V2_medical(end)-V1_elderly(end)-V2_elderly(end)-V2_others(end);

% [V1_medical_ori, V1_medical, V2_medical_ori, V2_medical,...
%     V1_elderly_ori, V1_elderly, V2_elderly_ori, V2_elderly, ...
%     vs_MT,vs_ET,vs_M1,vs_E1,vs_M2,vs_E2] ...
%     = ImportVaccineData(home,iPC,Data,pref,dateEN,ps,vaccine_disp_switch);

%--- Constructing the reference level of output ---%
[potentialGDP, referenceGDP, alpha] = construct_GDP(GDP,TdataGDP);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alpha(end-1)=alpha(end-2);%
% alpha(end)=alpha(end-2);%set alpha level for the first 2 weeks of July at the level of the last week of June
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%--- Regress mobility on alpha to estimate the elasticity h ---%
[Malt,h_all_ori,h_all_se,h_ori,h_se] = estimate_h(M,alpha,TdataGDP,RetroH,hconstant);

%--- Plot mobility data ---%
figname = string(['Mobility_GDP_' char(pref)]);
f = figure('Name',figname);
plot_mobility(Malt,alpha,Tdata,TdataGDP,MonthWeekJP,xtick1,fs,16)
if figure_save == 1
    %saveas(figure(2),[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
    saveas(f,[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
end

%--- Import ICU data ---%
ICU_nation = zeros(Tdata+1,1);
ICU_pref = zeros(Tdata+1,1);
BED = zeros(Tdata,1);
% if ICU_nation == 1
ICU_nation(2:Tdata+1,1) = Data(:,22);
BED(1:Tdata,1) = Data(:,23);
% else
ICU_pref(2:Tdata+1,1) = Data(:,21);
% end

%--- Import hospitalized patient data ---%
hospital = zeros(Tdata+1,1);
hospital(2:end) = Data(:,25);
hospital(isnan(hospital)) = 0;

%--- Plot ICU data ---%
figname = string(['ICU_nation_transition_' char(pref)]);
f = figure('Name',figname);
plot(ICU_nation, 'LineWidth', 1.5)
title('Transition of  ICU (National Definition)')
ytickformat('%,6.0f')
xticks(find(WeekNumber==1))
xticklabels(MonthWeekJP(WeekNumber==1))
lgd.NumColumns = 2;
xtickangle(45)

figname = string(['ICU_prefecture_transition_' char(pref)]);
f = figure('Name',figname);
plot(ICU_pref, 'LineWidth', 1.5)
title('Transition of  ICU (Prefecture-Specific Definition)')
ytickformat('%,6.0f')
xticks(find(WeekNumber==1))
xticklabels(MonthWeekJP(WeekNumber==1))
lgd.NumColumns = 2;
xtickangle(45)


%--- Compute the history of S, I, R, D in the data period ---%
% [S,I,R,D]...
%     = SIRD(Tdata,POP0,N,E1,E2,...
%     V1_elderly,V1_medical,V2_elderly,V2_medical,...
%     gamma,dD,TdataGDP,referenceGDP,alpha);
[S,I,R,D]...
    = SIRD(Tdata,POP0,N,E1,E2,...
    V1_elderly,V1_medical,V1_others,V2_elderly,V2_medical,V2_others,...
    gamma,dD,TdataGDP,referenceGDP,alpha);
InitialValues = [S(end),I(end),R(end),D(end),ICU_nation(end),ICU_pref(end),hospital(end)];





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Main analysis starts here %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ======================================================= %
% ==================== CHANGE HERE ====================== %
% ======================================================= %

% TH = [-0.5,0,0.5]; %[-0.5,0,0.5,-0.5,0,0.5];
TH = [-2,0,2];
SoE_TH_off=[1100,1400,1490]; %[1100,1400,1490,1100,1400,1490];
th_off_date_vec =  [5,5,5]; % [6,6,6,8,8,8];


TH_index = TH;

xline_ind = find(date == datetime(2021,6,17)); % when the state of emergency end

alpha_May = mean(alpha((dateEN >= datetime(2020, 5, 07)) & (datetime(2020, 5, 28) >= dateEN)));
alpha_Jan = mean(alpha((dateEN >= datetime(2021, 1, 07)) & (datetime(2021, 1, 28) >= dateEN)));
alpha_Feb = mean(alpha((dateEN >= datetime(2020, 2, 7)) & (datetime(2020, 2, 28) >= dateEN))); % output loss without the state of emergency
alpha_Nov = mean(alpha((dateEN >= datetime(2020, 11, 5)) & (datetime(2020, 11, 26) >= dateEN))); % output loss without the state of emergency

alpha_Jan_2020 = mean(alpha((dateEN >= datetime(2020, 1, 1)) & (datetime(2020, 1, 31) >= dateEN)));

ind_date = ind_date - 1;

% alpha_off = alpha_Feb; %0.5*alpha_Feb + 0.5*alpha_Nov;
alpha_off = alpha_Jan_2020; %0.5*alpha_Feb + 0.5*alpha_Nov;

if pindex == 1
    
    state = 1;
    %alpha_scale = 0.94;
    beta_scale = 0.85; %0.85 ... optimisitc
    alpha_scale = 1.0;
    alpha_on = alpha(end); %alpha_Jan;%(0.5 * alpha_May + 0.5 * alpha_Jan); %alpha_Jan; %alpha_scale * (0.5 * alpha_May + 0.5 * alpha_Jan);
    h_scale = 1;
    %     beta_shock_after_emergency = 0.0;
    beta_shock_after_emergency = 0.075;
    rho_after_emergency = 0.95;
    %alpha_jump = 0;
    alpha_jump = (alpha(end)-alpha_on)/(alpha_off - alpha_on);
    %     ind_date = ind_date - 1;
    % paces_ori = 9800000;
    % gradual_paces = 3;
    % VT3share = 0.55;
    % ind_date = find(date == datetime(2021,7,1));
elseif pindex == 5
    state = 1;
    alpha_on = 1 * alpha_May;
    h_scale = 1;
    beta_shock_after_emergency = 0.7;
    % beta_shock_after_emergency = 0;
    rho_after_emergency = 0.75;
    %alpha_jump = 0;
    alpha_jump = (alpha(end)-alpha_on)/(alpha_off - alpha_on);
end
h = h_ori;
h_all = h_all_ori;
h(2) = h_scale*h_ori(2);


% ======================================================= %
% ======================================================= %
% for iAlpha = 1:length(alpha_on_vector_sim)
%     alpha_on = alpha_on_vector_sim(iAlpha)
nTH = length(TH);
DMat = nan(1,length(TH));
AlphaMat = nan(1,length(TH));
SimData = nan(SimPeriod+1,length(InitialValues),length(TH));
AlphaPath = nan(SimPeriod,length(TH));
NPath = nan(SimPeriod,length(TH));
SimERN = nan(SimPeriod,length(TH));
THonPath = nan(SimPeriod,length(TH));
BackDataN = zeros(SimPeriod+8,length(TH_index));
BackDataAlpha = zeros(SimPeriod+8,length(TH_index));
BackDataERN = zeros(SimPeriod+8,length(TH_index));
BackDataDA = zeros(length(TH),3);
betaPath = zeros(SimPeriod,nTH);
betaAvgVec = zeros(nTH,1);
deltaAvgVec = zeros(nTH,1);
betaAlphaPath = zeros(SimPeriod,nTH);
betaTildePath = zeros(SimPeriod,nTH);
deltaPath = zeros(SimPeriod,nTH);
deltaICUPath = zeros(SimPeriod,nTH);
varSharePath = zeros(SimPeriod,nTH);
varDeltaSharePath = zeros(SimPeriod,nTH);
Sim_dDPath=zeros(SimPeriod,nTH);
SimICU_nation= zeros(SimPeriod+1,nTH);
SimICU_pref = zeros(SimPeriod+1,nTH);
SimHospital = zeros(SimPeriod+1,nTH);

for iTH = 1:length(TH)
    %%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%
    %var_infection2 = TH(iTH);
    %beta_scale=TH;
    th_off1=SoE_TH_off(iTH)*7;
    h = h_ori + h_se*TH(iTH);
    h_all = h_all_ori + h_all_se*TH(iTH);
    th_off_date = th_off_date_vec(iTH);
    %%%%%%%%%%%%%%%%% Figure Titles %%%%%%%%%%%%%%%%%%%%%%%%%
    figname_main = 'ERN_MainResults_';
    figname_beta = 'BetaPath';
    figname_beta_loop = ['BetaPath' '_' sprintf(ft, th_off1 / 7)];
    figname_beta_tilde_loop = ['BetaTildePath_', sprintf(ft, th_off1 / 7)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [delta,beta_tilde,ERN,beta,ICU_inflow_nation, ICU_inflow_pref, Hospital_inflow,...
        gammaT,delta_average,delta_ICU_average,ICU_inflow_avg_nation,ICU_inflow_avg_pref,Hospital_inflow_avg,delta_sample,beta_avg,beta_se,delta_se]...
        = Time_Series_Average(S,I,D,ICU_nation,ICU_pref,hospital,dD,N,Tdata,SimPeriod,...
        RetroPeriod,POP0,gamma,hconstant,h_all,alpha,k,...
        gamma_ICU,ICU_adjustment,gamma_Hospital,Hospital_adjustment,RetroPeriodDelta,RetroPeriodICU);
    
    
    %%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%
    delta_past_avg = delta_average; %Past 17 weeks average
    beta_bar = beta_avg; %Past weeks average
    beta_avg = beta_bar + beta_se*TH(iTH);
    delta_average = delta_past_avg + delta_se*TH(iTH);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elderly_total = ps*elderly_jp;
    medical_total = ps*medical_jp;
    ordinary_total = ps*ordinary_jp;
    %medical = medical_total*accept_share;
    medical = medical_total;
    elderly = elderly_total*accept_share;
    ordinary = ordinary_total*accept_share_ordinary;
    %         elderly = elderly - (sum(V1_elderly));
    %--- Eliminate the effects of vaccination from delta ---%
    delta_ss = delta_average*(0.1063/1.53); %Share of the death rate among youth to the death rate of all populaiton
    VD_elderly = D1*V1_elderly + (D2-D1)*V2_elderly;
    VD_ordinary = (D1*V1_medical + (D2-D1)*V2_medical) + (D1*V1_others + (D2-D1)*V2_others);
    share = ((1 - sum(VD_elderly(1:end-2))/elderly_total) * (delta_average - delta_ss) ...
        + (1-sum(VD_ordinary(1:end-2))/(ordinary_total))*delta_ss)/delta_average;
    delta_wo_vaccination = delta_average/share; %Eliminate the effects of vaccination in the past average value, Past 17 without vaccine effects
    
    %--- Eliminate the effects of vaccination from delta ---%
    ICU_ss = delta_ICU_average*(0.3916/1.62); %Share of the death rate among youth to the death rate of all populaiton
    share = ((1 - sum(VD_elderly(1:end-2))/elderly_total) * (delta_ICU_average - ICU_ss) ...
        + (1-sum(VD_ordinary(1:end-2))/(ordinary_total))*ICU_ss)/delta_ICU_average;
    ICU_wo_vaccination = delta_ICU_average/share; %Eliminate the effects of vaccination in the past average value
    
    
    
    %--- Eliminate the effects of alpha variant from delta and beta---%
    [var_share,var_prev,var_initial] = var_share_prev(Data(:,20),SimPeriod,var_ss,var_growth);
    var_infection_adjustment = (1-mean(var_prev(end-RetroPeriod+1:end)))*1 ...
        + mean(var_prev(end-RetroPeriod+1:end))*(1 + var_infection); %Relative increase of infectiousness (alpha varaint, past 17 weeks)
    var_infection_delta_adjustment = (1-mean(var_prev(end-RetroPeriodDelta+1:end)))*1 ...
        + mean(var_prev(end-RetroPeriodDelta+1:end))*(1 + var_infection_delta); %Relative increase of death rate (alpha varaint, past 17 weeks)
    var_infection_ICU_adjustment = (1-mean(var_prev(end-RetroPeriodICU+1:end)))*1 ...
        + mean(var_prev(end-RetroPeriodICU+1:end))*(1 + var_infection_delta); %Relative increase of death rate (alpha varaint, past 17 weeks)
    
    delta_wo_alpha = delta_wo_vaccination/var_infection_delta_adjustment; %Eliminate the effects of alpha variant in the past average value
    
    ICU_wo_alpha = ICU_wo_vaccination/var_infection_ICU_adjustment;
    
    beta_wo_alpha = beta_avg/var_infection_adjustment; %Eliminate the effects of alpha variatn
    
    %--- Eliminate the effects of delta variant from delta and beta---%
    [var_share2,var_prev2,var_initial2] = var_share_prev(Data(:,24),SimPeriod,var_ss,var_growth2);
    var_infection_adjustment2 = (1-mean(var_prev2(end-RetroPeriod+1:end)))*1 ...
        + mean(var_prev2(end-RetroPeriod+1:end))*(1+var_infection2); %Relative increase of infectiousness (alpha varaint, past 17 weeks)
    var_infection_delta_adjustment2 = (1-mean(var_prev2(end-RetroPeriodDelta+1:end)))*1 ...
        + mean(var_prev2(end-RetroPeriodDelta+1:end))*(1 + var_infection_delta2); %Relative increase of death rate (alpha varaint, past 17 weeks)
    var_infection_ICU_adjustment2 = (1-mean(var_prev2(end-RetroPeriodICU+1:end)))*1 ...
        + mean(var_prev2(end-RetroPeriodICU+1:end))*(1 + var_infection_delta2); %Relative increase of death rate (alpha varaint, past 17 weeks)
    delta_average = delta_wo_alpha/var_infection_delta_adjustment2; %Eliminate the effects of delta variant in the past average value
    deltaAvgVec(iTH) = delta_average;
    
    beta_avg = beta_wo_alpha/var_infection_adjustment2; %Eliminate the effects of delta variant
    betaAvgVec(iTH) = beta_avg;
    
    ICU_average = ICU_wo_alpha/var_infection_ICU_adjustment2;
    
    %--- Construct vaccine distribution and delta path---%
    % [V,deltaT,VT] = ...
    %     vaccine_distribution_medical(V1_medical,V2_medical,...
    %     V1_elderly,V2_elderly,...
    %     elderly_jp,medical_jp,ordinary_jp,accept_share,...
    %     delta_average,lag,medical_duration,...
    %     paces_ori,sw_vacpath,gradual_paces,...
    %     E1,E2,D1,D2,ps,POP0,SimPeriod);
    
    [V,deltaT,VT,delta_ICU] = ...
        vaccine_distribution(V1_medical,V2_medical,...
        V1_elderly,V2_elderly,V1_others,V2_others,...
        ind_date,ind_date2,date_slowdown,lagged_VT2share,VT3share,...
        elderly,elderly_total,medical,ordinary,delta_average,ICU_average,lag,...
        medical_duration,paces_ori,paces2,paces3,sw_vacpath,gradual_paces,gradual_paces2,...
        E1,E2,D1,D2,ps,POP0,SimPeriod,Tdata);
    
    %Assuming that alpha variant = 100% without delta variant
    deltaT = deltaT*(1+var_infection_delta); %Alpha variant adjsutment for death rate
    delta_ICU = delta_ICU*(1+var_infection_delta); %Alpha variant adjsutment for ICU rate
    
    %Delta variant adjustment
    deltaT = deltaT.*(1+var_infection_delta2*var_share2);
    delta_ICU = delta_ICU.*(1+var_infection_delta2*var_share2);
    %     [deltaT,delta_bar] = construct_delta_variant(RetroPeriodDelta,...
    %         retro_lb,retro_ub,deltaT,delta,delta_average,var_prev,var_share,var_infection_delta,I);
    %     [delta_ICU,delta_ICU_bar] = construct_delta_variant(RetroPeriodDelta,...
    %         retro_lb,retro_ub,delta_ICU,delta,delta_average,var_prev,var_share,var_infection_delta,I);
    
    
    % figure for vaccine path
    if vaccine_figure_loop == 0
        if iTH == 1
            plot_vaccinepath_ym(200,VT,V1_medical,V2_medical,V1_elderly,V2_elderly,V1_others,V2_others,date,ps,MonthWeekJP,WeekNumber,Tdata,fs,ldfs,fn);
            %             plot_deltapath(201,delta,deltaT,deltaT(1),MonthWeekJP,WeekNumber,Tdata,fs,fn,iTH);
        end
        plot_deltapath(240,delta,deltaT,deltaT(1),MonthWeekJP,WeekNumber,Tdata,fs,fn,iTH);
        title("致死率（現在のレベルで標準化）", 'FontSize',fs,'FontName',fn);
        plot_serious_ratepath(250,delta,deltaICUPath(:,iTH),deltaICUPath(1,iTH),MonthWeekJP,WeekNumber,Tdata,fs,fn,iTH);
        title("重症化率（現在のレベルで標準化）", 'FontSize',fs,'FontName',fn);
        if figure_save == 1 && iTH == nTH
            saveas(figure(200), [home 'Figures/' char(pref)  '/Vaccine_Path.png']);
        end
        if figure_save == 1 && iTH == nTH
            saveas(figure(240), [home 'Figures/' char(pref)  '/Death_rate.png']);
        end
        %         if figure_save == 1 && iTH == 3
        %             saveas(figure(250), [home 'Figures/' char(pref)  '/ICU_rate_DR' '_'  sprintf('%.0f', paces_ori),   '.png']);
        %         end
    elseif vaccine_figure_loop == 1
        plot_vaccinepath_ym(200+iTH,VT,V1_medical,V2_medical,V1_elderly,V2_elderly,V1_others,V2_others,date,ps,MonthWeekJP,WeekNumber,Tdata,fs,ldfs,fn);
        subplot(1,2,1)
        title(string(['新規ワクチン接種本数（週ごと）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
        subplot(1,2,2)
        title(string(['累計ワクチン接種本数（週ごと）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
    end
    
    %     if figure_save == 1 && iTH == 1
    %         saveas(figure(200), [home 'Figures/' char(pref) '/Vaccine_path_' num2str(paces_ori) '.png']);
    %     end
    %
    %     if figure_save == 1 && iTH == 3
    %         saveas(figure(240), [home 'Figures/' char(pref)  '/Death_rate_DR'  '.png']);
    %         saveas(figure(250), [home 'Figures/' char(pref)   '/ICU_rate_DR'  '.png']);
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     beta_r = 0;
    %     for retrop = retro_lb:retro_ub
    %         beta_r = beta_r + mean(beta(end-retrop+1:end));
    %     end
    %     beta_avg = beta_r/(retro_ub-retro_lb+1);
    %     [betaT,betaT_woAR1,beta_bar] = construct_beta(SimPeriod,retro_lb,retro_ub,...
    %         beta,beta_avg,var_prev,var_share,var_infection);
    
    %%%%  デルタ株　　　%%%%
    %%%%%%%%%%%%%%%%%%%%%​
    %     var_intercept2 = log(var_initial2/(var_ss-var_initial2));
    %     var_share2 = exp((1:SimPeriod)'*var_growth2+var_intercept2).*var_ss./(1+exp((1:SimPeriod)'*var_growth2+var_intercept2));
    %     %     var_power = 1+var_share2*var_growth2;
    %     %     var_boost = ones(SimPeriod,1);
    %     %     var_boost(var_start:end) = var_power(1:SimPeriod-var_start+1);
    %     var_shareInd = zeros(SimPeriod,1);
    %     var_shareInd(var_start,1) = var_initial2;
    %     var_shareInd(var_start+1:end,1) = var_share2(1:SimPeriod-var_start);
    %     varDeltaSharePath(:,iTH) = var_shareInd;
    %     betaT_Eng = betaT;
    %     betaT_Eng = beta_AR1(betaT_temp_ini, beta_rho, betaT_Eng, start_beta);
    %     betaT = betaT.*(1+var_shareInd*var_infection2);
    %     betaT = beta_AR1(betaT_temp_ini, beta_rho, betaT, start_beta);
    
    %Assuming that alpha variant = 100% without delta variant
    beta_Eng = beta_avg*(1+var_infection); %Alpha variant adjsutment for beta
    
    %Delta variant adjustment
    betaT = beta_Eng.*(1+var_infection2.*var_share2);
    
    %AR1 adjustment for betaT
    betaT_woAR = betaT;
    betaT = beta_AR1(betaT_temp_ini, beta_rho, betaT, start_beta);
    if iTH < 4 
        betaT(1)=betaT(1)*1.47;
        betaT(2)=betaT(2)*1.23;
        betaT(3:4)=beta_scale*betaT(3:4);
    else
        betaT(1)=betaT(1)*1.26; %1.32;
        betaT(2:4)=beta_scale*betaT(2:4);
    end
    
    %AR1 adjustment for deltaT
    deltaT_woAR = deltaT;
    deltaT = beta_AR1(deltaT_temp_ini, delta_rho, deltaT, start_delta);
%     
    %AR1 adjustment for delta_ICU
%     delta_ICU_woAR = delta_ICU;
%     delta_ICU = beta_AR1(deltaT_temp_ini, delta_rho, delta_ICU, start_delta);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [DMat(iTH),AlphaMat(iTH),AlphaPath(:,iTH),SimData(:,:,iTH),NPath(:,iTH),SimERN(:,iTH),THonPath(:,iTH),SimICU_nation(:,iTH),SimICU_pref(:,iTH),SimHospital(:,iTH),betaShock] ...
        = Covid_projection_pref_date...
        (InitialValues,alpha_on,alpha_off,...
        th_on1,th_on2,th_off1,th_off2,th_off3,...
        betaT,gammaT,deltaT,deltaT(1),delta_ICU,V,h,k,POP0,hconstant,...
        DRi,state,ICU_inflow_avg_nation,ICU_inflow_avg_pref,gamma_ICU,...
        Hospital_inflow_avg,gamma_Hospital,beta_shock_after_emergency,rho_after_emergency,alpha_jump,th_off_date);
    
    % = Covid_projection_ICU4...
    % (InitialValues,alpha_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,betaT,gammaT,deltaT,deltaT(1),V,h,k,POP0,hconstant,DRi,state,ICU_inflow_avg,gamma_ICU,beta_shock_after_emergency,rho_after_emergency);
    betaPath(:,iTH) = betaShock;
    betaShock_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaShock;
    betaTildePath(:,iTH) = betaShock_tilde;
    
    deltaPath(:,iTH) = deltaT;
    deltaICUPath(:,iTH) = delta_ICU;
    varSharePath(:,iTH) = var_share2;
    
    Sim_dDPath(:,iTH)=SimData(2:end,4,iTH)-SimData(1:end-1,4,iTH);
    %     % Plot betaT
    %     if beta_figure_loop == 0
    %         if iTH == 1
    %             %figure(400) % Figure for BetaT with Variant Share
    %             f=figure('Name',char(figname_beta));
    %             set(gcf,'Position',[100,100,1200,500])
    %             %plot_beta_India(var_initial,var_initial2,var_share,var_shareInd,beta,beta_avg,betaT,betaT_Eng,betaT_woAR1,dateD,SimDate,Tdata,MonthWeekJP,MonthNumber,WeekNumber,fn)
    %             plot_beta_India( var_initial,var_initial2,var_share,var_shareInd,beta,beta_avg,betaShock,betaT_Eng,betaT,  dateD,SimDate,Tdata,MonthWeekJP,MonthNumber,WeekNumber,fn)
    %             if figure_save == 1
    %                 saveas(f,[home 'Figures/' char(pref) '/' char(figname_beta) '.png']);
    %             end
    %             betaT_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT;
    %             betaT_woAR1_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT_woAR1;
    %             beta_tilde_avg = beta_avg*((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k);
    %
    %             %figure(401)
    %             figure('Name','BetaTildePath')
    %             set(gcf,'Position',[100,100,1200,500])
    %             plot_beta(var_initial,var_share,beta_tilde,beta_tilde_avg,betaShock_tilde,betaT_tilde,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
    %             title('β tildeの推移','FontSize',20,'FontWeight','normal','FontName',fn)
    %         end
    %     elseif beta_figure_loop == 1
    %         %figure(400+iTH) % Figure for BetaT with Variant Share
    %         figure('Name',figname_beta_loop)
    %         set(gcf,'Position',[100,100,1200,500])
    %         plot_beta_India( var_initial,var_initial2,var_share,var_shareInd,beta,beta_avg,betaShock,betaT_Eng,betaT,  dateD,SimDate,Tdata,MonthWeekJP,MonthNumber,WeekNumber,fn)
    %         %plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaT,betaT_woAR1,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
    %         %plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaShock,betaT,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
    %         title(string(['βの推移（',sprintf(ft,TH(iTH)),'）']),'FontSize',20,'FontWeight','normal','FontName',fn)
    %         betaT_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT;
    %         betaShock_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaShock;
    %         betaT_woAR1_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT_woAR1;
    %         beta_tilde_avg = beta_avg*((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k);
    %
    %         %figure(440 + iTH)
    %
    %         figure('Name',figname_beta_tilde_loop)
    %         set(gcf,'Position',[100,100,1200,500])
    %         plot_beta(var_initial,var_share,beta_tilde,beta_tilde_avg,betaShock_tilde,betaT_tilde,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
    %         title(string(['β tildeの推移（',sprintf(ft,TH(iTH)),'）']),'FontSize',20,'FontWeight','normal','FontName',fn)
    %     end
    
    
end

%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TH = [600,800,1000];
%TH_index = [600,800,1000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%
%minAlpha = min(minAlphaMat); %minimum alpha when variants have no effects.
minAlpha = alpha_off; % 経済損失0 = 2020年2月のGDP level

AlphaM = AlphaMat(~isnan(AlphaMat));
AlphaM = (AlphaM - minAlpha)*prefGDP*10000;
DM = DMat(~isnan(DMat));
BackDataDA(1:length(TH),:) = [round(AlphaM'),round(DM'),round(TH',1)];
%--- Record how many times on and off are triggered ---%
waves = zeros(1,length(TH));
for i = 1:length(TH)
    svec = zeros(SimPeriod-1,1);
    for t = 1:SimPeriod-1
        svec(t) = AlphaPath(t+1,i)-AlphaPath(t,i);
    end
    waves(i) = sum(svec>0);
end
%%
for l = 1:2 %1:2 when english version needed
    % Generate graphs for the website
    lng = language{l};
    %figname = 100 + l;
    %figure(figname);
    f = figure('Name',[char(figname_main) char(lng)]);
    if iPC == 1
        f.WindowState = 'maximized';
    else
        set(gcf,'Position',[100,100,1000,1000])
    end
    %--- Number of people who are in ICU ---%
    subplot(3, 3, 1)
    BackDataICU_pref = ...
        plot_ICU(TH,TH_index,ICU_pref,SimICU_pref,ICU_limit,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,fs,fn,ft,l);
    if l == 1
        title('ICU (Tokyo Standard)','FontSize',fs,'FontWeight','normal')
        xticklabels(MonthWeekEN(WeekNumber==1))
        %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
    elseif l == 2
        title('重症患者数（都基準）','FontSize',fs,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP(WeekNumber==1))
    end
    xlim([63, Tdata+24])
    ylim([0 ICU_limit*1.1])
    
    
    subplot(3, 3, 2)
    %     [BackDataICU_nation, BackDataICU_pref] = ...
    %         plot_ICU_both(TH,TH_index,ICU_nation,SimICU_nation,BED(Tdata),...
    %         ICU_pref,SimICU_pref,ICU_limit,MonthWeekEN,MonthWeekJP,WeekNumber,...
    %         Tdata,linecolor,fs,fn,ft,l);
    BackDataICU_nation=...
        plot_ICU(TH,TH_index,ICU_nation,SimICU_nation,BED(Tdata),MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,fs,fn,ft,l);
    if l == 1
        title('ICU (National Standard)','FontSize',fs,'FontWeight','normal')
        xticklabels(MonthWeekEN(WeekNumber==1))
        %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
    elseif l == 2
        title('重症患者数（国基準）','FontSize',fs,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP(WeekNumber==1))
    end
    xlim([63, Tdata+24])
    
    lineName = cell(nTH,1);
    for i = 1:nTH
        if l == 1
            lineName{i} = ['Relative Infection Rate of Delta Variants = ' sprintf(ft,1 + TH_index(i))];
        elseif l == 2
            lineName{i} = ['デルタ株の相対感染力:' sprintf(ft,1 + TH_index(i)) '倍'];
        end
    end
    lineStyle = {'-','-','-','-','-','-'};
    lineWidth = [0.5,1.5,0.5,0.5,1.5,0.5];
    lgfs = 12;
    show_other = 0; % Do not show lines other than TH index
    yft = '%.0f';
    column = 1;
    
    %--- Number of newly hospitalized ---%
    subplot(3,3,3)
    eTitle = 'Hospitalized Patients';
    jTitle = '入院患者数';
    plot_function(TH,TH_index,hospital(2:end),SimHospital(2:end,:),...
        MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineWidth,lineName,...
        eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'Northwest')
    %round(Sim_dDPath/7,0)
    xline(83,'LineWidth',1.5,'HandleVisibility','off');
    xlim([63, Tdata+24])
    %     ylim([0 20])
    Handle=legend;
    set(Handle, 'Visible','off');
    yline(Hospital_limit/2,'--k',"50%",'LineWidth',1.0,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
    yline(Hospital_limit,'--k',"100%",'LineWidth',1.0,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
    ylim([0 Hospital_limit*1.1])
    
    subplot(3,3,4)
    eTitle = 'New Deaths (Daily Average)';
    jTitle = '新規死亡者数（1日平均）';
    plot_function(TH,TH_index,dD/7,Sim_dDPath/7,...
        MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineWidth,lineName,...
        eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'Northwest')
    %round(Sim_dDPath/7,0)
    xline(83,'LineWidth',1.5,'HandleVisibility','off');
    xlim([63, Tdata+24])
    ylim([0 20])
    Handle=legend;
    set(Handle, 'Visible','off');
    
    %--- Number of new cases ---%
    eTitle = 'New Cases (Daily Average)';
    jTitle = '新規感染者数（1日平均）';
    subplot(3,3,5)
    plot_function(TH,TH_index,N/7,NPath/7,...
        MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineWidth,lineName,...
        eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'Southwest')
    Handle=legend;
    set(Handle, 'Visible','off');
    xline(83,'LineWidth',1.5,'HandleVisibility','off');
    xticks(find(WeekNumber==1))
    xtickangle(45)
    xlim([63, Tdata+24])
    ylim([0, 3500])
    
    %--- GDP Path ---%
    subplot(3,3,7)
    MonthWeek = [MonthWeekEN,MonthWeekJP];
    plot_Alpha(alpha,AlphaPath,TH,TH_index,MonthWeek(:,l),WeekNumber,MonthNumber,Tdata,{'blue','blue','blue','red','red','red'},ft,fs,fn,l)
    title('GDP','FontSize',fs,'FontWeight','normal')
    xline(83,'LineWidth',1.5,'HandleVisibility','off');
    
    %--- Trade-off Curve ---%
    subplot(3, 3, 8)
    plot_Tradeoff2(AlphaM, DM, waves, TH, TH_index, l, linecolor, fs, fn)
    
    
    
    if figure_save == 1
        saveas(f,[home 'Figures/' char(pref) '/' char(figname_main) char(lng) '.png']);
    end
end %End of language loop = figure loop
%%
th_wave = 1;
for i = 1:nTH
    if abs(TH(i) - TH_index(th_wave)) < 0.0001
        BackDataN(:,th_wave) = [N(end-7:end);NPath(:,i)];
        BackDataAlpha(:,th_wave) = [alpha(end-7:end);AlphaPath(:,i)];
        BackDataERN(:,th_wave) = [ERN(end-7:end);SimERN(:,i)];
        BackDatadD(:,th_wave) = [dD(end-7:end);Sim_dDPath(:,i)];
        BackDataHospital(:,th_wave) = [hospital(end-7:end);SimHospital(2:end,i)];
        th_wave = th_wave + 1;
        if th_wave > length(TH_index)
            th_wave = length(TH_index);
        end
    end
end

if data_save == 1
    titleN = strings(1,1+length(TH_index)*6);
    titleN(1) = "週";
    for ti = 1:length(TH_index)
        titleN(1,1+ti) = string(['新規感染者数（',sprintf(ft,TH_index(ti)),'）']);
        titleN(1,1+length(TH_index)+ti) = string(['経済活動（',sprintf(ft,TH_index(ti)),'）']);
        titleN(1,1+length(TH_index)*2+ti) = string(['実効再生産数（',sprintf(ft,TH_index(ti)),'）']);
        titleN(1,1+length(TH_index)*3+ti) = string(['重症者数_国基準（',sprintf(ft,TH_index(ti)),'）']);
        titleN(1,1+length(TH_index)*4+ti) = string(['重症者数_都基準（',sprintf(ft,TH_index(ti)),'）']);
        titleN(1,1+length(TH_index)*5+ti) = string(['新規死亡者数（',sprintf(ft,TH_index(ti)),'）']);
        titleN(1,1+length(TH_index)*6+ti) = string(['入院患者数（',sprintf(ft,TH_index(ti)),'）']);
    end
    TN = table([titleN;MonthWeekJP(Tdata-7:end-1),round(BackDataN(:,1:length(TH_index))/7),...
        round(100*(1-BackDataAlpha(:,1:length(TH_index))),1),round(BackDataERN(:,1:length(TH_index)),2),...
        round(BackDataICU_nation(:,1:length(TH_index))),round(BackDataICU_pref(:,1:length(TH_index))),...
        round(BackDatadD(:,1:length(TH_index))/7),round(BackDataHospital(:,1:length(TH_index)))]);
    titleAD = ["経済損失（億円）","死亡者数","ケース"];
    TAD = table([titleAD;BackDataDA(1:length(TH),:)]);
    writetable(TN,[home 'Figures/' char(pref) '/BackData_' char(figname_main)  char(pref) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
    writetable(TAD,[home 'Figures/' char(pref) '/BackData_' char(figname_main)  char(pref) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
end


%% Other Plot
lgfs = 20;
column = 3;
l = 2; %Japanese
show_other = 0; % Do not show lines other than TH index
yft = '%.0f';

% Transitions of N and I
figname_NI = 'NI_transitions';
lineName = cell(nTH,1);
lineName2 = cell(nTH,1);
for i = 1:nTH
    lineName{i} = ['I, 宣言解除 : ' sprintf(ft, TH(i))];
    lineName2{i} = ['N, 宣言解除 : ' sprintf(ft, TH(i))];
end
eTitle = 'Transitions of N and I';
jTitle = 'NとIの推移';
figure('Name',char(figname_NI));
set(gcf,'Position',[100,100,1200,800])
plot_function2(TH,TH_index,I(2:end),SimData(2:end,2,:),N/7,NPath/7,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineName,...
    lineName2,eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
xlim([Tdata-31 Tdata+41])

% Transitions of S and R
figname_SR = 'SR_transitions';
lineName = cell(nTH,1);
lineName2 = cell(nTH,1);
for i = 1:nTH
    lineName{i} = ['S, 宣言解除 :  ' sprintf(ft, TH(i))];
    lineName2{i} = ['R, 宣言解除 :  ' sprintf(ft, TH(i))];
end
eTitle = 'Transitions of S and R';
jTitle = 'SとRの推移';
figure('Name',char(figname_SR));
set(gcf,'Position',[100,100,1200,800])
plot_function2(TH,TH_index,S(2:end),SimData(2:end,1,:),R(2:end),SimData(2:end,3,:),...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineName,...
    lineName2,eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
xlim([Tdata-21 Tdata+41])

% Transitions of Beta
figname_beta = 'Beta transition';
figure('Name',char(figname_beta));
set(gcf,'Position',[100,100,1200,800])
lineName = cell(nTH,1);
for i = 1:nTH
    lineName{i} = ['宣言解除 :  ' sprintf(ft, TH(i))];
end
eTitle = 'Transitions of \beta';
jTitle = '\betaの推移';
yft = '%.3f';
plot_function(TH,TH_index,beta,betaPath,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineWidth,lineName,...
    eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
plot([NaN(Tdata-1,1);ones(SimPeriod+1,1)*beta_bar],'--k','LineWidth',1.5,'DisplayName','過去17週間平均')
plot([NaN(Tdata-1,1);ones(SimPeriod+1,1)*beta_wo_alpha],'--b','LineWidth',1.5,'DisplayName','過去17週間平均, アルファ株影響除去')
plot([NaN(Tdata-1,1);ones(SimPeriod+1,1)*beta_avg],'--r','LineWidth',1.5,'DisplayName','過去17週間平均, 変異株影響除去')
plot([NaN(Tdata,1);betaT_woAR],'-g','LineWidth',1.5,'DisplayName','Without AR(1) shock')
xline(Tdata-RetroPeriod-1,'HandleVisibility','off')
xlim([Tdata-31 Tdata+41])

% Transitions of Beta Tilde
figname_beta_tilde = 'Beta Tilde transition';
figure('Name',char(figname_beta_tilde));
set(gcf,'Position',[100,100,1200,800])
eTitle = 'Transitions of \beta tilde';
jTitle = '\beta tildeの推移';
yft = '%.3f';
plot_function(TH,TH_index,beta_tilde,betaTildePath,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineWidth,lineName,...
    eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
xline(Tdata-RetroPeriod-1,'HandleVisibility','off')
xlim([Tdata-31 Tdata+41])
if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref)  '/Beta_Tilde_Path.png']);
end


% Transitions of Variant Share Rate
figname_var = 'Variant Share Rate transition';
figure('Name',char(figname_var));
set(gcf,'Position',[100,100,1200,800])
eTitle = 'Transitions of Variant Rate';
jTitle = '変異株割合の推移';
plot_function(TH,TH_index,var_prev2,varSharePath,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineWidth,lineName,...
    eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
xlim([Tdata-31 Tdata+41])
if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref)  '/Variant_Share.png']);
end

% Transitions of Death Rate
column = 2;
figname_delta = 'Death Rate transition';
figure('Name',char(figname_delta));
set(gcf,'Position',[100,100,1200,800])
column = 2;
eTitle = 'Transitions of Death Rate';
jTitle = '死亡率の推移';
plot_function(TH, TH_index, delta, deltaPath, ...
    MonthWeekEN, MonthWeekJP, WeekNumber, Tdata, xline_ind, linecolor,lineWidth, lineName, ...
    eTitle, jTitle, fs, lgfs, fn, ft, yft, column, l, show_other, 'SouthOutside')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * delta_past_avg], '--k', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * delta_wo_vaccination], '--b', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均, ワクチン影響除去')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * delta_wo_alpha], ':r', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均, ワクチン影響除去+アルファ株影響除去')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * delta_average], '--r', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均, ワクチン影響除去+変異株影響除去')
xline(Tdata - RetroPeriodDelta - 1, 'HandleVisibility', 'off')
xlim([Tdata - 31 Tdata + 41])
if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref)  '/' char(figname_delta) '.png']);
end


% Transitions of ICU Rate
figname_ICU = 'ICU Rate transition Naitonal Standard';
figure('Name',char(figname_ICU));
set(gcf,'Position',[100,100,1200,800])
eTitle = 'Transitions of ICU Rate (Naitonal Standard)';
jTitle = '重症化率の推移(国基準)';
plot_function(TH,TH_index,delta.*ICU_inflow_nation,ICU_inflow_avg_nation*deltaICUPath,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineWidth,lineName,...
    eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * delta_ICU_average*ICU_inflow_avg_nation], '--k', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * ICU_wo_vaccination*ICU_inflow_avg_nation], '--b', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均, ワクチン影響除去')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * ICU_wo_alpha*ICU_inflow_avg_nation], ':r', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均, ワクチン影響除去+アルファ株影響除去')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * ICU_average*ICU_inflow_avg_nation], '--r', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均, ワクチン影響除去+変異株影響除去')
xline(Tdata-RetroPeriodICU-1,'HandleVisibility','off')
xlim([Tdata-31 Tdata+41])
if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref)  '/' char(figname_ICU) '.png']);
end

% Transitions of ICU Rate
figname_ICU = 'ICU Rate transition Local Standard';
figure('Name',char(figname_ICU));
set(gcf,'Position',[100,100,1200,800])
eTitle = 'Transitions of ICU Rate (Local Standard)';
jTitle = '重症化率の推移(地方基準)';
plot_function(TH,TH_index,delta.*ICU_inflow_pref,ICU_inflow_avg_pref*deltaPath,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineWidth,lineName,...
    eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * delta_past_avg*ICU_inflow_avg_pref], '--k', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * ICU_wo_vaccination*ICU_inflow_avg_pref], '--b', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均, ワクチン影響除去')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * ICU_wo_alpha*ICU_inflow_avg_pref], ':r', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均, ワクチン影響除去+アルファ株影響除去')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * ICU_average*ICU_inflow_avg_pref], '--r', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均, ワクチン影響除去+変異株影響除去')
xline(Tdata-RetroPeriodDelta-1,'HandleVisibility','off')
xlim([Tdata-31 Tdata+41])
if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref)  '/' char(figname_ICU) '.png']);
end


% Transitions of ERN
figname_ERN = 'ERN transition';
figure('Name',char(figname_ERN));
set(gcf,'Position',[100,100,1200,800])
eTitle = 'Transitions of ERN';
jTitle = '実効再生産数の推移';
plot_function(TH,TH_index,ERN,SimERN,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineWidth,lineName,...
    eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
xlim([Tdata-31 Tdata+41])
if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref)  '/' char(figname_ERN) '.png']);
end

% Transitions of New deaths
figname_dD = 'dD transition';
figure('Name',char(figname_dD));
set(gcf,'Position',[100,100,1200,800])
eTitle = 'Transitions of dD';
jTitle = '新期死亡者数の推移';
plot_function(TH,TH_index,dD,Sim_dDPath,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineWidth,lineName,...
    eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
xlim([Tdata-31 Tdata+41])
if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref)  '/' char(figname_dD) '.png']);
end
