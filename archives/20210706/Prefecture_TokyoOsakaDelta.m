% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata
%

clear variables
close all
iPC = 0; % 0 for Mac, 1 for Windows
if iPC == 1
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
    home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/';
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
ICU_nation = 1; % = 1 use national definition (NHK data), = 0 use data from Tokyo Keizai
% in the "Figure" folder
fs = 20; % common font size for many figures
ldfs = 12; % legend font size for vaccine path
ft = '%.1f';
if iPC == 1
    fn = 'Yu Gothic'; % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
linecolor = {'blue', 'black', 'red'};
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


% [V1_medical_ori, V1_medical, V2_medical_ori, V2_medical,...
%     V1_elderly_ori, V1_elderly, V2_elderly_ori, V2_elderly, ...
%     vs_MT,vs_ET,vs_M1,vs_E1,vs_M2,vs_E2] ...
%     = ImportVaccineData(home,iPC,Data,pref,dateEN,ps,vaccine_disp_switch);

%--- Constructing the reference level of output ---%
[potentialGDP, referenceGDP, alpha] = construct_GDP(GDP,TdataGDP);

%--- Regress mobility on alpha to estimate the elasticity h ---%
[Malt,h_all,h_all_se,h_ori,h_se] = estimate_h(M,alpha,TdataGDP,RetroH,hconstant);

%--- Plot mobility data ---%
figname = string(['Mobility_GDP_' char(pref)]);
f = figure('Name',figname);
plot_mobility(Malt,alpha,Tdata,TdataGDP,MonthWeekJP,xtick1,fs,16)
if figure_save == 1
    %saveas(figure(2),[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
    saveas(f,[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
end

%--- Import ICU data (Option 2) ---%
ICU = zeros(Tdata+1,1);
BED = zeros(Tdata,1);
if ICU_nation == 1
    ICU(2:Tdata+1,1) = Data(:,22);
    BED(1:Tdata,1) = Data(:,23);
else
    ICU(2:Tdata+1,1) = Data(:,21);
end
%--- Plot ICU data ---%
figname = string(['ICU_transition_' char(pref)]);
f = figure('Name',figname);
plot(ICU, 'LineWidth', 1.5)
title('Transition of  ICU')
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


%--- Compute the history of time-varying parameters ---%
[delta,beta_tilde,ERN,beta,ICU_inflow,...
    gammaT,delta_average,ICU_inflow_avg,delta_sample]...
    = Time_Series_Average(S,I,D,ICU,dD,N,Tdata,SimPeriod,...
    RetroPeriod,POP0,gamma,hconstant,h_all,alpha,k,...
    gamma_ICU,ICU_adjustment,RetroPeriodDelta);

elderly_total = ps*elderly_jp;
medical_total = ps*medical_jp;
ordinary_total = ps*ordinary_jp;
%medical = medical_total*accept_share;
medical = medical_total;
elderly = elderly_total*accept_share;
ordinary = ordinary_total*accept_share;
%         elderly = elderly - (sum(V1_elderly));
%--- Eliminate the effects of vaccination from delta ---%
delta_ss = delta_average*(0.1063/1.53); %Share of the death rate among youth to the death rate of all populaiton
VD_elderly = D1*V1_elderly + (D2-D1)*V2_elderly;
VD_ordinary = (D1*V1_medical + (D2-D1)*V2_medical) + (D1*V1_others + (D2-D1)*V2_others);
share = ((1 - sum(VD_elderly(1:end-2))/elderly_total) * (delta_average - delta_ss) ...
    + (1-sum(VD_ordinary(1:end-2))/(ordinary_total))*delta_ss)/delta_average;
delta_average = delta_average/share;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Main analysis starts here %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ======================================================= %
% ==================== CHANGE HERE ====================== %
% ======================================================= %

% TH = [17,21,25,29,33,37,40];
% TH_index = [17,25,33];

TH = [0.1, 0.2,0.3];
TH = (1+TH).*(1-TH .* 0.05) - 1; % TH = [1.1*0.995 - 1,1.2*0.99 - 1,1.3*0.985-1];

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
    
    state = 0;
    %alpha_scale = 0.94;
    alpha_scale = 1.0;
    alpha_on = (0.5 * alpha_May + 0.5 * alpha_Jan); %alpha_Jan; %alpha_scale * (0.5 * alpha_May + 0.5 * alpha_Jan);
    h_scale = 1;
    beta_shock_after_emergency = 0;
    rho_after_emergency = 0.5;
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
h(2) = h_scale*h_ori(2);


% ======================================================= %
% ======================================================= %
% for iAlpha = 1:length(alpha_on_vector_sim)
%     alpha_on = alpha_on_vector_sim(iAlpha)
nTH = length(TH);
DMat = nan(1,length(TH));
AlphaMat = nan(1,length(TH));
SimData = nan(SimPeriod+1,5,length(TH));
AlphaPath = nan(SimPeriod,length(TH));
NPath = nan(SimPeriod,length(TH));
SimERN = nan(SimPeriod,length(TH));
THonPath = nan(SimPeriod,length(TH));
BackDataN = zeros(SimPeriod+8,length(TH_index));
BackDataAlpha = zeros(SimPeriod+8,length(TH_index));
BackDataERN = zeros(SimPeriod+8,length(TH_index));
BackDataDA = zeros(length(TH),3);
betaPath = zeros(SimPeriod,nTH);
betaAlphaPath = zeros(SimPeriod,nTH);
betaTildePath = zeros(SimPeriod,nTH);
deltaPath = zeros(SimPeriod,nTH);
deltaICUPath = zeros(SimPeriod,nTH);
varSharePath = zeros(SimPeriod,nTH);
varDeltaSharePath = zeros(SimPeriod,nTH);

for iTH = 1:length(TH)
    
    %%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%
    var_infection2 = TH(iTH);
    
    %%%%%%%%%%%%%%%%% Figure Titles %%%%%%%%%%%%%%%%%%%%%%%%%
    figname_main = 'Delta_MainResults_';
    figname_beta = 'BetaPath';
    figname_beta_loop = ['BetaPath' '_' sprintf(ft, th_off1 / 7)];
    figname_beta_tilde_loop = ['BetaTildePath_', sprintf(ft, th_off1 / 7)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        elderly,elderly_total,medical,ordinary,delta_average,lag,...
        medical_duration,paces_ori,paces2,paces3,sw_vacpath,gradual_paces,gradual_paces2,...
        E1,E2,D1,D2,ps,POP0,SimPeriod,Tdata);
    
    [var_share,var_prev,var_initial] = var_share_prev(Data(:,20),SimPeriod,var_ss,var_growth);
%     var_prev = Data(:,20);
%     var_prev(isnan(var_prev))=0;
%     var_prev(71:end) = 1;
%     var_share = ones(SimPeriod,1);
%     [var_share_delta,var_prev_delta,var_initial_delta] = var_share_prev(Data(:,24),SimPeriod,var_ss,var_growth2);
    [deltaT,delta_bar] = construct_delta_variant(RetroPeriodDelta,...
        retro_lb,retro_ub,deltaT,delta,delta_average,var_prev,var_share,var_infection_delta,I);
    [delta_ICU,delta_ICU_bar] = construct_delta_variant(RetroPeriodDelta,...
        retro_lb,retro_ub,delta_ICU,delta,delta_average,var_prev,var_share,var_infection_delta,I);
    deltaPath(:,iTH) = deltaT;
    deltaICUPath(:,iTH) = delta_ICU;
    varSharePath(:,iTH) = var_share;
    
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
    beta_r = 0;
    for retrop = retro_lb:retro_ub
        beta_r = beta_r + mean(beta(end-retrop+1:end));
    end
    beta_avg = beta_r/(retro_ub-retro_lb+1);
    
    [betaT,betaT_woAR1,beta_bar] = construct_beta(SimPeriod,retro_lb,retro_ub,...
        beta,beta_avg,var_prev,var_share,var_infection);
    
    %%%%  デルタ株　　　%%%%
    %%%%%%%%%%%%%%%%%%%%%​
    var_intercept2 = log(var_initial2/(var_ss-var_initial2));
    var_share2 = exp((1:SimPeriod)'*var_growth2+var_intercept2).*var_ss./(1+exp((1:SimPeriod)'*var_growth2+var_intercept2));
    %     var_power = 1+var_share2*var_growth2;
    %     var_boost = ones(SimPeriod,1);
    %     var_boost(var_start:end) = var_power(1:SimPeriod-var_start+1);
    var_shareInd = zeros(SimPeriod,1);
    var_shareInd(var_start,1) = var_initial2;
    var_shareInd(var_start+1:end,1) = var_share2(1:SimPeriod-var_start);
    varDeltaSharePath(:,iTH) = var_shareInd;
    betaT_Eng = betaT;
    betaT_Eng = beta_AR1(betaT_temp_ini, beta_rho, betaT_Eng, start_beta);
    betaT = betaT.*(1+var_shareInd*var_infection2);
    
    betaT = beta_AR1(betaT_temp_ini, beta_rho, betaT, start_beta);
    
    
    InitialValues = [S(end),I(end),R(end),D(end),ICU(end)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [DMat(iTH),AlphaMat(iTH),AlphaPath(:,iTH),SimData(:,:,iTH),NPath(:,iTH),SimERN(:,iTH),THonPath(:,iTH),SimICU(:,iTH),betaShock] ...
        = Covid_projection_ICU5...
        (InitialValues,alpha_on,alpha_off,...
        th_on1,th_on2,th_off1,th_off2,th_off3,...
        betaT,gammaT,deltaT,deltaT(1),delta_ICU,V,h,k,POP0,hconstant,...
        DRi,state,ICU_inflow_avg,gamma_ICU,...
        beta_shock_after_emergency,rho_after_emergency,alpha_jump);
    % = Covid_projection_ICU4...
    % (InitialValues,alpha_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,betaT,gammaT,deltaT,deltaT(1),V,h,k,POP0,hconstant,DRi,state,ICU_inflow_avg,gamma_ICU,beta_shock_after_emergency,rho_after_emergency);
    betaPath(:,iTH) = betaShock;
    betaShock_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaShock;
    betaTildePath(:,iTH) = betaShock_tilde;
    
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

%minAlpha = min(minAlphaMat); %minimum alpha when variants have no effects.
minAlpha = alpha_off; % 経済損失0 = 2020年2月のGDP level

AlphaM = AlphaMat(~isnan(AlphaMat));
AlphaM = (AlphaM - minAlpha)*prefGDP*10000;
DM = DMat(~isnan(DMat));
BackDataDA(1:length(TH),:) = [round(AlphaM'),round(DM'),1+round(TH',1)];
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
        set(gcf,'Position',[100,100,1200,800])
    end
    subplot(2, 2, 1)
    lineName = cell(nTH,1);
    for i = 1:nTH
        if l == 1
            lineName{i} = ['Relative Infection Rate of Delta Variants = ' sprintf(ft,1 + TH_index(i))];
        elseif l == 2
            lineName{i} = ['デルタ株の相対感染力:' sprintf(ft,1 + TH_index(i)) '倍'];
        end
    end
    lineStyle = {'-','-','-'};
    lgfs = 12;
    eTitle = 'Projected Path of New Cases';
    jTitle = '新規感染者数の推移';
    show_other = 0; % Do not show lines other than TH index
    yft = '%.0f';
    column = 1;
    plot_function(TH,TH_index,N/7,NPath/7,...
        MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineName,...
        eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'Southwest')
    xticks(find(WeekNumber==1))
    xtickangle(45)
    xlim([Tdata-7 Tdata+25])
                                              
    %--- Number of cumulative deaths ---%
    subplot(2, 2, 2)
    plot_Tradeoff2(AlphaM, DM, waves, TH, TH_index, l, linecolor, fs, fn)
    
    %--- Number of people who are in ICU ---%
    subplot(2,2,3)
    if ICU_nation == 1
        BackDataICU = ...
            plot_ICU_GOV(TH,TH_index,ICU,SimICU,BED(Tdata),...
            MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,...
            linecolor,fs,fn,ft,l);
    else
        BackDataICU = plot_ICU(TH,TH_index,ICU,SimICU,ICU_limit,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,fs,fn,ft,l,th_off1);
    end
    
    subplot(2,2,4)
    MonthWeek = [MonthWeekEN,MonthWeekJP];
    plot_Alpha(alpha,AlphaPath,TH,TH_index,MonthWeek(:,l),WeekNumber,MonthNumber,Tdata,linecolor,ft,fs,fn,l)
    title('GDP','FontSize',fs,'FontWeight','normal')

    %xlim([Tdata-7 Tdata+33])
    %         subplot(2,2,4)
    %         plot_ICU_N(TH,TH_index,N,NPath,ICU,SimICU,ICU_limit,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,20,fn,ft,l,th_off1/7,th_off2/7,th_off3/7)
    
    %         %--- Plot ICU inflow ---%
    %         figname = 'ICU_inflow';
    %         f = figure('Name',figname);
    %         set(gcf,'Position',[100,100,1200,500])
    %         subplot(1,2,1)
    %         plot(ICU_inflow, 'LineWidth', 1.5)
    %         title('ICU inflow')
    %         ytickformat('%,6.0f')
    %         xticks(find(WeekNumber==1))
    %         xticklabels(MonthWeekJP(WeekNumber==1))
    %         lgd.NumColumns = 2;
    %         xtickangle(45)
    %         subplot(1,2,2)
    %         plot(ICU_inflow.*delta, 'LineWidth', 1.5)
    %         title('ICU inflow * \delta')
    %         xticks(find(WeekNumber==1))
    %         xticklabels(MonthWeekJP(WeekNumber==1))
    %         lgd.NumColumns = 2;
    %         xtickangle(45)
    %         xlim([20 Tdata])
    
    if figure_save == 1
        saveas(f,[home 'Figures/' char(pref) '/' char(figname_main) char(lng) '.png']);
    end
end %End of language loop = figure loop

for i = 1:length(TH_index)
    BackDataN(:,i) = [N(end-7:end);NPath(:,i)];
    BackDataAlpha(:,i) = [alpha(end-7:end);AlphaPath(:,i)];
    BackDataERN(:,i) = [ERN(end-7:end);SimERN(:,i)];
end

if data_save == 1
    titleN = strings(1,1+length(TH_index)*3);
    titleN(1) = "週";
    for ti = 1:length(TH_index)
        titleN(1,1+ti) = string(['新規感染者数（',sprintf('%.1f',1+TH_index(ti)),'）']);
        titleN(1,1+length(TH_index)+ti) = string(['経済活動（',sprintf('%.1f',1+TH_index(ti)),'）']);
        titleN(1,1+length(TH_index)*2+ti) = string(['実効再生産数（',sprintf('%.1f',1+TH_index(ti)),'）']);
        titleN(1,1+length(TH_index)*3+ti) = string(['重症者数（',sprintf('%.1f',1+TH_index(ti)),'）']);
    end
    TN = table([titleN;MonthWeekJP(Tdata-7:end-1),round(BackDataN(:,1:length(TH_index))/7),round(100*(1-BackDataAlpha(:,1:length(TH_index))),1),round(BackDataERN(:,1:length(TH_index)),2),round(BackDataICU(2:end,1:length(TH_index)))]);
    titleAD = ["経済損失（億円）","死亡者数","ケース"];
    TAD = table([titleAD;BackDataDA(1:length(TH),:)]);
    writetable(TN,[home 'Figures/' char(pref) '/BackData_' char(figname_main)  char(pref) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
    writetable(TAD,[home 'Figures/' char(pref) '/BackData_' char(figname_main)  char(pref) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
end

%figname = 140 + pindex; %Plotting New Cases + Trade Off + Alpha Path
%figure(figname)
% f2=figure('Name',[figname_main char("Alpha_") lng]);
% if iPC == 1
%     f.WindowState = 'maximized';
% else
%     set(gcf,'Position',[100,100,1200,800])
% end
% subplot(2, 2, 1)
% plot_SimN_DR_km(TH, TH_index, N, NPath, alpha, AlphaPath, ERN, SimERN, THonPath, MonthWeekEN, MonthWeekJP, WeekNumber, Tdata, linecolor, 20, fn, ft, l,20);
% % plot_SimN(TH, TH_index, N, NPath, alpha, AlphaPath, ERN, SimERN, THonPath, MonthWeekEN, MonthWeekJP, WeekNumber, Tdata, linecolor, 20, fn, ft, l);
% subplot(2, 2, 2)
% plot_Tradeoff2(AlphaM,DM,waves,TH,TH_index,l,linecolor,fs,fn)
% % legend('3ヶ月','5ヶ月','7ヶ月','FontSize',10,'FontName',fn,'Location','northwest');
% subplot(2,2,3)
% plot_ICU_GOV(TH,TH_index,ICU,SimICU,BED(Tdata),...
%     MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,...
%     linecolor,fs,fn,ft,l);
% subplot(2,2,4)
% plot_Alpha(alpha,AlphaPath,TH,TH_index,MonthWeekJP,WeekNumber,MonthNumber,Tdata,linecolor,ft,fs,fn,2)

%     end
% end %end of prefecture loop
%
%

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
    lineName{i} = ['I, 宣言解除 = ' sprintf(ft, TH(i))];
    lineName2{i} = ['N, 宣言解除 = ' sprintf(ft, TH(i))];
end
eTitle = 'Transitions of N and I';
jTitle = 'NとIの推移';
figure('Name',char(figname_NI));
set(gcf,'Position',[100,100,1200,800])
plot_function2(TH,TH_index,I(2:end)/7,SimData(2:end,2,:)/7,N/7,NPath/7,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineName,...
    lineName2,eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
xlim([Tdata-31 Tdata+41])

% Transitions of S and R
figname_SR = 'SR_transitions';
lineName = cell(nTH,1);
lineName2 = cell(nTH,1);
for i = 1:nTH
    lineName{i} = ['S, 宣言解除 = ' sprintf(ft, TH(i))];
    lineName2{i} = ['R, 宣言解除 = ' sprintf(ft, TH(i))];
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
    lineName{i} = ['\beta, 宣言解除 = ' sprintf(ft, TH(i))];
end
eTitle = 'Transitions of \beta';
jTitle = '\betaの推移';
yft = '%.3f';
plot_function(TH,TH_index,beta,betaPath,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineName,...
    eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
plot([NaN(Tdata-1,1);ones(SimPeriod+1,1)*beta_avg],':k','LineWidth',1.5,'DisplayName','過去17週間平均')
plot([NaN(Tdata-1,1);ones(SimPeriod+1,1)*beta_bar],'--k','LineWidth',1.5,'DisplayName','過去17週間平均, 変異株影響除去')
plot([NaN(Tdata,1);betaT],'--g','LineWidth',1.5,'DisplayName','\beta w/o AfterShock')
plot([NaN(Tdata,1);betaT_Eng],'--c','LineWidth',1.5,'DisplayName','\beta w/o \delta Variant Effects')
xline(Tdata-RetroPeriod-1,'HandleVisibility','off')
xlim([Tdata-31 Tdata+41])

% Transitions of Beta Tilde
figname_beta_tilde = 'Beta Tilde transition';
figure('Name',char(figname_beta_tilde));
set(gcf,'Position',[100,100,1200,800])
lineName = cell(nTH,1);
for i = 1:nTH
    lineName{i} = ['\beta tilde, 宣言解除 = ' sprintf(ft, TH(i))];
end
eTitle = 'Transitions of \beta tilde';
jTitle = '\beta tildeの推移';
yft = '%.3f';
plot_function(TH,TH_index,beta_tilde,betaTildePath,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineName,...
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
lineName = cell(nTH,1);
lineName2 = cell(nTH,1);
for i = 1:nTH
    lineName{i} = ['\alpha 変異株割合, 宣言解除 = ' sprintf(ft, TH(i))];
    lineName2{i} = ['\delta 変異株割合, 宣言解除 = ' sprintf(ft, TH(i))];
end
eTitle = 'Transitions of Variant Rate';
jTitle = '変異株割合の推移';
plot_function2(TH,TH_index,var_prev,varSharePath,zeros(Tdata,1),varDeltaSharePath,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineName,...
    lineName2,eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
xlim([Tdata-31 Tdata+41])
if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref)  '/Variant_Share.png']);
end

% Transitions of Death Rate
figname_delta = 'Death Rate transition';
figure('Name',char(figname_delta));
set(gcf,'Position',[100,100,1200,800])
lineName = cell(nTH,1);
lineName2 = cell(nTH,1);
column = 2;
for i = 1:nTH
    lineName{i} = ['死亡率, 宣言解除 = ' sprintf(ft, TH(i))];
    lineName2{i} = ['死亡率, 宣言解除 = ' sprintf(ft, TH(i))];
end
eTitle = 'Transitions of Death Rate';
jTitle = '死亡率の推移';
plot_function(TH, TH_index, delta, deltaPath, ...
    MonthWeekEN, MonthWeekJP, WeekNumber, Tdata, xline_ind, linecolor, lineName, ...
    eTitle, jTitle, fs, lgfs, fn, ft, yft, column, l, show_other, 'SouthOutside')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * delta_average], ':k', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均')
plot([NaN(Tdata - 1, 1); ones(SimPeriod + 1, 1) * delta_bar], '--k', 'LineWidth', 1.5, 'DisplayName', '過去17週間平均, 変異株影響除去')
xline(Tdata - RetroPeriodDelta - 1, 'HandleVisibility', 'off')
xlim([Tdata - 31 Tdata + 41])
%
% Transitions of ICU Rate
figname_deltaICU = 'ICU Rate transition';
figure('Name',char(figname_deltaICU));
set(gcf,'Position',[100,100,1200,800])
lineName = cell(nTH,1);
lineName2 = cell(nTH,1);
for i = 1:nTH
    lineName{i} = ['重症化率, ワクチン接種本数/週 = ' sprintf(ft,TH(i))];
    lineName2{i} = ['重症化率, ワクチン接種本数/週 = ' sprintf(ft,TH(i))];
end
eTitle = 'Transitions of ICU Rate';
jTitle = '重症化率の推移';
plot_function(TH,TH_index,delta,deltaICUPath,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineName,...
    eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,'SouthOutside')
xline(Tdata-RetroPeriodDelta-1,'HandleVisibility','off')
xlim([Tdata-31 Tdata+41])

