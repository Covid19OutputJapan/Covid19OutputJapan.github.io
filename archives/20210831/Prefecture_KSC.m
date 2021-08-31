% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
iPC=0; % 0 for Mac, 1 for Windows
if iPC==1
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
%        home =     'C:\Users\tak10\Dropbox\fujii_nakata\Policy_request\Cabinet_2021MAY07_newV\';
else
    %home ='/Users/machikohei/Dropbox/fujii_nakata/Website/Codes/Cabinet_2021APR28_DR/';
    home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/';
    %home = '/Users/ymaeda/Documents/ym_doc/Tokyo_Univrsity_MA/Reserach Assistant/fujii_nakata/Codes/';
    %home = '/Users/shotaro/Dropbox/fujii_nakata/Website/Codes/Cabinet_2021APR28_DR/';
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 0;      % save back data
vaccine_figure_loop = 0; % =0 appear only once; =1 appear every loop; 
beta_figure_loop = 0; % =0 appear only once; =1 appear every loop; 
vaccine_disp_switch = 1; % =0 not display the summry of # of the vaccinated ; =1 display 
% in the "Figure" folder
fs = 16;            % common font size for many figures
ldfs = 8;            % legend font size for vaccine path
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
linecolor = {'red','blue','black','green'};
language = {'EN','JP'};
%======================================================================%

%================== Model Fixed Parameter Values ============================%
SimPeriod = 52;        % simulation period in weeks
gamma = 7/12;          % recovery rate from Covid % Should change this to 7/12 (4/25 Kohei Machi)
k = 2;                 % exponent of (1-h*alpha)
hconstant = 1;         % 0 = without intercept, 1 = with intercept for h regression
medical_start_date = datetime(2021,3,18);
elderly_start_date = datetime(2021,5,13);
RetroPeriod = 17;      % retroactive periods used to estimate gamma and delta
tt = 12; % Showing previous t periods for the plot
% Parameters for beta
beta_rho = 0.85;
retro_ub = 17; % Control the moving average of beta (beta_avg = sum_{t = lb}^{ub} (1/(ub-lb + 1) sum_{x=1}^t (1/t) beta_t)
retro_lb = 17;
% Parameters for mobility estimation
retroH_switch = 1; %If retroH_switch == 1, retroH = TdataGDP - 4, else = retroH
RetroH = 15;
% Parameters for variants
var_infection = 0.5;
var_ss = 1.0;         % steady-state share of variant
var_growth = 0.47;
% Parameters for Vaccine Path
paces_ori = 3500000;
gradual_paces = 6;
sw_vacpath = 0;
% Population
POP_jp = 125710000;
medical_jp = 4700000;
elderly_jp = 36000000;
ordinary_jp = (POP_jp-elderly_jp-medical_jp);
accept_share = 0.8;
% vaccine pace
PF = 1; % 0 for AZ, 1 for PF
if PF == 0  % AZ
    E1 = 0.615;
    E2 = 0.64;
    D1 = 0.8;
    D2 = 0.85;
else   % PF
    E1 = 0.625;
    E2 = 0.895;
    D1 = 0.8;
    D2 = 0.94;
end

%================== Parameter Values (Prefecture Specific) ============================%
PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka','Aichi','Fukuoka','Hyogo'};
GDPVector = [106,36,23,21,40,41,20,20]; % 兆円, one trillion yen (chou-yen)
% 緊急事態宣言の発令基準
th_on_vector = [750,500,400,350,1000,350,350]; % present.
th_on_vector2 = [2000,1000,800,700,2000,350,350]; % 高齢者がうち終わったあとの基準

% 緊急事態宣言の解除基準
th_off_vector = [550,100,110,80,825,50,60]; %1回目の緊急事態宣言の解除基準
th_off_vector2 = [400,100,80,60,500,50,60]; %2回目の緊急事態宣言の解除基準
th_off_vector3 = [400,100,80,60,500,50,60]; %3回目の緊急事態宣言の解除基準
% 経済回復速度
DRi_vector = [6,6,6,6,6,6,6,6];
% 緊急事態宣言の強さの基準: Simulation開始時のERNをいくつにするか
%ERN_on_vector =  [0.55,0.5,0.5,0.5,0.65]; % 緊急事態宣言下のERN
alpha_on_vector = [0.2,0.1,0.1,0.1,0.2];
% Size of AR(1) shock for beta process
betaT_temp_ini_vec = [0.2,0,0,0,0.15,0,0,0];
% Initial Share of Variants (Need to change these values every week)
%var_initial_vector = [0.2531, 0.2909, 0.1335, 0.0909, 0.7721, 0, 0, 0,0]; % 4/26
% var_initial_vector = [0.556244, 0.450617, 0.514451, 0.422360, 0.824487, 0, 0, 0,0]; % 5/10 (4/19-4/25のデータ)
var_initial_vector = [0.631038, 0.596774, 0.578947, 0.56,     0.856401, 0, 0, 0, 0];% 5/17 (4/26-5/2のデータ)
% var_ss_vector = [1, 1, 1, 1, 0.85, 0, 0, 0,0]; % 5/10
var_ss_vector = [1, 1, 1, 1, 0.90, 0, 0, 0,0]; % 5/17
share_index = 2; % = 2 for Monday, = 1 after Wednesday
%=========== use this code to extrapolate var_initial_vector ==================%
var_initial_vector = ini2now_infection_rate(var_initial_vector,var_growth,var_ss_vector,SimPeriod,share_index); 
%0.7624    0.6774    0.7306    0.6518    0.8398         0         0         0         0
%=====================================================%


% ------------------------- Experiment Vectors ----------------------- %
% different threshold for declaring the state of emergency
THON_vector = {[900:300:1500,2000:1000:4000],400:100:800,200:100:600,300:50:500,300:100:1000,200:50:400,200:50:400,300:50:500};
THON_index_vector = {[1500,3000],[600,800],[300,600],[350,500],[400,600],[400,600],[300,600],[350,500]};

% Duration of Recovery
DR_vector = {4:2:16,0:2:16,0:2:16,0:2:16,6};
DR_index = [6,12]; %色付けして強調するものを決める

% Different ERN during the state of emergency
ERNON_vector = {0.3:0.02:0.7,[],[],[],0.3:0.02:0.7};

% Variant Infection
var_infection_vec = 0:0.1:0.6; %[0.3, 0.5]
var_infection_index = [0.3, 0.5, 0.7];

% Vaccine Paces
paces_vector = [3500000,7000000];

% Different threshold for lifting the state of emergency
TL_vector = {100:50:500,70:10:100,80:10:130,50:10:100,100:50:500,20:10:100,10:10:100,50:10:100}; % 解除基準分析をコントロールしている Cell array
TL_index_vector = {[250,500,100],[100,80],[130,110],[100,80],[250,500,100],[50,20],[60,10]}; %色付けして強調するものを決める

for pindex = 2:2 %:length(PrefVector) %change this parameter for prefecture
    %====================== Model parameter values ======================%
    pref = PrefVector{pindex};        % prefecture to be analyzed
    prefGDP = GDPVector(pindex);
    %====================================================================%
    betaT_temp_ini = betaT_temp_ini_vec(pindex);
    var_initial = var_initial_vector(pindex);
    var_ss = var_ss_vector(pindex);
    
    
    %--- Import data ---%
    % Covid data are recorded at weekly frequencies (Mon-Sun)
    % The first week start on January 20 (Mon), 2020
    import_prefecture
    
    %--- Construct weekly vaccine data ---%　
    [V1_medical_ori, V1_medical, V2_medical_ori, V2_medical,...
        V1_elderly_ori, V1_elderly, V2_elderly_ori, V2_elderly, ...
        vs_MT,vs_ET,vs_M1,vs_E1,vs_M2,vs_E2] ...
        = ImportVaccineData(home,iPC,pref,dateEN,ps,vaccine_disp_switch);
    %------------------------------------------------------------------------%
     
     %--- Constructing the reference level of output ---%
    [potentialGDP, referenceGDP, alpha] = construct_GDP(GDP,TdataGDP);
    
    %--- Regress mobility on alpha to estimate the elasticity h ---%　関数化？
    [Malt,h_all,h_all_se,h,h_se] = estimate_h(M,alpha,TdataGDP,RetroH,hconstant);

    %--- Plot mobility data ---%
    figname = 'Mobility_GDP';
    f = figure('Name',figname);
    plot_mobility(Malt,alpha,Tdata,TdataGDP,MonthWeekJP,xtick1,fs,16)
    if figure_save == 1
        %saveas(figure(2),[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
        saveas(f,[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Main analysis starts here %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Projection parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DRi = DRi_vector(pindex);
    th_on1 = th_on_vector(pindex)*7;         % threshold to place the state of emergency
    th_on2 = th_on_vector2(pindex)*7;         % threshold to place the state of emergency
    th_off1 = th_off_vector(pindex)*7;
    th_off2 = th_off_vector2(pindex)*7;
    th_off3 = th_off_vector3(pindex)*7;
    %ERN_on = ERN_on_vector(pindex);
    alpha_on = alpha_on_vector(pindex); %(((ERN_on*(POP0/S(end))*((gammaT(1)+delta_average)/beta_avg)).^(1/k))-1)*(h(1)/h(2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% CHANGE HERE %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if pindex == 1 || pindex == 5
        state = 1;
    else 
        state = 0;
    end 
    alpha_May = mean(alpha((dateEN >= datetime(2020,5,07)) & (datetime(2020,5,28)>= dateEN ))); 
    alpha_Jan = mean(alpha((dateEN >= datetime(2021,1,07)) & (datetime(2021,1,28)>= dateEN ))); 
    %alpha_on = 0.5*alpha_May + 0.5*alpha_Jan;
    alpha_on = alpha_May;
    %alpha_on_vector = [alpha_on1, alpha_on2];
    TH = var_infection_vec + 1;
    TH_index= [1.3,1.5];
    if max(TH)<3
        ft = '%.2f';
    else
        ft = '%.0f';
    end
    % Change both paces and infection rate at the same time
    %     TH = [3600000,7000000,3600000,7000000];  %TH = paces_vector;
    %     TH_index = [3600000,7000000,3600000,7000000]; %TH_index = paces_vector;
    %     TH2 = [0.3, 0.3 ,0.5, 0.5];
    %     TH2_index = [0.3, 0.5];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    DMat = nan(1,length(TH));
    AlphaMat = nan(1,length(TH));
    SimData = nan(SimPeriod+1,4,length(TH));
    AlphaPath = nan(SimPeriod,length(TH));
    NPath = nan(SimPeriod,length(TH));
    SimERN = nan(SimPeriod,length(TH));
    BackDataN = zeros(SimPeriod+8,length(TH_index));
    BackDataAlpha = zeros(SimPeriod+8,length(TH_index));
    BackDataERN = zeros(SimPeriod+8,length(TH_index));
    BackDataDA = zeros(length(TH),3);
%     for iAlpha = 1:length(alpha_on_vector)
%         alpha_on = alpha_on_vector(iAlpha)
    for iTH = 1:length(TH)
        %%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%
        % paces_ori = TH(iTH);
        % var_infection = TH2(iTH);
%         th_off1 = TH(iTH)*7
%         th_off2 = TH(iTH)*7;
%         th_off3 = TH(iTH)*7;
        var_infection = TH(iTH) - 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %--- Compute the history of S, I, R, D in the data period ---%
        S = zeros(Tdata+1,1);
        I = zeros(Tdata+1,1);
        R = zeros(Tdata+1,1);
        D = zeros(Tdata+1,1);
        S(1)=POP0;
        for i = 1:Tdata
            S(i+1)=S(i)-N(i)-E1*(V1_elderly(i)+V1_medical(i))-(E2-E1)*(V2_elderly(i)+V2_medical(i));
            I(i+1)=I(i)+N(i)-gamma*I(i)-dD(i);
            R(i+1)=R(i)+gamma*I(i)+E1*(V1_elderly(i)+V1_medical(i))+(E2-E1)*(V2_elderly(i)+V2_medical(i));
            D(i+1)=D(i)+dD(i);
            if i > TdataGDP
                GDP(i) = referenceGDP(i)*(1-alpha(i));
            end
        end
        
        %--- Compute the history of time-varying parameters ---%
        delta = (D(2:Tdata+1)-D(1:Tdata))./I(1:Tdata);                              % death rate
        beta_tilde = (POP0.*N(1:Tdata))./((S(1:Tdata).*I(1:Tdata)));   % overall infection rate
        ERN = (S(1:end-1)/POP0).*beta_tilde./(gamma+delta);                                        % effective reproduction number
        if hconstant == 0
            beta = beta_tilde./(1+h_all*alpha).^k;                                      % raw infection rate
        elseif hconstant == 1
            %     beta = beta_tilde./(h(1)+h(2)*alpha).^k;
            beta = beta_tilde./(1+(h_all(2)/h_all(1))*alpha).^k;
        end
        
        %--- Construct time series of parameters ---%
        gammaT = gamma*ones(SimPeriod,1);
        delta_sample = delta(end-RetroPeriod+1:end);
        delta_average = sum(delta_sample.*(I(end-RetroPeriod+1:end)/sum(I(end-RetroPeriod+1:end))));
        
        %--- Construct vaccine dstribution ---%
        paces = ps*paces_ori; %3600000;
        vacpath = zeros(SimPeriod,1);
        vacpath(1+sw_vacpath:gradual_paces) = (paces/(gradual_paces-sw_vacpath)):(paces/(gradual_paces-sw_vacpath)):paces;
        vacpath(gradual_paces+1:end) = paces*ones(SimPeriod-gradual_paces,1);
        elderly_total = ps*elderly_jp;
        medical_total = ps*medical_jp;
        ordinary_total = ps*ordinary_jp;
        medical = medical_total*accept_share;
        elderly = elderly_total*accept_share;
        ordinary = ordinary_total*accept_share;
        
        elderly = elderly - (sum(V1_elderly));
        % Original
%         [V,deltaT,VT] = vaccine_distribution_medical(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2,ps,POP0,3);
        [V,deltaT,VT] = vaccine_distribution_medical(vacpath,medical,V1_medical,V2_medical,elderly,V1_elderly,V2_elderly, ordinary,elderly_total,delta_average,E1,E2,D1,D2,ps,POP0,3,10);
        %[V,deltaT,VT] = vaccine_distribution_simple(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2);
        delta_ss = delta_average*(0.09/1.28);
        delta_ss = (delta_average - delta_ss)*(1-accept_share*D2)+delta_ss; %death rate after full vaccination of elderly = 2021/09/30
        
        
        
        
        
        %% figure for vaccine path 
        if vaccine_figure_loop == 0
            if iTH == 1
                 plot_vaccinepath(200,VT,V1_medical,V2_medical,V1_elderly,V2_elderly,SimPeriod,ps,MonthWeekJP,WeekNumber,Tdata,fs,ldfs,fn);
                
                plot_deltapath(201,delta,deltaT,delta_average,MonthWeekJP,WeekNumber,Tdata,fs,fn);
                
            end
        elseif vaccine_figure_loop == 1
             plot_vaccinepath(200,VT,V1_medical,V2_medical,V1_elderly,V2_elderly,SimPeriod,ps,MonthWeekJP,WeekNumber,Tdata,fs,ldfs,fn);
            subplot(1,2,1)
            title(string(['新規ワクチン接種本数（週ごと）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
            subplot(1,2,2)
            title(string(['累計ワクチン接種本数（週ごと）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);

            plot_deltapath(240+iTH,delta,deltaT,delta_average,MonthWeekJP,WeekNumber,Tdata,fs,fn);
            title(string(['致死率（現在のレベルで標準化）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        beta_r = 0;
        for retrop = retro_lb:retro_ub
            beta_r = beta_r + mean(beta(end-retrop+1:end));
        end
        beta_avg = beta_r/(retro_ub-retro_lb+1);
        [betaT,betaT_woAR1,var_share,var_share_prev] = construct_beta(home,SimPeriod,Tdata,dateEN,pref,...
            retro_lb,retro_ub,beta,beta_avg,var_initial,var_ss,var_infection,var_growth,betaT_temp_ini,beta_rho);
        
        alpha_off = mean(alpha((dateEN >= datetime(2020,10,1)) & (datetime(2020,11,26)>= dateEN ))); % output loss without the state of emergency
        InitialValues = [S(end),I(end),R(end),D(end)];
        ERNCheck = (S(end)/POP0).*(((1+(h(2)/h(1))*alpha_off).^k).*beta_avg)./(gammaT(1)+delta_average);
        
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if state == 0 %pindex < 5
            [DMat(iTH),AlphaMat(iTH),AlphaPath(:,iTH),SimData(:,:,iTH),NPath(:,iTH),SimERN(:,iTH),THonPath(:,iTH)] ...
                = Covid_projection_control_gradual_off_threshold_delta2...
                (InitialValues,alpha_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,betaT,gammaT,deltaT,delta_ss,delta_average,V,h,k,POP0,hconstant,DRi,alpha(end));
%               Construct no-variants-coutnerfactual alpha matrix
%             [~,minAlphaMat(iTH)] ...
%                 = Covid_projection_control_gradual_off_threshold_delta2...
%                 (InitialValues,alpha_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,beta_avg*ones(SimPeriod,1),gammaT,deltaT,delta_ss,V,h,k,POP0,hconstant,DRi,alpha(end));
        elseif state == 1
            [DMat(iTH),AlphaMat(iTH),AlphaPath(:,iTH),SimData(:,:,iTH),NPath(:,iTH),SimERN(:,iTH),THonPath(:,iTH)] ...
                = Covid_projection_control_gradual_threshold_delta2...
                (InitialValues,alpha_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,betaT,gammaT,deltaT,delta_ss,delta_average,V,h,k,POP0,hconstant,DRi);
%             [~,minAlphaMat(iTH)] ...
%                 = Covid_projection_control_gradual_threshold_delta2...
%                 (InitialValues,alpha_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,beta_avg*ones(SimPeriod,1),gammaT,deltaT,delta_ss,V,h,k,POP0,hconstant,DRi);
        end
        
        % Plot betaT 
        if beta_figure_loop == 0
            if iTH == 1
                %figure(400) % Figure for BetaT with Variant Share
                f=figure('Name','BetaPath');
                set(gcf,'Position',[100,100,1200,500])
                plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaT,betaT_woAR1,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
                if figure_save == 1
                    saveas(f,[home 'Figures/' char(pref) '/beta_path' '.png']);
                end
                
                betaT_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT;
                betaT_woAR1_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT_woAR1;
                beta_tilde_avg = beta_avg*((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k);
                
                %figure(401)
                figure('Name','BetaTildePath')
                set(gcf,'Position',[100,100,1200,500])
                plot_beta(var_initial,var_share,beta_tilde,beta_tilde_avg,betaT_tilde,betaT_woAR1_tilde,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
                title('β tildeの推移','FontSize',20,'FontWeight','normal','FontName',fn)
            end
        elseif beta_figure_loop == 1
            %figure(400+iTH) % Figure for BetaT with Variant Share
            figure('Name',string(['BetaPath_', sprintf(ft,TH(iTH))]))
            set(gcf,'Position',[100,100,1200,500])
            plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaT,betaT_woAR1,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
            title(string(['βの推移（',sprintf(ft,TH(iTH)),'）']),'FontSize',20,'FontWeight','normal','FontName',fn)
            
            betaT_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT;
            betaT_woAR1_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT_woAR1;
            beta_tilde_avg = beta_avg*((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k);
            
            %figure(440 + iTH)
            figure('Name',string(['BetaTildePath_', sprintf(ft,TH(iTH))]))
            set(gcf,'Position',[100,100,1200,500])
            plot_beta(var_initial,var_share,beta_tilde,beta_tilde_avg,betaT_tilde,betaT_woAR1_tilde,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
            title(string(['β tildeの推移（',sprintf(ft,TH(iTH)),'）']),'FontSize',20,'FontWeight','normal','FontName',fn)
        end
        
        
    end
    
    %minAlpha = min(minAlphaMat); %minimum alpha when variants have no effects.
    minAlpha = alpha_off; % 経済損失0 = 2020年10-11月のGDP level
    
    AlphaM = AlphaMat(~isnan(AlphaMat));
    AlphaM = (AlphaM - minAlpha)*prefGDP*10000;
    DM = DMat(~isnan(DMat));
    BackDataDA(1:length(TH),:) = [round(AlphaM'),round(DM'),TH'];
    %--- Record how many times on and off are triggered ---%
    waves = zeros(1,length(TH));
    for i = 1:length(TH)
        svec = zeros(SimPeriod-1,1);
        for t = 1:SimPeriod-1
            svec(t) = AlphaPath(t+1,i)-AlphaPath(t,i);
        end
        waves(i) = sum(svec>0);
    end
    %waves_th(1:length(waves),y,th_off_index) = waves;
    
    for l = 1:2 %1:2 when english version needed
        % Generate graphs for the website
        lng = language{l};
        %figname = 100 + l;
        %figure(figname);
        figname = string(['MainResults_' char(lng) ]);
        f = figure('Name',figname);
        set(gcf,'Position',[100,100,1200,500])
        subplot(1,2,1)
        [BackDataN,BackDataAlpha,BackDataERN] = plot_SimN(TH,TH_index,N,NPath,alpha,AlphaPath,ERN,SimERN,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,20,fn,ft,l);
        
        %--- Number of cumulative deaths ---%
        subplot(1,2,2)
        plot_Tradeoff(AlphaM,DM,waves,TH,TH_index,l,linecolor,fs,fn)
        if figure_save == 1
            saveas(f,[home 'Figures/' char(pref) '/Baseline_'  char(lng) '.png']);
            %saveas(figure(figname),[home 'Figures/' char(pref) '/MainResult_' char(lng) '.png']);
        end
    end %End of language loop = figure loop
    
    if data_save == 1
        titleN = strings(1,1+length(TH_index)*3);
        titleN(1) = "週";
        for ti = 1:length(TH_index)
            titleN(1,1+ti) = string(['新規感染者数（',sprintf(ft,TH_index(ti)),'）']);
            titleN(1,1+length(TH_index)+ti) = string(['経済活動（',sprintf(ft,TH_index(ti)),'）']);
            titleN(1,1+length(TH_index)*2+ti) = string(['実効再生産数（',sprintf(ft,TH_index(ti)),'）']);
        end
        TN = table([titleN;MonthWeekJP(Tdata-7:end-1),round(BackDataN(:,1:length(TH_index))/7),round(100*(1-BackDataAlpha(:,1:length(TH_index))),1),round(BackDataERN(:,1:length(TH_index)),2)]);
        titleAD = ["経済損失（億円）","死亡者数","ケース"];
        TAD = table([titleAD;BackDataDA(1:length(TH),:)]);
        
        writetable(TN,[home 'Figures/' char(pref) '/BackData_' char(pref)  '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
        writetable(TAD,[home 'Figures/' char(pref) '/BackData_' char(pref) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
    end
    
    %figname = 140 + pindex; %Plotting New Cases + Trade Off + Alpha Path
    %figure(figname)
    figname = string(['MainResult+Alpha_' char(pref)]);
    figure('Name',figname)
    set(gcf,'Position',[100,100,1400,600])
    subplot(1,3,1)
    plot_SimN(TH,TH_index,N,NPath,alpha,AlphaPath,ERN,SimERN,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,20,fn,ft,l);
    subplot(1,3,2)
    plot_Tradeoff(AlphaM,DM,waves,TH,TH_index,2,linecolor,fs,fn)
    subplot(1,3,3)
    plot_Alpha(alpha,AlphaPath,TH,TH_index,MonthWeekEN,WeekNumber,Tdata,linecolor,ft,fs,fn,2)

    
%     end
end %end of prefecture loop


% %% Figures added on 2021/03/25
% cumDataEffV = effectiveness*cumsum(V2_w);
% cumDataEffV = E1*cumsum(V1_w) + (E2-E1)*cumsum(V2_w);
% cumSimEffV = cumsum(V)+cumDataEffV(end);
% cumEffV  = [cumDataEffV; cumSimEffV];
%
% long_delta = [delta;deltaT];
% long_beta = [beta;betaT];
% long_gamma = [ones(Tdata,1)*gamma;gammaT];
% for i=1:length(TH)
%     long_alpha(:,i) = [alpha;AlphaPath(:,i,SimCases)];
%     long_beta_tilde(:,i) = [beta_tilde;betaT.*(1+h(2)/h(1).*AlphaPath(:,i,SimCases)).^k];
%     long_ERN(:,i) = [ERN;SimERN(:,i,SimCases)];
%     long_GDP(:,i) = [GDP;(1-AlphaPath(:,i,SimCases)).*referenceGDP(Tdata+1:Tdata+SimPeriod)];
%     long_N(:,i) = [N;NPath(:,i,SimCases)];
% end
% long_S = [S(1:end);SimData(2:end,1,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
% long_I = [I(1:end);SimData(2:end,2,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
% long_gamI = [gamma*I(2:end);gammaT.*SimData(2:end,2,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
% long_R = [R(1:end);SimData(2:end,3,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
% long_D = [D(1:end);SimData(2:end,4,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
% long_dD = [D(2:end)-D(1:end-1);SimData(2:end,4,(abs(TH - TH_index(1)) < 0.0001),SimCases)-SimData(1:end-1,4,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
%
% figure(1111)
% set(gcf,'Position',[100,100,1000,800])
% subplot(3,3,1)
% plot(long_S)
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('Suceptibles')
% subplot(3,3,2)
% plot(long_N(:,(abs(TH - TH_index(1)) < 0.0001)))
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('New Cases')
% subplot(3,3,3)
% plot(long_I)
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('Infected')
% subplot(3,3,4)
% plot(long_gamI)
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('New Recovery')
% subplot(3,3,5)
% plot(long_dD)
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('New Death')
% subplot(3,3,6)
% plot(long_R)
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('Cumulative Recovery')
% subplot(3,3,7)
% plot(long_D)
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('Cumulative Death')
% subplot(3,3,8)
% plot(cumEffV)
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('Cumulative Effective Vaccinated')
%
% figure(1112)
% set(gcf,'Position',[100,100,1000,800])
% subplot(3,3,1)
% plot(long_delta)
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('Delta')
%
% subplot(3,3,2)
% plot(long_beta)
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('Beta')
%
% subplot(3,3,3)
% for i=1:length(TH)
%     plot(long_alpha(:,i))
%     hold on
% end
% hold off
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('Alpha')
%
% subplot(3,3,4)
% for i = 1:length(TH)
%     plot(long_GDP(:,i))
%     hold on
% end
% hold off
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('GDP to Reference GDP')
%
% subplot(3,3,5)
% for i = 1:length(TH)
%     plot(long_ERN(:,i))
%     hold on
% end
% hold off
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('ERN')
%
% subplot(3,3,6)
% for i = 1:length(TH)
%     plot(long_N(:,i))
%     hold on
% end
% hold off
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('New Cases')
%
% subplot(3,3,7)
% plot(long_delta(30:end))
% xline(Tdata-30+1,'LineWidth',1.5,'HandleVisibility','off');
% title('Delta since 31st week')
%
% subplot(3,3,8)
% plot(long_beta(30:end))
% xline(Tdata-30+1,'LineWidth',1.5,'HandleVisibility','off');
% title('Beta since 31st week')
%
% subplot(3,3,9)
% for i = 1:length(TH)
%     plot(long_ERN(30:end,i))
%     hold on
% end
% hold off
% xline(Tdata-30+1,'LineWidth',1.5,'HandleVisibility','off');
% title('ERN since 31st week')
%
% figure(1114)
% Tbound = 30;
% set(gcf,'Position',[100,100,1200,800])
% subplot(4,2,1)
% plot(long_delta(1:Tbound))
% title('Delta (1-30 week)')
%
% subplot(4,2,2)
% plot(long_beta(1:Tbound))
% hold on
% plot(long_beta_tilde(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)))
% plot(long_N(1:Tbound,(abs(TH - TH_index(1)) < 0.0001))./long_I(1:Tbound))
% plot(POP0./long_S(1:Tbound,1))
% hold off
% legend('beta','beta tilde','N/I','S0/St')
% lgd = legend;
% lgd.Location = 'southeast';
% title('Beta (1-30 week)')
%
% subplot(4,2,3)
% yyaxis left
% plot(long_N(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)))
% hold on
% plot(long_I(1:Tbound))
% ylabel('# of People')
% hold on
% yyaxis right
% ylabel('New Death')
% plot(long_dD(1:Tbound))
% legend('N','I','D')
% title('New Cases / Infected / Cumlative Death (1-30 week)')
%
% subplot(4,2,4)
% yyaxis left
% plot(long_alpha(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)))
% hold on
% ylabel('Alpha')
% yyaxis right
% plot(long_ERN(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)))
% ylabel('ERN')
% legend('alpha','ERN')
% title('Alpha (1-30 week)')
%
% subplot(4,2,5)
% plot(long_delta(Tbound+1:end))
% xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
% title('Delta since 31st week')
%
% subplot(4,2,6)
% plot(long_beta(Tbound+1:end))
% hold on
% plot(long_beta_tilde(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)))
% plot(long_N(Tbound+1:Tdata,(abs(TH - TH_index(1)) < 0.0001))./long_I(Tbound+1:Tdata))
% plot(POP0./long_S(Tbound+1:Tdata,1))
% hold off
% legend('beta','beta tilde','N/I','S0/St')
% lgd = legend;
% lgd.Location = 'southeast';
% xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
% ylim([0.5 inf])
% title('Beta since 31st week')
%
% subplot(4,2,7)
% yyaxis left
% plot(long_N(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)))
% hold on
% plot(long_I(Tbound+1:end))
% ylabel('# of People')
% hold on
% yyaxis right
% ylabel('New Death')
% plot(long_dD(Tbound+1:end))
% legend('N','I','dD')
% xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
% title('New Cases / Infected / New Death since 31st week')

% subplot(4,2,8)
% yyaxis left
% plot(long_alpha(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)))
% hold on
% ylabel('Alpha')
% yyaxis right
% plot(long_ERN(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)))
% ylabel('ERN')
% legend('alpha','ERN')
% xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
% title('Alpha since 31st week')

%
% figure(1115)
% Tbound = 30;
% set(gcf,'Position',[100,100,1400,600])
% subplot(2,3,1)
% yyaxis left
% plot(long_delta(1:Tbound),'LineWidth',2)
% hold on
% yyaxis right
% plot(long_beta(1:Tbound),'LineWidth',2)
% plot(long_ERN(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)),'-k','LineWidth',2)
% legend('delta','beta','ERN')
% lgd = legend;
% lgd.Location = 'northwest';
% title('Delta (1-30 week)')
%
% subplot(2,3,2)
% yyaxis left
% plot(long_N(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)),'LineWidth',2)
% hold on
% plot(long_I(1:Tbound),'-k','LineWidth',2)
% ylabel('# of People')
% hold on
% yyaxis right
% ylabel('New Death')
% plot(long_dD(1:Tbound),'LineWidth',2)
% legend('N','I','D')
% lgd = legend;
% lgd.Location = 'northwest';
% title('New Cases / Infected / Cumlative Death (1-30 week)')
%
% subplot(2,3,3)
% yyaxis left
% plot(long_alpha(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)),'LineWidth',2)
% hold on
% yyaxis right
% plot(long_beta(1:Tbound),'LineWidth',2)
% hold on
% plot(long_beta_tilde(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)),'-k','LineWidth',2)
% legend('alpha','beta','beta tilde')
% lgd = legend;
% lgd.Location = 'southeast';
% title('Alpha (1-30 week)')
%
% subplot(2,3,4)
% yyaxis left
% plot(long_delta(Tbound+1:end),'LineWidth',2)
% hold on
% yyaxis right
% plot(long_beta(Tbound+1:end),'LineWidth',2)
% plot(long_ERN(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)),'-k','LineWidth',2)
% xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
% legend('delta','beta','ERN')
% lgd = legend;
% lgd.Location = 'northeast';
% title('Delta since 31st week')
%
% subplot(2,3,5)
% yyaxis left
% plot(long_N(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)),'LineWidth',2)
% hold on
% plot(long_I(Tbound+1:end),'-k','LineWidth',2)
% ylabel('# of People')
% hold on
% yyaxis right
% ylabel('New Death')
% plot(long_dD(Tbound+1:end),'LineWidth',2)
% legend('N','I','dD')
% xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
% title('New Cases / Infected / Cumlative Death since 31st week')
%
% subplot(2,3,6)
% yyaxis left
% plot(long_alpha(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)),'LineWidth',2)
% hold on
% yyaxis right
% plot(long_beta(Tbound+1:end),'LineWidth',2)
% hold on
% plot(long_beta_tilde(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)),'-k','LineWidth',2)
% xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
% legend('alpha','beta','beta tilde')
% title('Alpha since 31st week')










%alpha Tdata AverageAlpha CumD LagResults
% alpha = FileData0413.alpha;
% Tdata = FileData0413.Tdata;
% AverageAlpha = FileData0413.AverageAlpha;
% CumD = FileData0413.CumD;
% LagResults = FileData0413.LagResults;
% save('testJapan20210413.mat','alpha','Tdata','AverageAlpha','CumD','LagResults');