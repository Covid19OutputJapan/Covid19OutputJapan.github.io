% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
iPC=0; % 0 for Mac, 1 for Windows
if iPC==1
    %home = '\Users\masam\Dropbox\fujii_nakata\Policy_request\Cabinet_2021APR28_DR\';
    %home = 'C:\Users\tak10\Dropbox\fujii_nakata\Website\Codes\Cabinet_2021APR28_DR\';
else
    home ='/Users/machikohei/Dropbox/fujii_nakata/Website/Codes/Cabinet_2021APR28_DR/';
    %home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/Cabinet_2021APR28_DR/';
    %home = '/Users/ymaeda/Documents/ym_doc/Tokyo_Univrsity_MA/Reserach Assistant/fujii_nakata/Codes/';
    %home = '/Users/shotaro/Dropbox/fujii_nakata/Website/Codes/Cabinet_2021APR28_DR/';
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 0;      % save back data
% in the "Figure" folder
fs = 16;            % common font size for many figures
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
paces_ori = 3600000;
gradual_paces = 10;
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
th_off_vector = [500,100,110,80,825,50,60]; %1回目の緊急事態宣言の解除基準
th_off_vector2 = [400,100,110,80,500,50,60]; %2回目の緊急事態宣言の解除基準
th_off_vector3 = [400,100,110,80,500,50,60]; %3回目の緊急事態宣言の解除基準
% 経済回復速度
DRi_vector = [6,6,6,6,6,6,6,6];
% 緊急事態宣言の強さの基準: Simulation開始時のERNをいくつにするか
ERN_on_vector =  [0.55,0.5,0.5,0.5,0.65]; % 緊急事態宣言下のERN
% Size of AR(1) shock for beta process
betaT_temp_ini_vec = [0.2,0,0,0,0.15,0,0,0];
% Initial Share of Variants (Need to change these values every week)
var_initial_vector = [0.2531, 0.2909, 0.1335, 0.0909, 0.7721, 0, 0, 0,0]; % 4/26
share_index = 2;
%=========== use this code to extrapolate var_initial_vector ==================%
var_initial_vector = ini2now_infection_rate(var_initial_vector,var_growth,SimPeriod,share_index);
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
var_infection_vec = [0.3, 0.5];

% Vaccine Paces
paces_vector = [3600000,7000000];

% Different threshold for lifting the state of emergency
TL_vector = {100:50:500,70:10:100,80:10:130,50:10:100,100:50:500,20:10:100,10:10:100,50:10:100}; % 解除基準分析をコントロールしている Cell array
TL_index_vector = {[250,500,100],[100,80],[130,110],[100,80],[250,500,100],[50,20],[60,10]}; %色付けして強調するものを決める

for pindex = 1:1 %:length(PrefVector) %change this parameter for prefecture
    %====================== Model parameter values ======================%
    pref = PrefVector{pindex};        % prefecture to be analyzed
    prefGDP = GDPVector(pindex);
    %====================================================================%
    %--- Import data ---%
    % Covid data are recorded at weekly frequencies (Mon-Sun)
    % The first week start on January 20 (Mon), 2020
    if iPC==1
        covid = importdata([home '\Covid_weekly.csv']);  % Import weekly Covid data by prefecture
    else
        covid = importdata([home 'Covid_weekly.csv']);  % Import weekly Covid data by prefecture
    end
    Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
    % Columns: 1 = date, 2 = new positive, 3 = death, 4 = mobility,
    % 5 = GDP, 6 = population, 7 = GDO
    dateD = Data(:,1) + 21916;
    N = Data(:,2);
    dD = Data(:,3);
    M = Data(:,4);
    % GDP = Data(:,5);
    POP = Data(:,6);
    GDP = Data(:,7);
    GDP = GDP/GDP(1)*100;   % Normalization: 100 = 2020JAN
    Tdata= size(Data,1);    % Data period in weeks
    POP0 = POP(1);          % initial population
    ps = POP0/POP_jp;    % population share
    xtick1 = 1:13:Tdata;
    dateEN = datetime(dateD,'ConvertFrom','excel');
    SimDate = dateD(end)+7:7:dateD(end)+7*SimPeriod;
    SimDateEN = datetime(SimDate,'ConvertFrom','excel');
    %--- Create Month-Week labels ---%
    dateP = dateD(end)+7:7:dateD(end)+7*(SimPeriod+1);
    date = [dateD;dateP'];
    date = datetime(date,'ConvertFrom','excel');
    MonthNumber = month(date);
    WeekNumber = zeros(length(MonthNumber),1);
    WeekNumber(1:2) = [3;4];
    for i = 3:length(WeekNumber)
        if MonthNumber(i)~=MonthNumber(i-1)
            WeekNumber(i) = 1;
        else
            WeekNumber(i) = WeekNumber(i-1) + 1;
        end
    end
    MonthWeekJP = strings([length(MonthNumber),1]);
    MonthWeekEN = strings([length(MonthNumber),1]);
    for i = 1:length(MonthNumber)
        MonthWeekJP(i) = [num2str(date(i).Month) '月第' num2str(WeekNumber(i)) '週'];
        MonthWeekEN(i) = [datestr(date(i),'mmm') '-' num2str(WeekNumber(i)) 'w'];
    end
    M = 1+0.01*M;
    TdataGDP = Tdata-sum(isnan(GDP));
    if retroH_switch == 1
        RetroH = TdataGDP -4;
    end
    if isempty(find(SimDateEN == medical_start_date,1)) == 0
        medical_start = find(SimDateEN == medical_start_date);
    else
        medical_start = 1;
    end
    elderly_start = find(SimDateEN == elderly_start_date);
    VacStart = find(SimDateEN == datetime(2021,4,1));
    End2020 = find(dateEN == datetime(2021,1,7));
    %--- Construct weekly vaccine data ---%　関数化？
    if iPC == 1
        vaccine = importdata([home 'vaccine_daily.xls']);   % Import daily vaccine data
        dateV = datetime(vaccine.textdata(2:end,1),'InputFormat','yyyy/MM/dd');
    else
        vaccine = importdata([home 'vaccine_daily.xls']);
        dateV = datetime(vaccine.textdata(2:end,1),'InputFormat','yyyy/MM/dd');
    end
    %vaccine = importdata([home 'vaccine_daily.xls']);
    vaccine_data = vaccine.data;
    totalV_d = vaccine_data(:,1);
    V1_d = vaccine_data(:,2);
    V2_d = vaccine_data(:,3);
    V1_d(isnan(V1_d)) = 0;
    V2_d(isnan(V2_d)) = 0;
    Vdata = size(vaccine,1);
    V1_w = zeros(length(dateEN),1);
    V2_w = zeros(length(dateEN),1);
    for i = 1:1:length(dateEN)
        index_date = (dateV >= dateEN(i)-3 & dateV <= dateEN(i)+3);
        V1_w(i) = sum(V1_d(index_date));
        V2_w(i) = sum(V2_d(index_date));
    end
    V1_w = ps*V1_w;
    V2_w = ps*V2_w;
    %--- Constructing the reference level of output ---%　関数化？
    potentialGDP = zeros(52*3,1);
    %potentialGDP(1) = (100/(1.0122))*(1.0063^(1/12)) %GDP in 2019/12 = 100
    potentialGDP(1) = 100;
    for i = 2:length(potentialGDP)
        if i <= 13
            potentialGDP(i) = potentialGDP(i-1)*(1.0063^(1/52));
        elseif i <= 52
            potentialGDP(i) = potentialGDP(i-1)*(1.0021^(1/52));
        elseif i <= 104
            potentialGDP(i) = potentialGDP(i-1)*(1.0024^(1/52));
        elseif i <= 156
            potentialGDP(i) = potentialGDP(i-1)*(1.0021^(1/52));
        end
    end
    referenceGDP = potentialGDP.*(1+0.0166);
    referenceGDP(1:2) = [];
    
    alpha = (1 - GDP(1:TdataGDP)./referenceGDP(1:TdataGDP));   % output loss in percentage
    %--- Regress mobility on alpha to estimate the elasticity h ---%　関数化？
    Malt=M;
    Malt(50)=0.5*(Malt(49)+Malt(51));
    Y = Malt(4:TdataGDP);
    X = alpha(4:TdataGDP);
    if hconstant == 0
        Y = Y - 1;
        h_all = (X'*X)\X'*Y;              % OLS estimate of h
        reg = X*h_all;
        r = Y - reg;                % r is the residuals, which is the observed minus fitted values
        SSE = sum(r.^2);            % SSE is the sum of squared errors
        MSE=SSE/(length(Y)-1);      % mean squared error
        h_all_se=sqrt(MSE/sum(X.^2));   % standard error of h
    elseif hconstant == 1
        XC = [ones(length(X),1), X];
        h_all = (XC'*XC)\XC'*Y;         % OLS estimate of h with constant
        reg = XC*h_all;
        r = Y - reg;
        SSE = sum(r.^2);
        MSE=SSE/(length(Y)-1);      % mean squared error
        h_all_se = zeros(2,1);
        h_all_se(1)=sqrt(MSE/sum(XC(:,1).^2));   % standard error of h
        h_all_se(2)=sqrt(MSE/sum(XC(:,2).^2));
    end
    
    Y = Malt(TdataGDP-RetroH:TdataGDP);
    X = alpha(TdataGDP-RetroH:TdataGDP);
    if hconstant == 0
        Y = Y - 1;
        h = (X'*X)\X'*Y;              % OLS estimate of h
        reg = X*h;
        r = Y - reg;                % r is the residuals, which is the observed minus fitted values
        SSE = sum(r.^2);            % SSE is the sum of squared errors
        MSE=SSE/(length(Y)-1);      % mean squared error
        h_se=sqrt(MSE/sum(X.^2));   % standard error of h
    elseif hconstant == 1
        XC = [ones(length(X),1), X];
        h = (XC'*XC)\XC'*Y;         % OLS estimate of h with constant
        reg = XC*h;
        r = Y - reg;
        SSE = sum(r.^2);
        MSE=SSE/(length(Y)-1);      % mean squared error
        h_se = zeros(2,1);
        h_se(1)=sqrt(MSE/sum(XC(:,1).^2));   % standard error of h
        h_se(2)=sqrt(MSE/sum(XC(:,2).^2));
    end
    %--- Plot mobility data ---%
    figure(2)
    yyaxis left
    plot(Malt,'k','LineWidth',1.5)
    ylabel('Mobility')
    yyaxis right
    hold on
    plot((1-alpha)*100,'r-.','LineWidth',1.5)
    hold on
    plot((1-alpha(1:TdataGDP))*100,'b-.','LineWidth',1.5)
    ylabel('GDP')
    xlim([1 Tdata])
    xticks(xtick1)
    xticklabels(MonthWeekJP(xtick1))
    legend('Mobility (L axis)','Imputed GDP (R axis)', 'Data GDP (R axis)','FontSize',16,'Location','southeast');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'b';
    ax.YAxis(1).FontSize = fs;
    ax.YAxis(2).FontSize = fs;
    ax.XAxis.FontSize = fs;
    xtickangle(45)
    if figure_save == 1
        saveas(figure(2),[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
    end
    
    %--- Compute the history of S, I, R, D in the data period ---%
    S = zeros(Tdata+1,1);
    I = zeros(Tdata+1,1);
    R = zeros(Tdata+1,1);
    D = zeros(Tdata+1,1);
    S(1)=POP0;
    for i = 1:Tdata
        S(i+1)=S(i)-N(i)-E1*V1_w(i)-(E2-E1)*V2_w(i);
        I(i+1)=I(i)+N(i)-gamma*I(i)-dD(i);
        R(i+1)=R(i)+gamma*I(i)+E1*V1_w(i)+(E2-E1)*V2_w(i);
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% Projection starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------------------------------------------------ %
    % --------------- Main analysis starts here --------------- %
    % ------------------------------------------------ %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Projection parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DRi = DRi_vector(pindex);
    th_on1 = th_on_vector(pindex)*7;         % threshold to place the state of emergency
    th_on2 = th_on_vector2(pindex)*7;         % threshold to place the state of emergency
    th_off1 = th_off_vector(pindex)*7;
    th_off2 = th_off_vector2(pindex)*7;
    th_off3 = th_off_vector3(pindex)*7;
    ERN_on = ERN_on_vector(pindex);
    
    % Change both paces and infection rate at the same time
%     TH = [3600000,7000000,3600000,7000000];
%     TH_index = [3600000,7000000,3600000,7000000];
%     TH2 = [0.3, 0.3 ,0.5, 0.5];
%     TH2_index = [0.3, 0.5];
    
    TH = 0.3:0.1:0.8;
    TH_index = [0.5,0.6,0.8];
    
    
    %TH = paces_vector;
    %TH_index = paces_vector;

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
    
    for iTH = 1:length(TH)
        %%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%
%         paces_ori = TH(iTH);
%         var_infection = TH2(iTH);
          ERN_on = TH(iTH)

    betaT_temp_ini = betaT_temp_ini_vec(pindex);
    var_initial = var_initial_vector(pindex);
    
    alpha_off = mean(alpha((dateEN >= datetime(2020,10,1)) & (datetime(2020,11,26)>= dateEN ))); % output loss without the state of emergency
    InitialValues = [S(end),I(end),R(end),D(end)];
    
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
    % Original
    [V,deltaT,VT] = vaccine_distribution_medical(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2,ps,POP0,3);
    %[V,deltaT,VT] = vaccine_distribution_simple(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2);
    delta_ss = delta_average*(0.09/1.28);
    delta_ss = (delta_average - delta_ss)*(1-accept_share*D2)+delta_ss; %death rate after full vaccination of elderly = 2021/09/30
    
    
    
   
    
    %% figure for vaccine path (関数化)
    
    figure(200)
    set(gcf,'Position',[100,100,1200,800])
    Data_vac = VT(:,:);
    medStuff1_nv = [V1_w; Data_vac(:,5)]; % In future, you should change V1_w to V1past_med.
    medStuff2_nv = [V2_w; Data_vac(:,6)]; % In future, you should change V2_w to V2past_med.
    elderly1_nv = [zeros(Tdata,1); Data_vac(:,1)];
    elderly2_nv = [zeros(Tdata,1); Data_vac(:,2)];
    others1_nv = [zeros(Tdata,1); Data_vac(:,3)];
    others2_nv = [zeros(Tdata,1); Data_vac(:,4)];
    % Find the cumulative number of vaccines for each agent
    medStuff1_cv = cumsum(medStuff1_nv);
    medStuff2_cv = cumsum(medStuff2_nv);
    elderly1_cv = cumsum(elderly1_nv);
    elderly2_cv = cumsum(elderly2_nv);
    others1_cv = cumsum(others1_nv);
    others2_cv = cumsum(others2_nv);
    Area_nv = [medStuff1_nv,medStuff2_nv,elderly1_nv,elderly2_nv,others1_nv,others2_nv];
    Area_cv = [medStuff1_cv,medStuff2_cv,elderly1_cv,elderly2_cv,others1_cv,others2_cv];
    Area_nv = Area_nv / (ps*10000);
    Area_cv = Area_cv / (ps*100000000);
    % Area_nv = [medStuff1_nv;medStuff2_nv;elderly1_nv;elderly2_nv;others1_nv;others2_nv];
    % Area_cv = [medStuff1_cv;medStuff2_cv;elderly1_cv;elderly2_cv;others1_cv;others2_cv];
    
    for ifig = 1:1:2
        if ifig == 1
            subplot(3,2,1)
            %         figure(71)
            Areagraph = bar(Area_nv,'stacked','LineStyle','none');
        else
            subplot(3,2,2)
            %         figure(72)
            Areagraph = area(Area_cv,'LineStyle','none');
            %         Areagraph = area(transpose(Area_cv),'LineStyle','none');
            Areagraph = bar(Area_cv,'stacked','LineStyle','none');
        end
        xtickangle(45)
        xticks(find(WeekNumber==1))
        xlim([Tdata-7 Tdata+56])
        xticklabels(MonthWeekEN(WeekNumber==1))
        
        ldfs = fs/3;
        if ifig == 1
            ylim([0 800])
            ylabel('ワクチン本数（万本）','FontName',fn)
            %             legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3),Areagraph(2),Areagraph(1)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','医療従事者２本目','医療従事者１本目','FontSize',ldfs,'Location','NorthEast','FontName',fn);
            legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','FontSize',ldfs,'Location','NorthEast','FontName',fn);
            title('新規ワクチン接種本数（週ごと）', 'FontSize',fs,'FontName',fn);
        else
            ylim([0 2.5])
            ylabel('ワクチン本数（億本）','FontName',fn)
            %             legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3),Areagraph(2),Areagraph(1)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','医療従事者２本目','医療従事者１本目','FontSize',ldfs,'Location','NorthWest','FontName',fn);
            legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','FontSize',ldfs,'Location','NorthWest','FontName',fn);
            title('累計ワクチン接種本数（週ごと）', 'FontSize',fs,'FontName',fn);
        end
        
        ax = gca;
        ax.YAxis.Color = 'k';
        ax.YAxis.FontSize = fs;
        ax.XAxis.FontSize = fs;
        Areagraph(3).FaceColor = [0.2 0.6 0.8];
        Areagraph(1).FaceColor = [0.6 0.6 0.6];
        Areagraph(5).FaceColor = [0.4 0.4 0.8];
        Areagraph(4).FaceColor = [0.2 0.6 0.8];
        Areagraph(2).FaceColor = [0.6 0.6 0.6];
        Areagraph(6).FaceColor = [0.4 0.4 0.8];
        Areagraph(1).FaceAlpha = 1;
        Areagraph(3).FaceAlpha = 0.7;
        Areagraph(5).FaceAlpha = 1;
        Areagraph(2).FaceAlpha = 0.8;
        Areagraph(4).FaceAlpha = 0.5;
        Areagraph(6).FaceAlpha = 0.8;
    end
    subplot(3,2,3)
    plot([delta; deltaT]./delta_average,'-k','LineWidth',1.5);
    xtickangle(45)
    xticks(find(WeekNumber==1))
    xlim([Tdata+1 Tdata+56])
    xticklabels(MonthWeekEN(WeekNumber==1))
    title('致死率（現在のレベルで標準化）','FontName',fn)
    ax = gca;
    ax.YAxis.Color = 'k';
    ax.YAxis.FontSize = fs;
    ax.XAxis.FontSize = fs;
    legend;
    
    %compute_minALPHA;
    
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        beta_r = 0;
        for retrop = retro_lb:retro_ub
            beta_r = beta_r + mean(beta(end-retrop+1:end));
        end
        beta_avg = beta_r/(retro_ub-retro_lb+1);
        betaT = beta_avg*ones(SimPeriod,1);
        
        %for iS = 1:length(D1_vec)
        
        var_intercept = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
        % var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept
        
        % step 1: extrapolation in the variant share in the last 17 weeks     -mean(retro_lb,retro_ub)+1 = -17 + 1
        var_share = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
        var_share_prev = exp((-mean(retro_lb,retro_ub)+1:0)'*var_growth+var_intercept)./(1+exp((-mean(retro_lb,retro_ub)+1:0)'*var_growth+var_intercept));
        % step 2: relative infection rate in the last 17 weeks
        relative_var_infection_prev = 1 + var_share_prev * var_infection;
        % step 3:
        no_var_beta = beta(end-mean(retro_lb,retro_ub)+1:end) ./ relative_var_infection_prev;
        beta_bar = mean(no_var_beta);
        betaT = beta_bar*(1+var_infection*var_share);
        betaT_woAR1 = betaT;
        
        betaT_temp = betaT_temp_ini;
        betaT(1,1) = (1+betaT_temp)*betaT(1,1);
        for i = 2:length(betaT)
            betaT_temp = betaT_temp * beta_rho;
            betaT(i) = betaT(i) * (1+betaT_temp);
        end
                
        altA_on = (((ERN_on*(POP0/S(end))*((gammaT(1)+delta_average)/beta_avg)).^(1/k))-1)*(h(1)/h(2));
        ERNCheck = (S(end)/POP0).*(((1+(h(2)/h(1))*alpha_off).^k).*beta_avg)./(gammaT(1)+delta_average);
                
        %---- Simulation ---%
        if pindex < 5
            [DMat(iTH),AlphaMat(iTH),AlphaPath(:,iTH),SimData(:,:,iTH),NPath(:,iTH),SimERN(:,iTH)] ...
                = Covid_projection_control_gradual_off_threshold_delta2...
                (InitialValues,altA_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,betaT,gammaT,deltaT,delta_ss,V,h,k,POP0,hconstant,DRi,alpha(end));
        else
            [DMat(iTH),AlphaMat(iTH),AlphaPath(:,iTH),SimData(:,:,iTH),NPath(:,iTH),SimERN(:,iTH)] ...
                = Covid_projection_control_gradual_threshold_delta2...
                (InitialValues,altA_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,betaT,gammaT,deltaT,delta_ss,V,h,k,POP0,hconstant,DRi);
        end
        
%         
%         plot_var_share(1,1) = var_initial;
%         plot_var_share(2:length(SimDate)+1,1) = var_share(:,1);
%         tt = 12; %Showing previous t periods for the plot
%         
%         plot_betaT(1:tt,1) = beta(end-tt+1:end);
%         plot_betaT(tt+1:length(SimDate)+tt,1) = betaT(:,1);
%         
%         plot_betaT_woAR1(1:tt,1) = beta(end-tt+1:end);
%         plot_betaT_woAR1(tt+1:length(SimDate)+tt,1) = betaT_woAR1(:,1);
%         
%         plot_betaT3(1:tt,1) = beta(end-tt+1:end);
%         plot_betaT3(tt+1:length(SimDate)+tt,1) = beta_avg*ones(length(SimDate),1);
%         
%         xaxis_vec = 0:1:length(SimDate);
%         xaxis_vec2 = 1:1:length(SimDate)+tt;
%         xaxis_vec3 = 1:1:Tdata+tt;
%         
%         MonthWeekEN_Sim = MonthWeekEN(length(dateD)+1:end,1);
%         MonthNumber_Sim = MonthNumber(length(dateD)+1:end,1);
%         WeekNumber_Sim = WeekNumber(length(dateD)+1:end,1);
%         MonthWeekJP_Sim = MonthWeekJP(length(dateD)+1:end,1);
%         
%         MonthWeekEN_Sim2 = MonthWeekEN(length(dateD)-tt+1:end,1);
%         MonthNumber_Sim2 = MonthNumber(length(dateD)-tt+1:end,1);
%         WeekNumber_Sim2 = WeekNumber(length(dateD)-tt+1:end,1);
%         MonthWeekJP_Sim2 = MonthWeekJP(length(dateD)-tt+1:end,1);
%         
%         figure(4000)
%         set(gcf,'Position',[100,100,1200,500])
%         subplot(1,2,1)
%         plot(xaxis_vec,plot_var_share(:,1)*100,'-r','LineWidth',2)
%         ax = gca;
%         ax.YAxis.FontSize = 20;
%         ax.XAxis.FontSize = 20;
%         ax.YAxis.Exponent = 0;
%         ytickformat('%,3.0f')
%         xticks(find(WeekNumber_Sim==1))
%         title('新規感染者数の推移','FontSize',20,'FontWeight','normal')
%         title('変異株シェアの推移','FontSize',20,'FontWeight','normal','FontName',fn)
%         xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
%         xlabel('週','FontSize',20,'FontName',fn)
%         ylabel('変異株割合(%)','FontSize',20,'FontName',fn)
%         xtickangle(45)
%         xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
%         xlim([0,length(SimDate)]);
%         % lgd = legend;
%         % lgd.Location = 'northwest';
%         subplot(1,2,2)
%         plot(xaxis_vec2(tt+1:end), plot_betaT(tt+1:end,1),'-r','LineWidth',2)
%         hold on
%         plot(xaxis_vec2(tt+1:end), plot_betaT_woAR1(tt+1:end,1),'-r','LineWidth',2.0)
%         plot(xaxis_vec2,plot_betaT3,'-k','LineWidth',2)
%         ax = gca;
%         ax.YAxis.FontSize = 20;
%         ax.XAxis.FontSize = 20;
%         ax.YAxis.Exponent = 0;
%         ytickformat('%,0.2f')
%         xticks(find(WeekNumber_Sim2==1))
%         xtickangle(45)
%         legend('感染率(変異株成長ケース)','過去の平均','FontSize',10,'FontName',fn)
%         xticklabels(MonthWeekJP_Sim2(WeekNumber_Sim2==1))
%         title('βの推移','FontSize',20,'FontWeight','normal','FontName',fn)
%         xlabel('週','FontSize',20)
%         ylabel('感染率','FontSize',20,'FontName',fn)
%         xline(tt+1,'--','LineWidth',1.5,'HandleVisibility','off');
%         xlim([0,length(SimDate)+tt]);
%         lgd = legend;
%         lgd.Location = 'southeast';
%         if figure_save == 1
%             saveas(figure(4000),[home 'Figures/' char(pref) '/beta_path' '.png']);
%         end

        
    end

    minAlpha = min(AlphaMat); %関数化
        
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
        figname = 100 + l;
        figure(figname)
        set(gcf,'Position',[100,100,1200,500])
        
        subplot(1,2,1)
        if max(TH)<3
            ft = '%.2f';
        else
            ft = '%.0f';
        end
        th_wave = 1;
        for i = 1:length(TH)
            if abs(TH(i) - TH_index(th_wave)) < 0.0001
                plot([N;NPath(:,i)]/7,linecolor{th_wave},'LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
                BackDataN(:,th_wave) = [N(end-7:end);NPath(:,i)];
                BackDataAlpha(:,th_wave) = [alpha(end-7:end);AlphaPath(:,i)];
                BackDataERN(:,th_wave) = [ERN(end-7:end);SimERN(:,i)];
                th_wave = th_wave + 1;
                if th_wave > length(TH_index)
                    th_wave = length(TH_index);
                end
            else
                plot([N;NPath(:,i)]/7,'--','LineWidth',0.3,'DisplayName',sprintf(ft,TH(i)))
            end
            hold on      
        end
        plot(N/7,'k','LineWidth',2,'HandleVisibility','off')
        xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
        ax = gca;
        ax.YAxis.FontSize = 20;
        ax.XAxis.FontSize = 20;
        ax.YAxis.Exponent = 0;
        ytickformat('%,6.0f')
        xticks(find(WeekNumber==1))
        if l == 1
            title('Projected path of new cases','FontSize',20,'FontWeight','normal')
            xticklabels(MonthWeekEN(WeekNumber==1))
            %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
        elseif l == 2
            title('新規感染者数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
            xticklabels(MonthWeekJP(WeekNumber==1))
        end
        lgd = legend;
        lgd.NumColumns = 2;
        xtickangle(45)
        xlim([Tdata-7 Tdata+36])
        
        %--- Number of cumulative deaths ---%
        subplot(1,2,2)
        plot(AlphaM(waves==0),DM(waves==0),'-go','LineWidth',2,'MarkerSize',10);
        hold on
        plot(AlphaM(waves==1),DM(waves==1),'-mo','LineWidth',2,'MarkerSize',10);
        hold on
        plot(AlphaM(waves==2),DM(waves==2),'-bo','LineWidth',2,'MarkerSize',10);
        hold on
        plot(AlphaM(waves==3),DM(waves==3),'-ko','LineWidth',2,'MarkerSize',10);
        text(AlphaM,DM,string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
        hold on
        for i = 1:length(TH_index)
            scatter(AlphaM(abs(TH - TH_index(i)) < 0.0001),DM(abs(TH - TH_index(i)) < 0.0001),250,linecolor{i},'filled'); %for loop にする
            hold on
        end
        if l == 1
            xlabel('Output Loss (hundred million yen)','FontSize',20)
            ylabel('Cumulative Deaths','FontSize',20)
            title('Relationship between Covid-19 and output','FontSize',20,'FontWeight','normal')
        elseif l == 2
            xlabel('経済損失 (億円)','FontSize',20,'FontName',fn)
            ylabel('累計死亡者数','FontSize',20,'FontName',fn)
            title('コロナ感染と経済の関係','FontSize',20,'FontWeight','normal','FontName',fn)
        end
        xlim([0,inf])
        xtickangle(45)
        grid on
        ax = gca;
        ax.YAxis.FontSize = 20;
        ax.XAxis.FontSize = 20;
        ax.YAxis.Exponent = 0;
        ax.XAxis.Exponent = 0;
        ytickformat('%,6.0f')
        if figure_save == 1
            saveas(figure(figname),[home 'Figures/' char(pref) '/MainResult_' char(lng) '.png']);
        end
        
    end %End of language loop = figure loop
    
    if data_save == 1
        titleN = strings(1,1+length(TH_index)*3);
        titleN(1) = "週";
        for ti = 1:length(TH_index)
            titleN(1,1+ti) = string(['新規感染者数（',sprintf('%.0f',TH_index(ti)),'）']);
            titleN(1,1+length(TH_index)+ti) = string(['経済活動（',sprintf('%.0f',TH_index(ti)),'）']);
            titleN(1,1+length(TH_index)*2+ti) = string(['実効再生産数（',sprintf('%.0f',TH_index(ti)),'）']);
        end
        TN = table([titleN;MonthWeekJP(Tdata-7:end-1),round(BackDataN(:,1:length(TH_index))/7),round(100*(1-BackDataAlpha(:,1:length(TH_index))),1),round(BackDataERN(:,1:length(TH_index)),2)]);
        titleAD = ["経済損失（億円）","死亡者数","ケース"];
        TAD = table([titleAD;BackDataDA(1:length(TH),:)]);
        
        writetable(TN,[home 'Figures/' char(pref) '/BackData_' char(pref)  '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
        writetable(TAD,[home 'Figures/' char(pref) '/BackData_' char(pref) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
        
        
    end

end %end of prefecture loop


% %% Figures added on 2021/03/25
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
