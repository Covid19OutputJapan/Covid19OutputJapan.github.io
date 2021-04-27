% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
iPC=0; % 0 for Mac, 1 for Windows
if iPC==1
    %     home = '\Users\shcor\Dropbox\Website\Codes\';
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
    %     home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata\Website\Codes\';
    %home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/codes/';
    home = '/Users/shotaro/Dropbox/fujii_nakata/Website/codes/';
    %     home = '/Users/machikohei/Dropbox/fujii_nakata/Website/Codes/';
    %      home = '/Users/ymaeda/Documents/ym_doc/Tokyo_Univrsity_MA/Reserach Assistant/fujii_nakata/Codes/';
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 1;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 1;      % save back data
% in the "Figure" folder
fs = 16;            % common font size for many figures
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
%======================================================================%

SimCases = 4; % 1 = Different thresholds(ON); 2 = Different durations for economic recovery ; 3 = different ERN ; 4 = different infection rate % 5 = Different thresholds (OFF)
beta_option_index = [1,42]; %{1,3,42} 1 ... baseline, 3 ... 気の緩み,
% 41 ... 変異株シナリオ1(アメリカのgrowth rate), 42 ... 変異株シナリオ2(イギリスのgrowth rate)

SimPeriod = 52;        % simulation period in weeks
gamma = 7/12;          % recovery rate from Covid % Should change this to 7/12 (4/25 Kohei Machi)
k = 2;                 % exponent of (1-h*alpha)
hconstant = 1;         % 0 = without intercept, 1 = with intercept for h regression
medical_start_date = datetime(2021,3,18);
elderly_start_date = datetime(2021,5,13);
RetroPeriod = 17;      % retroactive periods used to estimate gamma and delta

PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka','Aichi','Fukuoka','Hyogo'};

GDPVector = [106,36,23,21,40,41,20,20]; % 兆円, one trillion yen (chou-yen)
% th_on_vector = [1250,700,400,350,1000,350,350]; % threshold for another emergency
% th_on_vector2 = th_on_vector*2;
th_on_vector = [750,500,400,350,750,350,350]; % 1回目の緊急時代宣言をする基準 %present.
th_on_vector2 = [2000,1000,800,700,2000,350,350]; % 2回目の緊急時代宣言をする基準
th_on_vector3 = [2000,700,400,350,2000,350,350]; % 3回目の緊急時代宣言をする基準 %full vaccination of elderly.

th_off_vector = {[250,350],[200,80],[260,110],[200,80],[400,40],[100,20],[120,10]};
th_off_vector_2 = {[250,350],[200,80],[260,110],[200,80],[400,40],[100,20],[120,10]};
th_off_vector_3 = {[250,350],[200,80],[260,110],[200,80],[400,40],[100,20],[120,10]};
th_off_vector42 = {[500,350],[100,80],[80,110],[60,80],[100,40],[50,20],[60,10]}; %Also change TL_index
th_off_vector42_2 = {[400,350],[100,80],[80,110],[60,80],[100,40],[50,20],[60,10]}; %Also change TL_index
th_off_vector42_3 = {[400,350],[100,80],[80,110],[60,80],[100,40],[50,20],[60,10]}; %Also change TL_index

% Used in Covid_projection_...m to control ERN (economic restriction) during the state of emergency
ERN_now_vector = [0.99,0.94,0.96,0.955,0.9,0.9,0.99]; % currently not used.
ERN_on_vector = [0.9,0.99,1.00,0.99,0.85,0.97,1.02]; % 1回目の緊急時代宣言した後のERN

% Control the imputation for alpha (AR1 process)
rho_vector= [0.85,0.9,0.9,0.8,0.9,0.8,0];
impute_periods_vector =[45,22,25,20,45,30,30];

% Control the moving average of beta (beta_avg = sum_{t = lb}^{ub} (1/(ub-lb + 1) sum_{x=1}^t (1/t) beta_t)
retro_ub = 17;
retro_lb = 17;

% Control AR(1) process for beta
%betaT_temp_ini = 0.289; %so that 6 weeks average = betaT(1)*1.2
betaT_temp_ini = -0.289/2; %Osaka ... twice as much as this number
beta_rho = 0.85;

% Showing previous t periods for the plot
tt = 12;

% Parameters for variants
%var_infection = 0.3;  % variant's infection rate (compared to nomral infection rate) [0.5, 0.7]
%var_infection = 0.5;
%var_infection_vector = [0.3,0.5];
var_infection_vector = 0.5;
nvi = length(var_infection_vector);
if SimCases == 4
    var_infection_vector = 0:0.1:0.7;
end
var_infection_index = [0.3,0.5,0.7]; %used in simcase == 4
var_ss = 1.0;         % steady-state share of variant
var_growth_vec = [0.47, 0.47]; %April 12nd, 2021
%=========== use this code to extrapolate var_initial_vector ==================%
%var_initial_vector = [0.099, 0.101, 0.027, 0.122, 0.633, 0, 0, 0,0]; % 4/19
var_initial_vector = [0.2531, 0.2909, 0.1335, 0.0909, 0.7721, 0, 0, 0,0]; % 4/26
for pindex = 1:5
    var_initial = var_initial_vector(pindex);
    var_growth = 0.47;
    logit_initial = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
    var_intercept = logit_initial;
    var_share = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
    var_initial_vector(pindex) = var_share(2); %var_share(2)
end
%=====================================================%
disp(var_initial_vector) %[0.4645 0.5122 0.2829 0.2038 0.8966  0  0  0  0]; % 4/26

v_scale_vec = ones(1,length(PrefVector)); %4/12 -
scaleA = 1; %Scenario A ... initial vector is multiplied by scaleA
scaleB = 1; %Scenario B ... initial vector is multiplied by scaleB

% vaccine pace
Vsimple = 0; % 0 for vaccine_distribution; 1 for vaccine_distribution_simple
PF = 1; % 0 for AZ, 1 for PF
Vgradual = 1; % 0 for flat Vpath, 1 for gradually increasing Vpath
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


% ------------------------- Experiment Vectors ----------------------- %
% Used for SimCases = 1: different threshold for declaring the state of emergency
THON_vector1 = {500:100:1000,400:100:800,200:100:600,300:50:500,800:50:1200,200:50:400,200:50:400,300:50:500};
THON_index_vector1 = {[900,700],[600,800],[300,600],[350,500],[900,1200],[250,400],[250,400],[350,500]};
THON_vector41 = {500:250:2000,400:100:800,200:100:600,300:50:500,200:50:600,200:50:400,200:50:400,300:50:500};
THON_index_vector41 = {[1250,2000],[600,800],[300,600],[350,500],[250,400],[400,600],[250,400],[350,500]};
THON_vector42 = {[900:300:1500,2000:1000:4000],400:100:800,200:100:600,300:50:500,300:100:1000,200:50:400,200:50:400,300:50:500};
THON_index_vector42 = {[1500,3000],[600,800],[300,600],[350,500],[400,600],[400,600],[300,600],[350,500]};

% Used for SimCases = 2: Duration of recovery
DR_vector = {4:2:16,0:2:16,0:2:16,0:2:16,0:2:16,0:12,0:2:16,0:2:16};
DR_index = [6,12]; %0->4,8->12 gradual 3/15

% Used for SimCases = 3: different ERN (economic restriction) during the state of emergency
ERNON_vector1 = {0.5:0.1:1.0,400:100:800,200:100:600,300:50:500,170:20:270,200:50:400,200:50:400,300:50:500};
ERNON_index_vector1 = {[1.0,0.7],[600,800],[300,600],[350,500],[190,250],[250,400],[250,400],[350,500]};
ERNON_vector42_1 = {[0.3,0.4,0.5,0.55,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74],[],[],[],...
    [0.3,0.4,0.5,0.55,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74]};
ERNON_vector42_2 = {[0.4,0.5,0.55,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74],[],[],[],...
    [0.4,0.5,0.55,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74]};
ERNON_vector42_3 = {[0.5,0.55,0.6,0.62,0.64,0.66,0.68,0.69,0.7,0.72,0.74,0.76,0.78],[],[],[],...
    [0.5,0.55,0.6,0.62,0.64,0.66,0.68,0.69,0.7,0.72,0.74,0.76,0.78]};
ERNON_index_vector42 = {[0.7,0.5],[],[],[],[0.7,0.5]};

% Used for SimCases = 5: different threshold for lifting the state of emergency
ERN_on_scenario2 = {0.55,0.5,0.5,0.5,0.65};
%ERN_on_scenario2 = {0.55,0,0,0,[0.75,0.65]};
% ERN_on_scenario42 = {[0.7,0.55,0.40],0,0,0,[0.7,0.65,0.57]};

TL_vector = {100:50:500,70:10:100,80:10:130,50:10:100,100:50:500,20:10:100,10:10:100,50:10:100}; % 解除基準分析をコントロールしている Cell array
TL_index_vector = {[250,500,100],[100,80],[130,110],[100,80],[250,500,100],[50,20],[60,10]}; %色付けして強調するものを決める

% TL_vector = {100:50:500,70:10:100,80:10:130,50:10:100,50:10:100,20:10:100,10:10:100,50:10:100};
% TL_index_vector = {[300,500],[200,80],[260,110],[200,80],[300,40],[100,20],[120,10]};

% ---------------------------------------------------%
sizemat = 20;

for pindex = 2:4  %:length(PrefVector) %change this parameter for prefecture
    %====================== Model parameter values ======================%
    pref = PrefVector{pindex};        % prefecture to be analyzed
    prefGDP = GDPVector(pindex);
    %====================================================================%
    if pindex == 5
        betaT_temp_ini = -0.289;
    end
    sizeth1 = length(th_off_vector{pindex});
    sizeth42 = length(th_off_vector42{pindex});
    sizeth = max([sizeth1, sizeth42]);
    NPath_beta = NaN(SimPeriod,sizemat,SimCases,length(beta_option_index));
    AlphaM_beta = NaN(sizemat,SimCases,length(beta_option_index));
    DM_beta = NaN(sizemat,SimCases,length(beta_option_index));
    waves_beta = NaN(sizemat,SimCases,length(beta_option_index));
    AlphaM_th = NaN(sizemat,sizeth);
    DM_th = NaN(sizemat,sizeth);
    TH_th = NaN(sizemat,SimCases,sizeth);
    waves_th= NaN(sizemat,SimCases,sizeth);
    
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
    GDP = GDP/GDP(1)*100;
    Tdata= size(Data,1);    % Data period in weeks
    POP0 = POP(1);          % initial population
    ps = POP0/125710000;    % population share
    xtick1 = 1:13:Tdata;
    dateEN = datetime(dateD,'ConvertFrom','excel');
    
    SimDate = dateD(end)+7:7:dateD(end)+7*SimPeriod;
    SimDateEN = datetime(SimDate,'ConvertFrom','excel');
    %--- Create Month-Week labels ---%
    dateP = dateD(end)+7:7:dateD(end)+7*(SimPeriod+1);
    date = [dateD;dateP'];
    date = datetime(date,'ConvertFrom','excel');
    MonthNumber = month(date);
    xaxis_vec = 0:1:length(SimDate); % used for plotting
    xaxis_vec2 = 1:1:length(SimDate)+tt; % used for plotting
    xaxis_vec3 = 1:1:Tdata+tt; % used for plotting
    
    
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
    RetroH = TdataGDP-4;         %     RetroH = 15;
    if isempty(find(SimDateEN == medical_start_date,1)) == 0
        medical_start = find(SimDateEN == medical_start_date);
    else
        medical_start = 1;
    end
    elderly_start = find(SimDateEN == elderly_start_date);
    VacStart = find(SimDateEN == datetime(2021,4,1));
    End2020 = find(dateEN == datetime(2021,1,7));
    
    %--- Construct weekly vaccine data ---%
    %             if iPC == 1
    %                 vaccine = importdata([home 'vaccine_daily.xls']);   % Import daily vaccine data
    %             else
    %                 vaccine = importdata([home 'vaccine_daily.xls']);
    %                 %vaccine = importdata([home 'vaccine_daily.csv']);   % Import daily vaccine data
    %             end
    %             dateV = datetime(vaccine.textdata(2:end,1),'InputFormat','yyyy/MM/dd');
    
    vaccine = importdata([home 'vaccine_daily.xls']);
    dateV = datetime(vaccine.textdata(2:end,1),'InputFormat','MM/dd/yyyy');
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
    
    
    %--- Constructing the reference level of output ---%
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
    % gapGDP = zeros(length(potentialGDP),1);
    % gapGDP(1) = 0.0166;
    % for i = 2:length(gapGDP)
    %     gapGDP(i) = gapGDP(i-1)*(0.975^(12/52));
    % end
    referenceGDP = potentialGDP.*(1+0.0166);
    referenceGDP(1:2) = [];
    
    
    %--- Impute alpha (regress alpha on M)---%
    
    Malt=M;
    Malt(50)=0.5*(Malt(49)+Malt(51));
    alpha = (1 - GDP(1:TdataGDP)./referenceGDP(1:TdataGDP));   % output loss in percentage
    X = Malt(TdataGDP-impute_periods_vector(pindex):TdataGDP);
    Y = alpha(TdataGDP-impute_periods_vector(pindex):TdataGDP);
    XC = [ones(length(X),1), X];
    s = (XC'*XC)\XC'*Y;         % OLS estimate of h with constant
    reg = XC*s;
    r_p = Y - reg; %residual
    SSE_p = sum(r_p.^2);
    eps_p = zeros(Tdata-TdataGDP,1);
    
    
    %Variable eps_p
    eps_p(1) = r_p(end);
    for i = 1:Tdata-TdataGDP-1
        eps_p(i+1) = rho_vector(pindex)*eps_p(i);
    end
    
    %     %Constant eps_p
    %     eps_p_adjust= 0; %r_p(end);
    %     for i = 1:Tdata-TdataGDP
    %         eps_p(i)=eps_p(i)+eps_p_adjust;
    %     end
    
    alpha_pred = s(1)+s(2)*Malt(TdataGDP+1:Tdata)+eps_p;
    
    alpha = [alpha;alpha_pred];
    
    
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
    %             if iPC==1
    %                 xticklabels(MonthWeekJP(xtick1))
    %             else
    %                 xticklabels(MonthWeekJP(xtick1))
    %             end
    legend('Mobility (L axis)','Imputed GDP (R axis)', 'Data GDP (R axis)','FontSize',16,'Location','southeast');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'b';
    ax.YAxis(1).FontSize = fs;
    ax.YAxis(2).FontSize = fs;
    ax.XAxis.FontSize = fs;
    xtickangle(45)
    
    
    %--- Regress mobility on alpha to estimate the elasticity h ---%
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
    
    
    %%%%%%%%%%%%%%%%% Projection starts here %%%%%%%%%%%%%%%%%
    %Normal
    alpha_off = mean(alpha((dateEN >= datetime(2020,10,1)) & (datetime(2020,11,26)>= dateEN ))); % output loss without the state of emergency
    InitialValues = [S(end),I(end),R(end),D(end)];
    
    %--- Construct time series of parameters ---%
    gammaT = gamma*ones(SimPeriod,1);
    delta_sample = delta(end-RetroPeriod+1:end);
    delta_average = sum(delta_sample.*(I(end-RetroPeriod+1:end)/sum(I(end-RetroPeriod+1:end))));
    
    %--- Construct vaccine dstribution ---%
    pace = ps*3600000;
    vacpath = zeros(SimPeriod,1);
    if Vgradual == 1
        % case 1: gradually increasing Vpath
        vacpath(1:10) = (pace/10):(pace/10):pace;
        vacpath(11:end) = pace*ones(42,1);
    else
        % case 2: flat Vpath
        vacpath(5:end) = pace*ones(48,1);
    end
    medical = ps*4700000*0.8;
    elderly = ps*36000000*0.8;
    ordinary = (125710000-36000000-4700000)*ps*0.8;
    elderly_total = ps*36000000;
    [V,deltaT,VT] = vaccine_distribution_medical(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2,ps,POP0,3);
    %             [V,deltaT,VT] = vaccine_distribution_simple(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2);
    delta_ss = delta_average*(0.09/1.28);
    delta_th1 = (delta_average - delta_ss)*0.9+delta_ss;
    delta_ss = (delta_average - delta_ss)*(1-0.8*D2)+delta_ss; %death rate after full vaccination of elderly = 2021/09/30
    beta_r = 0;
    for retrop = retro_lb:retro_ub
        beta_r = beta_r + mean(beta(end-retrop+1:end));
    end
    beta_avg = beta_r/(retro_ub-retro_lb+1);
    betaT = beta_avg*ones(SimPeriod,1);
    
    % ------------------------------------------------ %
    % --------------- Beta loop starts --------------- %
    % ------------------------------------------------ %
    for ivi = 1:nvi
        var_infection = var_infection_vector(ivi);
        for bb = 1:length(beta_option_index)%[beta_option_index]
            close all
            beta_option = beta_option_index(bb);
            if beta_option == 42
                th_off_vec = th_off_vector42{pindex};
                th_off_vec2 = th_off_vector42_2{pindex};
                th_off_vec3 = th_off_vector42_3{pindex};
            elseif beta_option == 41
                th_off_vec = th_off_vector41{pindex};
            else
                th_off_vec = th_off_vector{pindex};
                th_off_vec2 = th_off_vector_2{pindex};
                th_off_vec3 = th_off_vector_3{pindex};
            end
            
            
            for th_off_index = 1:1 %1:length(th_off_vec)
                if th_off_index == 1
                    ERNON_vector42 = ERNON_vector42_1;
                elseif th_off_index == 2
                    ERNON_vector42 = ERNON_vector42_2;
                elseif th_off_index == 3
                    ERNON_vector42 = ERNON_vector42_3;
                end
                
                
                if beta_option == 3
                    betaT_temp = betaT_temp_ini;
                    betaT(1,1) = (1+betaT_temp)*betaT(1,1);
                    for i = 2:length(betaT)
                        betaT_temp = betaT_temp * beta_rho;
                        betaT(i) = betaT(i) * (1+betaT_temp);
                    end
                elseif beta_option == 42
                    var_initial = var_initial_vector(pindex); %scaleB *var_initial_vector(pindex) / v_scale_vec(pindex);
                    var_growth = var_growth_vec(2);
                    logit_initial = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
                    %         var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept
                    var_intercept = logit_initial;
                    
                    % step 1: extrapolation in the variant share in the last 17 weeks     -mean(retro_lb,retro_ub)+1 = -17 + 1
                    var_share = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
                    var_share_prev = exp((-mean(retro_lb,retro_ub)+1:0)'*var_growth+var_intercept)./(1+exp((-mean(retro_lb,retro_ub)+1:0)'*var_growth+var_intercept));
                    % step 2: relative infection rate in the last 17 weeks
                    relative_var_infection_prev = 1 + var_share_prev * var_infection;
                    % step 3:
                    no_var_beta = beta(end-mean(retro_lb,retro_ub)+1:end) ./ relative_var_infection_prev;
                    beta_bar = mean(no_var_beta);
                    
                    if SimCases ~= 4
                        betaT = beta_bar*(1+var_infection*var_share);
                    end
                    
                    betaT2 = betaT;
                    %             betaT2(1:5) =  1.2*betaT2(1:5); % 1.3880 = mean (Dec10 - Jan7, 5 weeks)
                    
                    plot_var_share(1,1) = var_initial;
                    plot_var_share(2:length(SimDate)+1,1) = var_share(:,1);
                    % plot_betaT(1,1) = beta_avg;
                    tt = 12; %Showing previous t periods for the plot
                    
                    plot_betaT(1:tt,1) = beta(end-tt+1:end);
                    plot_betaT(tt+1:length(SimDate)+tt,1) = betaT(:,1);
                    
                    % plot_betaT2(1,1) = beta_avg*(1+betaT_temp);
                    %                 plot_betaT2(1:tt,1) = beta(end-tt+1:end);
                    %                 plot_betaT2(tt+1:length(SimDate)+tt,1) = betaT2(:,1);
                    
                    plot_betaT3(1:tt,1) = beta(end-tt+1:end);
                    plot_betaT3(tt+1:length(SimDate)+tt,1) = beta_avg*ones(length(SimDate),1);
                    
                    xaxis_vec = 0:1:length(SimDate);
                    xaxis_vec2 = 1:1:length(SimDate)+tt;
                    xaxis_vec3 = 1:1:Tdata+tt;
                    
                    MonthWeekEN_Sim = MonthWeekEN(length(dateD)+1:end,1);
                    MonthNumber_Sim = MonthNumber(length(dateD)+1:end,1);
                    WeekNumber_Sim = WeekNumber(length(dateD)+1:end,1);
                    MonthWeekJP_Sim = MonthWeekJP(length(dateD)+1:end,1);
                    
                    MonthWeekEN_Sim2 = MonthWeekEN(length(dateD)-tt+1:end,1);
                    MonthNumber_Sim2 = MonthNumber(length(dateD)-tt+1:end,1);
                    WeekNumber_Sim2 = WeekNumber(length(dateD)-tt+1:end,1);
                    MonthWeekJP_Sim2 = MonthWeekJP(length(dateD)-tt+1:end,1);
                    
                    
                    for l = 1:2
                        if l == 1
                            figure(4201)
                        elseif l == 2
                            figure(4202)
                        end
                        set(gcf,'Position',[100,100,1200,500])
                        subplot(1,2,1)
                        plot(xaxis_vec,plot_var_share(:,1)*100,'-r','LineWidth',2)
                        ax = gca;
                        ax.YAxis.FontSize = 20;
                        ax.XAxis.FontSize = 20;
                        ax.YAxis.Exponent = 0;
                        ytickformat('%,3.0f')
                        xticks(find(WeekNumber_Sim==1))
                        if l == 1
                            % legend('Variant Share','FontSize',10)
                            title('Projected path of variants','FontSize',20,'FontWeight','normal')
                            xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
                            xlabel('Weeks','FontSize',20)
                            ylabel('Share of Variants(%)','FontSize',20)
                        elseif l == 2
                            % legend('変異株割合','FontName',fn,'FontSize',10)
                            title('新規感染者数の推移','FontSize',20,'FontWeight','normal')
                            title('変異株シェアの推移','FontSize',20,'FontWeight','normal','FontName',fn)
                            xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
                            xlabel('週','FontSize',20,'FontName',fn)
                            ylabel('変異株割合(%)','FontSize',20,'FontName',fn)
                        end
                        xtickangle(45)
                        xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
                        xlim([0,length(SimDate)]);
                        % lgd = legend;
                        % lgd.Location = 'northwest';
                        
                        subplot(1,2,2)
                        plot(xaxis_vec2(tt+1:end), plot_betaT(tt+1:end,1),'-r','LineWidth',2)
                        hold on
                        %                     plot(xaxis_vec2(tt+1:end), plot_betaT2(tt+1:end,1),'-r','LineWidth',2.0)
                        plot(xaxis_vec2,plot_betaT3,'-k','LineWidth',2)
                        ax = gca;
                        ax.YAxis.FontSize = 20;
                        ax.XAxis.FontSize = 20;
                        ax.YAxis.Exponent = 0;
                        ytickformat('%,0.2f')
                        xticks(find(WeekNumber_Sim2==1))
                        xtickangle(45)
                        if l == 1
                            legend('Raw Infection Rate with Variants','Past Average','FontSize',10)
                            title('Projected path of beta','FontSize',20,'FontWeight','normal')
                            xticklabels(MonthWeekEN_Sim2(WeekNumber_Sim2==1))
                            xlabel('Weeks','FontSize',20)
                            ylabel('Infection Rate','FontSize',20,'FontName',fn)
                        elseif l == 2
                            legend('感染率(変異株成長ケース)','過去の平均','FontSize',10,'FontName',fn)
                            xticklabels(MonthWeekJP_Sim2(WeekNumber_Sim2==1))
                            title('βの推移','FontSize',20,'FontWeight','normal','FontName',fn)
                            xlabel('週','FontSize',20)
                            ylabel('感染率','FontSize',20,'FontName',fn)
                        end
                        xline(tt+1,'--','LineWidth',1.5,'HandleVisibility','off');
                        xlim([0,length(SimDate)+tt]);
                        lgd = legend;
                        lgd.Location = 'southeast';
                    end
                    if figure_save == 1
                        saveas(figure(4201),[home 'Figures/' char(pref) '/beta_path' sprintf('%.0f', beta_option) '.png']);
                        saveas(figure(4202),[home 'Figures/' char(pref) '/beta_path' sprintf('%.0f', beta_option) '_jp.png']);
                    end
                    
                    %                 betaT = betaT2;
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Projection parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                th_on = th_on_vector(pindex)*7;         % threshold to place the state of emergency
                th_on2 = th_on_vector2(pindex)*7;         % threshold to place the state of emergency
                th_on3 = th_on_vector3(pindex)*7;         % threshold to place the state of emergency
                th_off = th_off_vec(th_off_index)*7;
                th_off2 = th_off_vec2(th_off_index)*7;
                th_off3 = th_off_vec3(th_off_index)*7;
                ERN_on = ERN_on_vector(pindex);
                if beta_option == 42
                    ERN_on = ERN_on_scenario2{pindex};
                end
                
                ERN_now = ERN_now_vector(pindex);
                altA_now = (((ERN_now*(POP0/S(end))*((gammaT(1)+deltaT(1))/beta_avg)).^(1/k))-1)*(h(1)/h(2));
                altA_on = (((ERN_on*(POP0/S(end))*((gammaT(1)+deltaT(1))/beta_avg)).^(1/k))-1)*(h(1)/h(2));
                ERNCheck = (S(end)/POP0).*(((1+(h(2)/h(1))*alpha_off).^k).*beta_avg)./(gammaT(1)+deltaT(1));
                
                y = SimCases % Choose which analysis you want
                if y == 4 && beta_option < 41
                    break
                end
                TL = TL_vector{pindex};
                TL_index = TL_index_vector{pindex};
                if beta_option == 1 || beta_option == 3
                    ON = THON_vector1{pindex}*7;
                    ON_index = THON_index_vector1{pindex};
                elseif beta_option == 42
                    ON = THON_vector42{pindex}*7;
                    ON_index = THON_index_vector42{pindex};
                end
                if beta_option == 1 || beta_option == 3
                    AON = (((ERNON_vector1{pindex}.*(POP0/S(end)).*((gammaT(1)+deltaT(1))./beta_avg)).^(1/k))-1).*(h(1)/h(2));
                    AON_index = ERNON_index_vector1{pindex};
                elseif beta_option == 42
                    AON = (((ERNON_vector42{pindex}.*(POP0/S(end)).*((gammaT(1)+deltaT(1))./beta_avg)).^(1/k))-1).*(h(1)/h(2));
                    AON_index = ERNON_index_vector42{pindex};
                end
                DR = DR_vector{pindex};
                if y == 1
                    TH = ON/7;
                    TH_index = ON_index;
                elseif y == 2
                    TH = DR;
                    TH_index = DR_index;
                elseif y == 3
                    if beta_option == 1 || beta_option == 3
                        TH = ERNON_vector1{pindex};
                    elseif beta_option == 42
                        TH = ERNON_vector42{pindex};
                    end
                    TH_index = AON_index;
                elseif y == 4
                    TH = 1+var_infection_vector;
                    TH_index = 1+var_infection_index;
                elseif y == 5
                    TH = TL;
                    TH_index = TL_index;
                end
                DMat = nan(SimCases,length(TH));
                AlphaMat = nan(SimCases,length(TH));
                SimData = nan(SimPeriod+1,length(InitialValues),length(TH),SimCases);
                AlphaPath = nan(SimPeriod,length(TH),SimCases);
                NPath = nan(SimPeriod,length(TH),SimCases);
                SimERN = nan(SimPeriod,length(TH),SimCases);
                BackDataN = zeros(SimPeriod+8,length(TH_index),SimCases);
                BackDataAlpha = zeros(SimPeriod+8,length(TH_index),SimCases);
                BackDataERN = zeros(SimPeriod+8,length(TH_index),SimCases);
                BackDataDA = zeros(length(TH),3,SimCases);
                
                %---- 1. Different thresholds to initiate the state of emergency ---%
                if y == 1
                    for i = 1:length(ON)
                        if pindex < 5
                            [DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData(:,:,i,1),NPath(:,i,1),SimERN(:,i,1)] = Covid_projection_control_gradual_off_threshold_delta2(InitialValues,altA_on,alpha_off,ON(i),ON(i),th_off,th_off2,betaT,gammaT,deltaT,delta_ss,delta_th1,V,h,k,POP0,hconstant,1+DR_index(1),alpha(end));
                        else
                            [DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData(:,:,i,1),NPath(:,i,1),SimERN(:,i,1)] = Covid_projection_control_gradual_threshold_delta2(InitialValues,altA_on,alpha_off,ON(i),ON(i),th_off,th_off2,betaT,gammaT,deltaT,delta_ss,delta_th1,V,h,k,POP0,hconstant,1+DR_index(1));
                        end
                    end
                    ON = ON/7;
                    
                    %---- 2. Different durations for economic recovery ---%
                elseif y == 2
                    for i = 1:length(DR)
                        if pindex < 5
                            [DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData(:,:,i,2),NPath(:,i,2),SimERN(:,i,2)] = Covid_projection_control_gradual_off_threshold_delta2(InitialValues,altA_on,alpha_off,th_on,th_on2,th_off,th_off2,betaT,gammaT,deltaT,delta_ss,delta_th1,V,h,k,POP0,hconstant,1+DR(i),alpha(end));
                        else
                            [DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData(:,:,i,2),NPath(:,i,2),SimERN(:,i,2)] = Covid_projection_control_gradual_threshold_delta2(InitialValues,altA_on,alpha_off,th_on,th_on2,th_off,th_off2,betaT,gammaT,deltaT,delta_ss,delta_th1,V,h,k,POP0,hconstant,1+DR(i));
                        end
                    end
                    
                    %---- 3. Different ERN on ---%
                elseif y == 3
                    for i = 1:length(AON)
                        if pindex < 5
                            [DMat(3,i),AlphaMat(3,i),AlphaPath(:,i,3),SimData(:,:,i,3),NPath(:,i,3),SimERN(:,i,3)] = Covid_projection_control_gradual_off_threshold_delta2(InitialValues,AON(i),alpha_off,th_on,th_on2,th_off,th_off2,betaT,gammaT,deltaT,delta_ss,delta_th1,V,h,k,POP0,hconstant,1+DR_index(1),alpha(end));
                        else
                            [DMat(3,i),AlphaMat(3,i),AlphaPath(:,i,3),SimData(:,:,i,3),NPath(:,i,3),SimERN(:,i,3)] = Covid_projection_control_gradual_threshold_delta2(InitialValues,AON(i),alpha_off,th_on,th_on2,th_off,th_off2,betaT,gammaT,deltaT,delta_ss,delta_th1,V,h,k,POP0,hconstant,1+DR_index(1));
                        end
                    end
                    
                    %-----------------------------------------------------%
                    
                    
                    %---- 4  Different infection rate (only for variant scenario) ---%
                elseif y == 4
                    for i = 1:length(var_infection_vector)
                        var_infection = var_infection_vector(i);
                        betaT = beta_bar*(1+var_infection*var_share);
                        if pindex < 5
                            [DMat(4,i),AlphaMat(4,i),AlphaPath(:,i,4),SimData(:,:,i,4),NPath(:,i,4),SimERN(:,i,4)] = Covid_projection_control_gradual_off_threshold_delta2(InitialValues,altA_on,alpha_off,th_on,th_on2,th_off,th_off2,betaT,gammaT,deltaT,delta_ss,delta_th1,V,h,k,POP0,hconstant,1+DR_index(1),alpha(end));
                        else
                            [DMat(4,i),AlphaMat(4,i),AlphaPath(:,i,4),SimData(:,:,i,4),NPath(:,i,4),SimERN(:,i,4)] = Covid_projection_control_gradual_threshold_delta2(InitialValues,altA_on,alpha_off,th_on,th_on2, th_off,th_off2,betaT,gammaT,deltaT,delta_ss,delta_th1,V,h,k,POP0,hconstant,1+DR_index(1));
                        end
                    end
                    %---- 5  Different Thresholds to Lift  ---%
                elseif y ==5
                    for i = 1:length(TL)
                        if pindex < 5
                            [DMat(5,i),AlphaMat(5,i),AlphaPath(:,i,5),SimData(:,:,i,5),NPath(:,i,5),SimERN(:,i,5)] = Covid_projection_control_gradual_off_threshold_delta2(InitialValues,altA_on,alpha_off,th_on,th_on2,TL(i)*7,TL(i)*7,betaT,gammaT,deltaT,delta_ss,delta_th1,V,h,k,POP0,hconstant,1+DR_index(1),alpha(end));
                        else
                            [DMat(5,i),AlphaMat(5,i),AlphaPath(:,i,5),SimData(:,:,i,5),NPath(:,i,5),SimERN(:,i,5)] = Covid_projection_control_gradual_threshold_delta2(InitialValues,altA_on,alpha_off,th_on,th_on2,TL(i)*7,TL(i)*7,betaT,gammaT,deltaT,delta_ss,delta_th1,V,h,k,POP0,hconstant,1+DR_index(1));
                        end
                    end
                end
                
                %             minAlpha = min(AlphaMat,[],'all');
                
                if y == 1
                    minAlpha1 = min(AlphaMat(y,:));
                    if beta_option == 1 && ivi == 1
                        save('alpha_min1.mat', 'minAlpha1')
                    end
                    minAlphaMat1 = load('alpha_min1.mat','minAlpha1');
                    minAlpha = minAlphaMat1.minAlpha1;
                elseif y == 2
                    minAlpha2 = min(AlphaMat(y,:));
                    if beta_option == 1 && ivi == 1
                        save('alpha_min2.mat', 'minAlpha2')
                    end
                    minAlphaMat2 = load('alpha_min2.mat','minAlpha2');
                    minAlpha = minAlphaMat2.minAlpha2;
                elseif y == 3
                    minAlpha3 = min(AlphaMat(y,:));
                    if beta_option == 1 && ivi == 1
                        save('alpha_min3.mat', 'minAlpha3')
                    end
                    minAlphaMat3 = load('alpha_min3.mat','minAlpha3');
                    minAlpha = minAlphaMat3.minAlpha3;
                elseif y == 4
                    minAlpha4 = min(AlphaMat(y,:));
                    save('alpha_min4.mat', 'minAlpha4')
                    minAlphaMat4 = load('alpha_min4.mat','minAlpha4');
                    minAlpha = minAlphaMat4.minAlpha4;
                elseif y == 5
                    minAlpha5 = min(AlphaMat(y,:));
                    if beta_option == 1 && ivi == 1
                        save('alpha_min5.mat', 'minAlpha5')
                    end
                    minAlphaMat5 = load('alpha_min5.mat','minAlpha5');
                    minAlpha = minAlphaMat5.minAlpha5;
                end
                
                TH_th(1:length(TH),y,th_off_index) = TH';
                AlphaM = AlphaMat(y,:);
                AlphaM_th(1:length(TH),th_off_index) = AlphaM;
                AlphaM = AlphaM(~isnan(AlphaM));
                DM = DMat(y,:);
                DM_th(1:length(TH),th_off_index) = DM;
                DM = DM(~isnan(DM));
                AlphaM = (AlphaM - minAlpha)*prefGDP*10000;
                AlphaM_th(1:length(TH),th_off_index) = (AlphaM_th(1:length(TH),th_off_index) - minAlpha)*prefGDP*10000;
                
                BackDataDA(1:length(TH),:,y) = [round(AlphaM'),round(DM'),TH'];
                %--- Record how many times on and off are triggered ---%
                waves = zeros(1,length(TH));
                for i = 1:length(TH)
                    svec = zeros(SimPeriod-1,1);
                    for t = 1:SimPeriod-1
                        svec(t) = AlphaPath(t+1,i,y)-AlphaPath(t,i,y);
                    end
                    waves(i) = sum(svec>0);
                end
                waves_th(1:length(waves),y,th_off_index) = waves;
                
                if y == 1
                    ATOFF = AlphaM(waves==0);
                    DTOFF = DM(waves==0);
                    % THTL = TH(waves==0);
                    THON = TH(waves==0);
                    ATONblue = AlphaM(abs(TH - TH_index(2)) < 0.0001);
                    DTONblue = DM(abs(TH - TH_index(2)) < 0.0001);
                elseif y == 2
                    ADR = AlphaM(waves==0);
                    DDR = DM(waves==0);
                    THON = TH(waves==0);
                    ADRblue = AlphaM(abs(TH - TH_index(2)) < 0.0001);
                    DDRblue = DM(abs(TH - TH_index(2)) < 0.0001);
                elseif y == 3
                    AERN = AlphaM(waves==0);
                    DERN = DM(waves==0);
                    % THTL = TH(waves==0);
                    THON = TH(waves==0);
                    AERNblue = AlphaM(abs(TH - TH_index(2)) < 0.0001);
                    DERNblue = DM(abs(TH - TH_index(2)) < 0.0001);
                elseif y == 4
                    Ainfec = AlphaM(waves==0);
                    Dinfec = DM(waves==0);
                    % THTL = TH(waves==0);
                    THON = TH(waves==0);
                    Ainfecblue = AlphaM(abs(TH - TH_index(2)) < 0.0001);
                    Dinfecblue = DM(abs(TH - TH_index(2)) < 0.0001);
                elseif y == 5
                    ATL = AlphaM(waves==0);
                    DTL = DM(waves==0);
                    % THTL = TH(waves==0);
                    THON = TH(waves==0);
                    ATLblue = AlphaM(abs(TH - TH_index(2)) < 0.0001);
                    DTLblue = DM(abs(TH - TH_index(2)) < 0.0001);
                end
                
                for l = 1:2 %1:2 when english version needed
                    if l == 1
                        if beta_option == 41
                            figure((41*1000+11+y)*100)
                        elseif beta_option == 42
                            figure((42*1000+11+y)*100)
                        else
                            figure((11+y)*100)
                        end
                    elseif l == 2
                        if beta_option == 41
                            figure((41*1000+(11+y)*10+1)*100)
                        elseif beta_option == 42
                            figure((42*1000+(11+y)*10+1)*100)
                        else
                            figure(((11+y)*10+1)*100)
                        end
                    end
                    set(gcf,'Position',[100,100,1000,800])
                    subplot(3,2,1)
                    if max(TH)<3
                        ft = '%.2f';
                    else
                        ft = '%.0f';
                    end
                    %             if y == 1 || y == 2 || y ==5
                    for i = 1:length(TH)
                        if abs(TH(i) - TH_index(1)) < 0.0001
                            plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
                            BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                            BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                            BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                        elseif abs(TH(i) - TH_index(2)) < 0.0001
                            plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
                            BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                            BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                            BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                        else
                            plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf(ft,TH(i)))
                        end
                        hold on
                    end
                    %             elseif y == 3 || y == 4
                    %                 for i = 1:length(TH)
                    %                     if abs(TH(i) - TH_index(1)) < 0.0001
                    %                         plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.2f',TH(i)))
                    %                         BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                    %                         BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                    %                         BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                    %                     elseif abs(TH(i) - TH_index(2)) < 0.0001
                    %                         plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.2f',TH(i)))
                    %                         BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                    %                         BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                    %                         BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                    %                     else
                    %                         plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.2f',TH(i)))
                    %                     end
                    %                     hold on
                    %                 end
                    %             end
                    plot(N/7,'k','LineWidth',2,'HandleVisibility','off')
                    xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
                    ax = gca;
                    ax.YAxis.FontSize = 10;
                    ax.XAxis.FontSize = 10;
                    ax.YAxis.Exponent = 0;
                    ytickformat('%,6.0f')
                    xticks(find(WeekNumber==1))
                    if l == 1
                        title('Projected path of new cases','FontSize',10,'FontWeight','normal')
                        xticklabels(MonthWeekEN(WeekNumber==1))
                        %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
                    elseif l == 2
                        %title('新規感染者数の推移','FontSize',20,'FontWeight','normal')
                        title('新規感染者数の推移','FontSize',10,'FontWeight','normal','FontName',fn)
                        xticklabels(MonthWeekJP(WeekNumber==1))
                    end
                    lgd = legend;
                    lgd.NumColumns = 2;
                    xtickangle(45)
                    if beta_option == 41 || beta_option == 42
                        xlim([Tdata-7 Tdata+41])
                    else
                        xlim([Tdata-7 Tdata+28])
                    end
                    
                    %--- Number of cumulative deaths ---%
                    subplot(3,2,2)
                    plot(AlphaM(waves==0),DM(waves==0),'-go','LineWidth',2,'MarkerSize',10);
                    hold on
                    plot(AlphaM(waves==1),DM(waves==1),'-mo','LineWidth',2,'MarkerSize',10);
                    hold on
                    plot(AlphaM(waves==2),DM(waves==2),'-bo','LineWidth',2,'MarkerSize',10);
                    hold on
                    text(AlphaM,DM,string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
                    hold on
                    scatter(AlphaM(abs(TH - TH_index(1)) < 0.0001),DM(abs(TH - TH_index(1)) < 0.0001),250,'red','filled');
                    hold on
                    scatter(AlphaM(abs(TH - TH_index(2)) < 0.0001),DM(abs(TH - TH_index(2)) < 0.0001),250,'blue','filled');
                    
                    if l == 1
                        xlabel('Output Loss (hundred million yen)','FontSize',10)
                        ylabel('Cumulative Deaths','FontSize',10)
                        title('Relationship between Covid-19 and output','FontSize',10,'FontWeight','normal')
                    elseif l == 2
                        xlabel('経済損失 (億円)','FontSize',10,'FontName',fn)
                        ylabel('累計死亡者数','FontSize',10,'FontName',fn)
                        title('コロナ感染と経済の関係','FontSize',10,'FontWeight','normal','FontName',fn)
                    end
                    xlim([0,inf])
                    xtickangle(45)
                    grid on
                    ax = gca;
                    ax.YAxis.FontSize = 10;
                    ax.XAxis.FontSize = 10;
                    ax.YAxis.Exponent = 0;
                    ax.XAxis.Exponent = 0;
                    ytickformat('%,6.0f')
                    
                    % alpha path
                    subplot(3,2,3)
                    for i = 1:length(TH)
                        if abs(TH(i) - TH_index(1)) < 0.0001
                            plot([1-alpha;1-AlphaPath(:,i,y)]*100,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                        elseif abs(TH(i) - TH_index(2)) < 0.0001
                            plot([1-alpha;1-AlphaPath(:,i,y)]*100,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                        else
                            plot([1-alpha;1-AlphaPath(:,i,y)]*100,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                        end
                        hold on
                    end
                    grid on
                    xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
                    ax = gca;
                    ax.YAxis.FontSize = 10;
                    ax.XAxis.FontSize = 10;
                    xticks(find(WeekNumber==1))
                    xtickangle(45)
                    ytickformat('%,3.0f')
                    xticklabels(MonthWeekJP(WeekNumber==1))
                    xlim([0 inf])
                    ylim([80 100])
                    if l == 1
                        xlabel('Weeks)','FontSize',10)
                        ylabel('1-alpha','FontSize',10)
                        title('1-alpha Path','FontSize',10,'FontWeight','normal')
                    elseif l == 2
                        xlabel('Weeks','FontSize',10,'FontName',fn)
                        ylabel('alpha','FontSize',10,'FontName',fn)
                        title('1-alpha path','FontSize',10,'FontWeight','normal','FontName',fn)
                    end
                    
                    % ERN path
                    subplot(3,2,4)
                    for i = 1:length(TH)
                        if abs(TH(i) - TH_index(1)) < 0.0001
                            plot([ERN;SimERN(:,i,y)],'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                        elseif abs(TH(i) - TH_index(2)) < 0.0001
                            plot([ERN;SimERN(:,i,y)],'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                        else
                            plot([ERN;SimERN(:,i,y)],'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                        end
                        hold on
                    end
                    grid on
                    xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
                    xticks(find(WeekNumber==1))
                    xtickangle(45)
                    ytickformat('%,2.4f')
                    xticklabels(MonthWeekJP(WeekNumber==1))
                    xlim([Tdata-12 inf])
                    ax = gca;
                    ax.YAxis.FontSize = 10;
                    ax.XAxis.FontSize = 10;
                    if l == 1
                        xlabel('Weeks)','FontSize',10)
                        ylabel('ERN','FontSize',10)
                        title('ERN Path','FontSize',10,'FontWeight','normal')
                    elseif l == 2
                        xlabel('Weeks','FontSize',10,'FontName',fn)
                        ylabel('ERN','FontSize',10,'FontName',fn)
                        title('ERN path','FontSize',10,'FontWeight','normal','FontName',fn)
                    end
                    
                    %GDP loss
                    subplot(3,2,5)
                    for i = 1:length(TH)
                        if abs(TH(i) - TH_index(1)) < 0.0001
                            plot([(1-alpha).*referenceGDP(1:Tdata)./100;(1-AlphaPath(:,i,y)).*referenceGDP(Tdata+1:Tdata+length(SimDate))./100]*100,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                        elseif abs(TH(i) - TH_index(2)) < 0.0001
                            plot([(1-alpha).*referenceGDP(1:Tdata)./100;(1-AlphaPath(:,i,y)).*referenceGDP(Tdata+1:Tdata+length(SimDate))./100]*100,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                        else
                            plot([(1-alpha).*referenceGDP(1:Tdata)./100;(1-AlphaPath(:,i,y)).*referenceGDP(Tdata+1:Tdata+length(SimDate))./100]*100,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                        end
                        hold on
                    end
                    grid on
                    xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
                    ax = gca;
                    ax.YAxis.FontSize = 10;
                    ax.XAxis.FontSize = 10;
                    xticks(find(WeekNumber==1))
                    xtickangle(45)
                    ytickformat('%,3.0f')
                    xticklabels(MonthWeekJP(WeekNumber==1))
                    xlim([0 inf])
                    ylim([80 100])
                    if l == 1
                        xlabel('Weeks','FontSize',10)
                        ylabel('GDP','FontSize',10)
                        title('GDP Path','FontSize',10,'FontWeight','normal')
                    elseif l == 2
                        xlabel('Weeks','FontSize',10,'FontName',fn)
                        ylabel('GDP','FontSize',10,'FontName',fn)
                        title('GDP path','FontSize',10,'FontWeight','normal','FontName',fn)
                    end
                    
                    
                    % Generate graphs for the website
                    if l == 1
                        if beta_option == 41
                            figure(41*1000+11+y)
                        elseif beta_option == 42
                            figure(42*1000+11+y + ivi*10000)
                        else
                            figure(11+y)
                        end
                    elseif l == 2
                        if beta_option == 41
                            figure(41*1000+(11+y)*10+1)
                        elseif beta_option == 42
                            figure(42*1000+(11+y)*10+1 + ivi*10000)
                        else
                            figure((11+y)*10+1)
                        end
                    end
                    set(gcf,'Position',[100,100,1200,500])
                    
                    subplot(1,2,1)
                    if max(TH)<3
                        ft = '%.2f';
                    else
                        ft = '%.0f';
                    end
                    for i = 1:length(TH)
                        if abs(TH(i) - TH_index(1)) < 0.0001
                            plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
                            BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                            BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                            BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                        elseif abs(TH(i) - TH_index(2)) < 0.0001
                            plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
                            BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                            BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                            BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                        elseif abs(TH(i) - TH_index(3)) < 0.0001
                            plot([N;NPath(:,i,y)]/7,'-k','LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
                            BackDataN(:,3,y) = [N(end-7:end);NPath(:,i,y)];
                            BackDataAlpha(:,3,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                            BackDataERN(:,3,y) = [ERN(end-7:end);SimERN(:,i,y)];
                        else
                            plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf(ft,TH(i)))
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
                    % uncomment if you are doing the ERN comparison.
%                     if y == 2 || y == 3 || y == 4
%                         %yline(th_on/7,'-.r',sprintf('%.0f',th_on/7),'LineWidth',0.5,'HandleVisibility','off');
%                         if th_on ~= th_on2 && max(waves) > 1
%                             yline(th_on2/7,'--r',sprintf('%.0f',th_on2/7),'LineWidth',0.5,'HandleVisibility','off');
%                         end
%                         yline(th_off/7,'-.b',sprintf('%.0f',th_off/7),'LineWidth',0.5,'HandleVisibility','off');
%                         %                 yticks(th_off/7);
%                         %                 ytickformat('%,6.0f')
%                     end
                    xticks(find(WeekNumber==1))
                    if l == 1
                        title('Projected path of new cases','FontSize',20,'FontWeight','normal')
                        xticklabels(MonthWeekEN(WeekNumber==1))
                        %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
                    elseif l == 2
                        title('新規感染者数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
                        xticklabels(MonthWeekJP(WeekNumber==1))
                    end
                    if y ~= 3 || Vsimple == 0
                        lgd = legend;
                    end
                    lgd.NumColumns = 2;
                    if pindex == 1 || pindex == 5
                        lgd.Location = 'northwest';
                    else
                        lgd.Location = 'northeast';
                    end
                    xtickangle(45)
                    if beta_option == 41 || beta_option == 42
                        xlim([Tdata-7 Tdata+36])
                    else
                        xlim([Tdata-7 Tdata+28])
                    end
                    
                    %--- Number of cumulative deaths ---%
                    subplot(1,2,2)
                    plot(AlphaM(waves==0),DM(waves==0),'-go','LineWidth',2,'MarkerSize',10);
                    hold on
                    plot(AlphaM(waves==1),DM(waves==1),'-mo','LineWidth',2,'MarkerSize',10);
                    hold on
                    plot(AlphaM(waves==2),DM(waves==2),'-bo','LineWidth',2,'MarkerSize',10);
                    hold on
                    plot(AlphaM(waves==3),DM(waves==3),'-ko','LineWidth',2,'MarkerSize',10);
                    %                 if y == 4
                    %                     text(AlphaM,DM,string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
                    %                 elseif y ~= 3 || Vsimple == 0
                    text(AlphaM,DM,string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
                    %                 end
                    hold on
                    scatter(AlphaM(abs(TH - TH_index(1)) < 0.0001),DM(abs(TH - TH_index(1)) < 0.0001),250,'r','filled');
                    hold on
                    scatter(AlphaM(abs(TH - TH_index(2)) < 0.0001),DM(abs(TH - TH_index(2)) < 0.0001),250,'b','filled');
                    if length(TH_index) >=3
                        scatter(AlphaM(abs(TH - TH_index(3)) < 0.0001),DM(abs(TH - TH_index(3)) < 0.0001),250,'k','filled');
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
                    if beta_option ~= 3 || betaT_temp_ini>=0
                        xlim([0,inf])
                    end
                    xtickangle(45)
                    grid on
                    ax = gca;
                    ax.YAxis.FontSize = 20;
                    ax.XAxis.FontSize = 20;
                    ax.YAxis.Exponent = 0;
                    ax.XAxis.Exponent = 0;
                    ytickformat('%,6.0f')
                end %End of language loop = figure loop
                
                
                if data_save == 1
                    titleN = strings(1,1+length(TH_index)*3);
                    titleN(1) = "週";
                    for ti = 1:length(TH_index)
                        titleN(1,1+ti) = string(['新規感染者数（',sprintf('%.0f',TH_index(ti)),'）']);
                        titleN(1,1+length(TH_index)+ti) = string(['経済活動（',sprintf('%.0f',TH_index(ti)),'）']);
                        titleN(1,1+length(TH_index)*2+ti) = string(['実効再生産数（',sprintf('%.0f',TH_index(ti)),'）']);
                    end
                    TN = table([titleN;MonthWeekJP(Tdata-7:end-1),round(BackDataN(:,1:length(TH_index),y)/7),round(100*(1-BackDataAlpha(:,1:length(TH_index),y)),1),round(BackDataERN(:,1:length(TH_index),y),2)]);
                    titleAD = ["経済損失（億円）","死亡者数","ケース"];
                    TAD = table([titleAD;BackDataDA(1:length(TH),:,y)]);
                    if beta_option == 42
                        if y == 1
                            %                 writetable(TN,[home 'Figures/' char(pref) '/BackData_Thresholds' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            %                 writetable(TAD,[home 'Figures/' char(pref) '/BackData_Thresholds' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                            writetable(TN,[home 'Figures/' char(pref) '/BackData_ThresholdsON' char(pref) '_' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            writetable(TAD,[home 'Figures/' char(pref) '/BackData_ThresholdsON' char(pref) '_' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi))  '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                        elseif y == 2
                            writetable(TN,[home 'Figures/' char(pref) '/BackData_GradualRecovery' char(pref) '_' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi))  '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            writetable(TAD,[home 'Figures/' char(pref) '/BackData_GradualRecovery' char(pref) '_' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi))  '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                        elseif y == 3
                            %                 writetable(TN,[home 'Figures/' char(pref) '/BackData_VaccineVariation' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            %                 writetable(TAD,[home 'Figures/' char(pref) '/BackData_VaccineVariation' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                            writetable(TN,[home 'Figures/' char(pref) '/BackData_EmergencyState' char(pref) '_' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi))  '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            writetable(TAD,[home 'Figures/' char(pref) '/BackData_EmergencyState' char(pref) '_' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi))  '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                        elseif y == 4
                            writetable(TN,[home 'Figures/' char(pref) '/BackData_VarInfection' char(pref) '_' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi))  '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            writetable(TAD,[home 'Figures/' char(pref) '/BackData_VarInfection' char(pref) '_' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi))  '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                        elseif y == 5
                            writetable(TN,[home 'Figures/' char(pref) '/BackData_Thresholds' char(pref) '_' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi))  '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            writetable(TAD,[home 'Figures/' char(pref) '/BackData_Thresholds' char(pref) '_' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi))  '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                        end
                    else
                        if y == 1
                            %                 writetable(TN,[home 'Figures/' char(pref) '/BackData_Thresholds' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            %                 writetable(TAD,[home 'Figures/' char(pref) '/BackData_Thresholds' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                            writetable(TN,[home 'Figures/' char(pref) '/BackData_ThresholdsON' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            writetable(TAD,[home 'Figures/' char(pref) '/BackData_ThresholdsON' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                        elseif y == 2
                            writetable(TN,[home 'Figures/' char(pref) '/BackData_GradualRecovery' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            writetable(TAD,[home 'Figures/' char(pref) '/BackData_GradualRecovery' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                        elseif y == 3
                            %                 writetable(TN,[home 'Figures/' char(pref) '/BackData_VaccineVariation' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            %                 writetable(TAD,[home 'Figures/' char(pref) '/BackData_VaccineVariation' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                            writetable(TN,[home 'Figures/' char(pref) '/BackData_EmergencyState' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            writetable(TAD,[home 'Figures/' char(pref) '/BackData_EmergencyState' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                        elseif y == 4
                            writetable(TN,[home 'Figures/' char(pref) '/BackData_VarInfection' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            writetable(TAD,[home 'Figures/' char(pref) '/BackData_VarInfection' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                        elseif y == 5
                            writetable(TN,[home 'Figures/' char(pref) '/BackData_Thresholds' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                            writetable(TAD,[home 'Figures/' char(pref) '/BackData_Thresholds' char(pref) '_' sprintf('%.0f', beta_option) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
                        end
                    end
                end
                
                for i = 1:length(TH)
                    NPath_beta(:,i,y,bb) = NPath(:,i,y);
                    AlphaM_beta(i,y,bb) = AlphaM(i);
                    DM_beta(i,y,bb) = DM(i);
                    waves_beta(i,y,bb) = waves(i);
                end
                
                if beta_option == 42
                    plot_ERN(1:tt,1) = ERN(end-tt+1:end);
                    plot_ERN(tt+1:length(SimDate)+tt,1) = SimERN(:,(abs(TH - TH_index(1)) < 0.0001),y);
                    for l = 1:2
                        if l == 1
                            figure(5101)
                        elseif l == 2
                            figure(5102)
                        end
                        set(gcf,'Position',[100,100,1200,500])
                        
                        subplot(1,2,1)
                        plot(xaxis_vec,plot_var_share(:,1)*100,'-r','LineWidth',2)
                        ax = gca;
                        ax.YAxis.FontSize = 20;
                        ax.XAxis.FontSize = 20;
                        ax.YAxis.Exponent = 0;
                        ytickformat('%,3.0f')
                        xticks(find(WeekNumber_Sim==1))
                        if l == 1
                            % legend('Variant Share','FontSize',10)
                            title('Projected path of variants','FontSize',20,'FontWeight','normal')
                            xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
                            xlabel('Weeks','FontSize',20)
                            ylabel('Share of Variants(%)','FontSize',20)
                        elseif l == 2
                            % legend('変異株割合','FontName',fn,'FontSize',10)
                            title('新規感染者数の推移','FontSize',20,'FontWeight','normal')
                            title('変異株シェアの推移','FontSize',20,'FontWeight','normal','FontName',fn)
                            xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
                            xlabel('週','FontSize',20,'FontName',fn)
                            ylabel('変異株割合(%)','FontSize',20,'FontName',fn)
                        end
                        xtickangle(45)
                        xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
                        xlim([0,length(SimDate)]);
                        ylim([0, 100])
                        
                        subplot(1,2,2)
                        if max(TH)<3
                            ft = '%.2f';
                        else
                            ft = '%.0f';
                        end
                        for i = 1:length(TH)
                            if abs(TH(i) - TH_index(1))<0.001
                                plot([ERN;SimERN(:,i,y)],'-r','LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
                            elseif abs(TH(i) - TH_index(2))<0.001
                                plot([ERN;SimERN(:,i,y)],'-b','LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
                            else
                                plot([ERN;SimERN(:,i,y)],'--','LineWidth',0.3,'DisplayName',sprintf(ft,TH(i)))
                            end
                            hold on
                        end
                        if l == 1
                            plot(xaxis_vec3(Tdata-tt+1:Tdata+1),plot_ERN(1:tt+1),'-k','LineWidth',2,'DisplayName','Actual')
                        elseif l == 2
                            plot(xaxis_vec3(Tdata-tt+1:Tdata+1),plot_ERN(1:tt+1),'-k','LineWidth',2,'DisplayName','これまで')
                        end
                        grid on
                        xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
                        yline(1.0,'LineWidth',1.5,'HandleVisibility','off');
                        xticks(find(WeekNumber==1))
                        xtickangle(45)
                        ytickformat('%,2.4f')
                        % xticklabels(MonthWeekJP(WeekNumber==1))
                        xlim([Tdata-11 inf])
                        ax = gca;
                        ax.YAxis.FontSize = 10;
                        ax.XAxis.FontSize = 10;
                        if l == 1
                            xticklabels(MonthWeekEN(WeekNumber==1))
                            xlabel('Weeks','FontSize',10)
                            ylabel('ERN','FontSize',10)
                            title('ERN Path','FontSize',10,'FontWeight','normal')
                        elseif l == 2
                            xticklabels(MonthWeekJP(WeekNumber==1))
                            title('実効再生産数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
                            xlabel('週','FontSize',20)
                            ylabel('ERN','FontSize',20,'FontName',fn)
                        end
                        legend
                    end
                end
                
                %%
                
                %--- Save figures ---%
                if figure_save == 1
                    %                 saveas(figure(5101),[home 'Figures/' char(pref) '/ERN_path' sprintf('%.0f', beta_option) '.png']);
                    %                 saveas(figure(5102),[home 'Figures/' char(pref) '/ERN_path' sprintf('%.0f', beta_option) '_jp.png']);
                    saveas(figure(2),[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
                    if beta_option == 41
                        if y == 1
                            saveas(figure(41012),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '.png']);
                            saveas(figure(41121),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '_jp.png']);
                        elseif y == 2
                            saveas(figure(41013),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '.png']);
                            saveas(figure(41131),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_jp.png']);
                        elseif y == 3
                            saveas(figure(41014),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '.png']);
                            saveas(figure(41141),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '_jp.png']);
                        elseif y == 5
                            saveas(figure(41016),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '.png']);
                            saveas(figure(41161),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_jp.png']);
                        end
                    elseif beta_option == 42
                        ff1 = 42*1000+11+y + ivi*10000;
                        ff2 = 42*1000+(11+y)*10+1 + ivi*10000;
                        if y ==1
                            saveas(figure(ff1),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '.png']);
                            saveas(figure(ff2),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '_jp.png']);
                        elseif y == 2
                            saveas(figure(ff1),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '.png']);
                            saveas(figure(ff2),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '_jp.png']);
                        elseif y == 3
                            saveas(figure(ff1),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '.png']);
                            saveas(figure(ff2),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '_jp.png']);
                        elseif y ==4
                            saveas(figure(ff1),[home 'Figures/' char(pref) '/VarInfection' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '.png']);
                            saveas(figure(ff2),[home 'Figures/' char(pref) '/VarInfection' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '_jp.png']);
                        elseif y == 5
                            saveas(figure(ff1),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '.png']);
                            saveas(figure(ff2),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_' 'ivi_' sprintf('%.1f', var_infection_vector(ivi)) '_jp.png']);
                        end
                    else
                        if y == 1
                            saveas(figure(12),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option)  '.png']);
                            saveas(figure(121),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '_jp.png']);
                        elseif y == 2
                            saveas(figure(13),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '.png']);
                            saveas(figure(131),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_jp.png']);
                        elseif y == 3
                            saveas(figure(14),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '.png']);
                            saveas(figure(141),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '_jp.png']);
                        elseif y == 5
                            saveas(figure(16),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '.png']);
                            saveas(figure(161),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_jp.png']);
                            %         saveas(figure(14),[home 'Figures/' char(pref) '/VaccineVariation' sprintf('%.0f', beta_option) '.png']);
                            %         saveas(figure(141),[home 'Figures/' char(pref) '/VaccineVariation' sprintf('%.0f', beta_option) '_jp.png']);
                            %         saveas(figure(200),[home 'Figures/' char(pref) '/ThresholdsAndGradual' sprintf('%.0f', beta_option) '_jp.png']);
                        end
                    end
                end
                
                
                
            end %end of th_off vector
            
        end %end of beta_option loop
        
    end %end of nvi(infection_vector) loop
    
end %end of prefecture loop

%% Figures added on 2021/03/25
% cumDataEffV = effectiveness*cumsum(V2_w);
cumDataEffV = E1*cumsum(V1_w) + (E2-E1)*cumsum(V2_w);
cumSimEffV = cumsum(V)+cumDataEffV(end);
cumEffV  = [cumDataEffV; cumSimEffV];

long_delta = [delta;deltaT];
long_beta = [beta;betaT];
long_gamma = [ones(Tdata,1)*gamma;gammaT];
for i=1:length(TH)
    long_alpha(:,i) = [alpha;AlphaPath(:,i,SimCases)];
    long_beta_tilde(:,i) = [beta_tilde;betaT.*(1+h(2)/h(1).*AlphaPath(:,i,SimCases)).^k];
    long_ERN(:,i) = [ERN;SimERN(:,i,SimCases)];
    long_GDP(:,i) = [GDP;(1-AlphaPath(:,i,SimCases)).*referenceGDP(Tdata+1:Tdata+SimPeriod)];
    long_N(:,i) = [N;NPath(:,i,SimCases)];
end
long_S = [S(1:end);SimData(2:end,1,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
long_I = [I(1:end);SimData(2:end,2,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
long_gamI = [gamma*I(2:end);gammaT.*SimData(2:end,2,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
long_R = [R(1:end);SimData(2:end,3,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
long_D = [D(1:end);SimData(2:end,4,(abs(TH - TH_index(1)) < 0.0001),SimCases)];
long_dD = [D(2:end)-D(1:end-1);SimData(2:end,4,(abs(TH - TH_index(1)) < 0.0001),SimCases)-SimData(1:end-1,4,(abs(TH - TH_index(1)) < 0.0001),SimCases)];

figure(1111)
set(gcf,'Position',[100,100,1000,800])
subplot(3,3,1)
plot(long_S)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('Suceptibles')
subplot(3,3,2)
plot(long_N(:,(abs(TH - TH_index(1)) < 0.0001)))
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('New Cases')
subplot(3,3,3)
plot(long_I)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('Infected')
subplot(3,3,4)
plot(long_gamI)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('New Recovery')
subplot(3,3,5)
plot(long_dD)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('New Death')
subplot(3,3,6)
plot(long_R)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('Cumulative Recovery')
subplot(3,3,7)
plot(long_D)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('Cumulative Death')
subplot(3,3,8)
plot(cumEffV)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('Cumulative Effective Vaccinated')

figure(1112)
set(gcf,'Position',[100,100,1000,800])
subplot(3,3,1)
plot(long_delta)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('Delta')

subplot(3,3,2)
plot(long_beta)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('Beta')

subplot(3,3,3)
for i=1:length(TH)
    plot(long_alpha(:,i))
    hold on
end
hold off
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('Alpha')

subplot(3,3,4)
for i = 1:length(TH)
    plot(long_GDP(:,i))
    hold on
end
hold off
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('GDP to Reference GDP')

subplot(3,3,5)
for i = 1:length(TH)
    plot(long_ERN(:,i))
    hold on
end
hold off
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('ERN')

subplot(3,3,6)
for i = 1:length(TH)
    plot(long_N(:,i))
    hold on
end
hold off
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('New Cases')

subplot(3,3,7)
plot(long_delta(30:end))
xline(Tdata-30+1,'LineWidth',1.5,'HandleVisibility','off');
title('Delta since 31st week')

subplot(3,3,8)
plot(long_beta(30:end))
xline(Tdata-30+1,'LineWidth',1.5,'HandleVisibility','off');
title('Beta since 31st week')

subplot(3,3,9)
for i = 1:length(TH)
    plot(long_ERN(30:end,i))
    hold on
end
hold off
xline(Tdata-30+1,'LineWidth',1.5,'HandleVisibility','off');
title('ERN since 31st week')

figure(1114)
Tbound = 30;
set(gcf,'Position',[100,100,1200,800])
subplot(4,2,1)
plot(long_delta(1:Tbound))
title('Delta (1-30 week)')

subplot(4,2,2)
plot(long_beta(1:Tbound))
hold on
plot(long_beta_tilde(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)))
plot(long_N(1:Tbound,(abs(TH - TH_index(1)) < 0.0001))./long_I(1:Tbound))
plot(POP0./long_S(1:Tbound,1))
hold off
legend('beta','beta tilde','N/I','S0/St')
lgd = legend;
lgd.Location = 'southeast';
title('Beta (1-30 week)')

subplot(4,2,3)
yyaxis left
plot(long_N(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)))
hold on
plot(long_I(1:Tbound))
ylabel('# of People')
hold on
yyaxis right
ylabel('New Death')
plot(long_dD(1:Tbound))
legend('N','I','D')
title('New Cases / Infected / Cumlative Death (1-30 week)')

subplot(4,2,4)
yyaxis left
plot(long_alpha(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)))
hold on
ylabel('Alpha')
yyaxis right
plot(long_ERN(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)))
ylabel('ERN')
legend('alpha','ERN')
title('Alpha (1-30 week)')

subplot(4,2,5)
plot(long_delta(Tbound+1:end))
xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
title('Delta since 31st week')

subplot(4,2,6)
plot(long_beta(Tbound+1:end))
hold on
plot(long_beta_tilde(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)))
plot(long_N(Tbound+1:Tdata,(abs(TH - TH_index(1)) < 0.0001))./long_I(Tbound+1:Tdata))
plot(POP0./long_S(Tbound+1:Tdata,1))
hold off
legend('beta','beta tilde','N/I','S0/St')
lgd = legend;
lgd.Location = 'southeast';
xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
ylim([0.5 inf])
title('Beta since 31st week')

subplot(4,2,7)
yyaxis left
plot(long_N(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)))
hold on
plot(long_I(Tbound+1:end))
ylabel('# of People')
hold on
yyaxis right
ylabel('New Death')
plot(long_dD(Tbound+1:end))
legend('N','I','dD')
xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
title('New Cases / Infected / New Death since 31st week')

subplot(4,2,8)
yyaxis left
plot(long_alpha(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)))
hold on
ylabel('Alpha')
yyaxis right
plot(long_ERN(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)))
ylabel('ERN')
legend('alpha','ERN')
xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
title('Alpha since 31st week')


figure(1115)
Tbound = 30;
set(gcf,'Position',[100,100,1400,600])
subplot(2,3,1)
yyaxis left
plot(long_delta(1:Tbound),'LineWidth',2)
hold on
yyaxis right
plot(long_beta(1:Tbound),'LineWidth',2)
plot(long_ERN(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)),'-k','LineWidth',2)
legend('delta','beta','ERN')
lgd = legend;
lgd.Location = 'northwest';
title('Delta (1-30 week)')

subplot(2,3,2)
yyaxis left
plot(long_N(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)),'LineWidth',2)
hold on
plot(long_I(1:Tbound),'-k','LineWidth',2)
ylabel('# of People')
hold on
yyaxis right
ylabel('New Death')
plot(long_dD(1:Tbound),'LineWidth',2)
legend('N','I','D')
lgd = legend;
lgd.Location = 'northwest';
title('New Cases / Infected / Cumlative Death (1-30 week)')

subplot(2,3,3)
yyaxis left
plot(long_alpha(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)),'LineWidth',2)
hold on
yyaxis right
plot(long_beta(1:Tbound),'LineWidth',2)
hold on
plot(long_beta_tilde(1:Tbound,(abs(TH - TH_index(1)) < 0.0001)),'-k','LineWidth',2)
legend('alpha','beta','beta tilde')
lgd = legend;
lgd.Location = 'southeast';
title('Alpha (1-30 week)')

subplot(2,3,4)
yyaxis left
plot(long_delta(Tbound+1:end),'LineWidth',2)
hold on
yyaxis right
plot(long_beta(Tbound+1:end),'LineWidth',2)
plot(long_ERN(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)),'-k','LineWidth',2)
xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
legend('delta','beta','ERN')
lgd = legend;
lgd.Location = 'northeast';
title('Delta since 31st week')

subplot(2,3,5)
yyaxis left
plot(long_N(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)),'LineWidth',2)
hold on
plot(long_I(Tbound+1:end),'-k','LineWidth',2)
ylabel('# of People')
hold on
yyaxis right
ylabel('New Death')
plot(long_dD(Tbound+1:end),'LineWidth',2)
legend('N','I','dD')
xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
title('New Cases / Infected / Cumlative Death since 31st week')

subplot(2,3,6)
yyaxis left
plot(long_alpha(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)),'LineWidth',2)
hold on
yyaxis right
plot(long_beta(Tbound+1:end),'LineWidth',2)
hold on
plot(long_beta_tilde(Tbound+1:end,(abs(TH - TH_index(1)) < 0.0001)),'-k','LineWidth',2)
xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
legend('alpha','beta','beta tilde')
title('Alpha since 31st week')




%% figure for vaccine path

figure(200)
set(gcf,'Position',[100,100,1200,800])
Data_vac = VT(:,:);
medStuff1_nv = [V1_w; Data_vac(:,5)]; % In future, you should change V1_w to V1past_med.
medStuff2_nv = [V2_w; Data_vac(:,6)]; % In future, you should change V2_w to V2past_med.
elderly1_nv = [zeros(Tdata,1); Data_vac(:,1)];
elderly2_nv = [zeros(Tdata,1); Data_vac(:,2)];
others1_nv = [zeros(Tdata,1); Data_vac(:,3)];
others2_nv = [zeros(Tdata,1); Data_vac(:,4)];
% medStuff1_nv = [V1_w; Data_vac(:,1)]; % In future, you should change V1_w to V1past_med.
% medStuff2_nv = [V2_w; Data_vac(:,2)]; % In future, you should change V2_w to V2past_med.
% elderly1_nv = [zeros(Tdata,1); Data_vac(:,3)];
% elderly2_nv = [zeros(Tdata,1); Data_vac(:,4)];
% others1_nv = [zeros(Tdata,1); Data_vac(:,5)];
% others2_nv = [zeros(Tdata,1); Data_vac(:,6)];
% medStuff1_nv = zeros(Tdata+size(Data_vac,1),1);
% medStuff2_nv = zeros(Tdata+size(Data_vac,1),1);
% elderly1_nv = [zeros(Tdata,1); Data_vac(:,1)];
% elderly2_nv = [zeros(Tdata,1); Data_vac(:,2)];
% others1_nv = [zeros(Tdata,1); Data_vac(:,3)];
% others2_nv = [zeros(Tdata,1); Data_vac(:,4)];



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

% ----------- effect of vaccines on mortality rate  ----------------------------
subplot(3,2,3)
plot([delta; deltaT(:)]./delta_average,'LineWidth',1.5);

xtickangle(45)
xticks(find(WeekNumber==1))
xlim([Tdata+1 Tdata+56])
xticklabels(MonthWeekEN(WeekNumber==1))

title('致死率（現在のレベルで標準化）','FontName',fn)
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;

% ------------ Compare mean(delta) and weighted(delta) ------------- %%

subplot(3,2,4)
% figure(172)
sw_delta = 31;
delta_path_1 = zeros(Tdata+SimPeriod,1);
delta_path_2 = zeros(Tdata+SimPeriod,1);
for i_delta = sw_delta:1:size(delta,1)
    delta_sample = delta(i_delta-RetroPeriod+1:i_delta);
    delta_path_1(i_delta) = mean(delta_sample);
    delta_path_2(i_delta) = sum(delta_sample.*(I(i_delta-RetroPeriod+1:i_delta)/sum(I(i_delta-RetroPeriod+1:i_delta))));
end
d(1)=plot(delta_path_1,'LineWidth',1);
hold on
d(2)=plot(delta_path_2,'LineWidth',1);
ylim([0 0.05])
xlim([sw_delta Tdata])
xtickangle(45)
xticks(find(WeekNumber==1))
xticklabels(MonthWeekEN(WeekNumber==1))
legend([d(1),d(2)],["Mean Delta", "Average Delta weighted by I"])
title('\delta')
xline(Tdata,'LineWidth',1,'HandleVisibility','off');

subplot(3,2,5)
% figure(173)
sw_delta = 25;
for i = 11:2:23
    delta_path_1 = zeros(Tdata+SimPeriod,1);
    for i_delta = sw_delta:1:size(delta,1)
        delta_sample = delta(i_delta-i+1:i_delta);
        delta_path_1(i_delta) = mean(delta_sample);
    end
    if i == 17
        d(i) = plot(delta_path_1,'-r','LineWidth',1);
    else
        d(i) = plot(delta_path_1,'--','LineWidth',0.5);
    end
    hold on
end
% ylim([0 0.05])
xlim([sw_delta Tdata])
xtickangle(45)
xticks(find(WeekNumber==1))
xticklabels(MonthWeekEN(WeekNumber==1))
legend([d(11:2:23)],string([11:2:23]),'FontSize',6,'Location','NorthEast')
title('Mean \delta')
xline(Tdata,'LineWidth',1,'HandleVisibility','off');


subplot(3,2,6)
% figure(174)
sw_delta = 25;
for i = 11:2:23
    delta_path_2 = zeros(Tdata+SimPeriod,1);
    for i_delta = sw_delta:1:size(delta,1)
        delta_sample = delta(i_delta-i+1:i_delta);
        delta_path_2(i_delta) = sum(delta_sample.*(I(i_delta-i+1:i_delta)/sum(I(i_delta-i+1:i_delta))));
    end
    if i == 17
        d(i) = plot(delta_path_2,'-r','LineWidth',1);
    else
        d(i) = plot(delta_path_2,'--','LineWidth',0.5);
    end
    hold on
end
% ylim([0 0.05])
xlim([sw_delta Tdata])
xtickangle(45)
xticks(find(WeekNumber==1))
xticklabels(MonthWeekEN(WeekNumber==1))
legend([d(11:2:23)],string([11:2:23]),'FontSize',6,'Location','NorthEast')
title("Average \delta weighted by I")
xline(Tdata,'LineWidth',1,'HandleVisibility','off');
