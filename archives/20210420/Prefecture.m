% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
iPC=0;
if iPC==1
    %     home = '\Users\shcor\Dropbox\Website\Codes\';
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
%     home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata\Website\Codes\';
%     home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/codes/';
%     home = '/Users/shotaro/Dropbox/fujii_nakata/Website/codes/';
      home = '/Users/machikohei/Dropbox/fujii_nakata/Website/Codes/';
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 1;      % save back data
% in the "Figure" folder
fs = 16;            % common font size for many figures
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
%======================================================================%

SimCases=4;
SimPeriod = 52;        % simulation period in weeks

PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka','Aichi','Fukuoka','Hyogo'};
GDPVector = [106,36,23,21,40,41,20,20];
% th_on_vector = [1250,600,300,350,400,250,250,350];
% th_on_vector = [1250,700,400,350,1000,350,350];   4/12
th_on_vector = [1250,700,400,350,1000,350,350];
th_on_vector2 = [1250,700,400,350,1000,350,350]*2; %[1750,700,400,350,500,350,350];
%If delta = Terminal state ... 3000? linearly decreasing
%th_on_vector2 = th_on_vector * 1.6;
%ERN_now_vector = [0.99,0.94,0.96,0.955,0.9,0.9,0.99,0.99];
%ERN_on_vector = [0.99,0.99,1.00,0.99,0.97,0.97,1.02,0.99];
%ERN_on_vector = [0.99,0.99,0.99,0.99,0.99,0.99,0.99];
ERN_now_vector = [0.99,0.94,0.96,0.955,0.9,0.9,0.99];
ERN_on_vector = [0.9,0.99,1.00,0.99,0.85,0.97,1.02]; %
ERN_on_scenario1 = [0.8,0.8,0.8,0.8,0.8,0.8,0.8]; %0.8
ERN_on_scenario2 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]; %0.7
TL_vector = {250:10:350,70:10:100,80:10:130,50:10:100,50:10:100,20:10:100,10:10:100,50:10:100};
% TL_index_vector = {[250,350],[200,80],[260,110],[200,80],[300,40],[100,20],[120,10]};
TL_index_vector = {[150,350],[200,80],[260,110],[200,80],[300,40],[100,20],[120,10]};
ERNON_vector1 = {0.5:0.1:1.0,400:100:800,200:100:600,300:50:500,170:20:270,200:50:400,200:50:400,300:50:500};
% ERNON_vector1 = {0.99};
ERNON_index_vector1 = {[1.0,0.7],[600,800],[300,600],[350,500],[190,250],[250,400],[250,400],[350,500]};
ERNON_vector41 = {0.5:0.1:0.9,400:100:800,200:100:600,300:50:500,200:50:600,200:50:400,200:50:400,300:50:500};
ERNON_index_vector41 = {[0.8, 0.5],[600,800],[300,600],[350,500],[250,400],[400,600],[250,400],[350,500]};
ERNON_vector42_1 = {[0.3,0.4,0.5,0.55,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74],400:100:800,200:100:600,300:50:500,300:100:1000,200:50:400,200:50:400,300:50:500};
ERNON_vector42_2 = {[0.4,0.5,0.55,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74],400:100:800,200:100:600,300:50:500,300:100:1000,200:50:400,200:50:400,300:50:500};
ERNON_vector42_3 = {[0.5,0.55,0.6,0.62,0.64,0.66,0.68,0.69,0.7,0.72,0.74,0.76,0.78],400:100:800,200:100:600,300:50:500,300:100:1000,200:50:400,200:50:400,300:50:500};
% If th_off = 150 ... use the middle one
% If th_off = 50 ... use the bottom one
ERNON_index_vector42 = {[0.7,0.5],[600,800],[300,600],[350,500],[400,600],[400,600],[300,600],[350,500]};
THON_vector1 = {500:100:1000,400:100:800,200:100:600,300:50:500,170:20:270,200:50:400,200:50:400,300:50:500};
THON_index_vector1 = {[900,700],[600,800],[300,600],[350,500],[190,250],[250,400],[250,400],[350,500]};
THON_vector41 = {500:250:2000,400:100:800,200:100:600,300:50:500,200:50:600,200:50:400,200:50:400,300:50:500};
THON_index_vector41 = {[1250,2000],[600,800],[300,600],[350,500],[250,400],[400,600],[250,400],[350,500]};
THON_vector42 = {[900:300:1500,2000:1000:4000],400:100:800,200:100:600,300:50:500,300:100:1000,200:50:400,200:50:400,300:50:500};
THON_index_vector42 = {[1500,3000],[600,800],[300,600],[350,500],[400,600],[400,600],[300,600],[350,500]};
% TL_vector = {240:10:300,70:10:100,80:10:130,50:10:100,0:0.1:1,0:0.1:1,0:0.1:1,0:0.1:1};
% TL_index_vector = {[240,300],[100,80],[130,110],[100,80],[100,40],[50,20],[60,10],[60,10]};
% th_off_vector = [250,80,110,80,40,20,10,10]; %Also change TL_index
th_off_vector = {[500,350],[200,80],[260,110],[200,80],[400,40],[100,20],[120,10]};
%th_off_vector41 = {[250,350],[100,80],[130,110],[100,80],[100,40],[50,20],[60,10]};
% th_off_vector42 = {[250,350],[100,80],[80,110],[60,80],[100,40],[50,20],[60,10]}; %Also change TL_index
th_off_vector42 = {[150,350],[100,80],[80,110],[60,80],[100,40],[50,20],[60,10]}; %Also change TL_index
% th_off_vector = {[280,150,250,280],[50,150,250],[50,150,250],[50,150,250],[50,150,250],[50,150,250],[50,150,250],[50,150,250]};%Also change TL_index
% Tokyo 280, 250, 200, 150, 100, 50
DR_vector = {0:2:16,0:2:16,0:2:16,0:2:16,0:2:16,0:12,0:2:16,0:2:16};
DR_index = [4,12]; %0->4,8->12 gradual 3/15
% rho_vector= [0.85,0.9,0.9,0.8,0.9,0.8,0,0];
% impute_periods_vector =[45,22,25,20,45,30,30,30];
% alpha_off_vector = [0.024,0.0290,0.0280,0.0449,0.0290,0.0380,0.028,0.025];
rho_vector= [0.85,0.9,0.9,0.8,0.9,0.8,0];
impute_periods_vector =[45,22,25,20,45,30,30];
alpha_off_vector = [0.024, 0.0290, 0.0280, 0.0449, 0.0290, 0.0380, 0.028];
%rho_vector= [1.0,1.0,1.0,1.0,1.0,1.0,1.0];
%DR_index line 367
retro_ub = 17;
retro_lb = 17;
%betaT_temp_ini = 0.289; %so that 6 weeks average = betaT(1)*1.2
betaT_temp_ini = -0.289/2; %Osaka ... twice as much as this number
% 0.2 ... 35人, 0.3 ... 60人, 0.25 ... 50人, 0.185 & 0.5 ... 30人
beta_rho = 0.85;
% 0.3, 0.5 のときTL_vector_index [280,420]
% 0.25, 0.5 のときTL_vector_index [280,380]
% 0.2, 0.6 のときTL_vector_index [280,360]
tt = 12; %Showing previous t periods for the plot

%Parameters of variants
var_infection = 0.5;  % variant's infection rate (compared to nomral infection rate) [0.5, 0.7]
var_infection_vector = [0.0:0.1:0.7];
var_infection_index = [0.3,0.5];
var_ss = 1.0;         % steady-state share of variant
% v_scale = 10;
%var_growth_vec = [0.1695,0.5182]; %Feb 2021
%var_growth_vec = [0.3526,0.4480]; %April 9 2021
var_growth_vec = [0.47, 0.47]; %April 12nd, 2021
% var_growth_vec = [0, 0.597]; %April 15, 2021 Tokyo, kansenken
% var_growth_vec = [0, 0.366]; %April 15, 2021 Osaka, kansenken

% var_initial_vector = [0.000282,0.008786,0.007703,0.000395,0.000967,0,0,0.110650,0.011697]; %Up to 3/16
%  var_initial_vector = [0.000685,0.009607,0.008621,0.008905,0.047151,0,0,0.133262,0.015965]; % Up to 3/23
%   var_initial_vector = [0.000614,0.005794,0.007463,0.009519,0.020688,0,0,0,0.014124]; % Up to 3/30
%var_initial_vector = [0.000596,0.005711,0.007347,0.009359,0.017937,0,0,0,0.014124]; % Up to 3/30
%var_initial_vector = [0.0083,0.0319,0.0089,0.0364,0.0825,0,0,0.133262,0.015965]; % 4/12
% var_initial_vector = [0.0210, 0.0778, 0.0225, 0.0882, 0.1871, 0, 0, 0,0]; % 4/12
% var_initial_vector = [0.1495, 0, 0, 0, 0.7340, 0, 0, 0,0]; % 4/12
var_initial_vector = [0.2195, 0.2234, 0.0663, 0.2624, 0.8153, 0, 0, 0,0]; % 4/19
% var_initial_vector = [0.2886, 0, 0, 0, 0.8301, 0, 0, 0,0]; % 4/12, kansenken

%=========== use this code to extrapolate var_initial_vector ==================%
% var_initial_vector = [0.099, 0.101, 0.027, 0.122, 0.633, 0, 0, 0,0]; % 4/19
% pindex = 5;
% var_initial = var_initial_vector(pindex);
% var_growth = 0.47;
% logit_initial = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
% var_intercept = logit_initial;
% var_share = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
% var_share(2)
%=====================================================%

% v_scale_vec = [0.0270,0.0929,0.2053,0.0877,0.330,0.1,0.1,0.4404,0.1]; %Based on 3/1-3/7 reported values; Aichi and Fukuoka are set to 10%
%  v_scale_vec = [0.03101,0.07192,0.15111,0.16345,0.19290,0.36209,0.38316,0.38657,0.21943]; %Three weeks average up to 3/7 reported values
%v_scale_vec = [0.03458,0.08280,0.14784,0.15309,0.22642,0,0,0,0.23680]; %Three weeks average up to 3/15
v_scale_vec = ones(1,length(PrefVector)); %4/12 -
%scaleA = 0.25; %Scenario A ... initial vector is multiplied by scaleA
%scaleB = 0.5; %Scenario B ... initial vector is multiplied by scaleB
scaleA = 1;
scaleB = 1;
vec1 = [1,5,10,20,50];
vec2 = [1,2.5,4,6.5,8];

% vaccine pace
Vsimple = 0; % 0 for vaccine_distribution; 1 for vaccine_distribution_simple
PF = 1; % 0 for AZ, 1 for PF
Vgradual = 1; % 0 for flat Vpath, 1 for gradually increasing Vpath
% if Vsimple == 0
%     VP = [3:1:17]; %[0];
%     VP_index = [3,8]; %[0,0];
% else
%     VP = [3]; %[0];
%     VP_index = [3,8]; %[0,0];
% end
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


beta_option_index = [42];
% if Vsimple == 0
%     beta_option_index = [1,42];
% else
%     beta_option_index = [42];
% end
%beta_option = 42; %{1,3,41,42} 1 ... baseline, 3 ... 気の緩み,
% 41 ... 変異株シナリオ1(アメリカのgrowth rate), 42 ... 変異株シナリオ2(イギリスのgrowth rate)


for pindex = 1:5  %:length(PrefVector)
    if pindex == 5
        betaT_temp_ini = -0.289;
    end
    NPath_beta = NaN(SimPeriod,40,SimCases,length(beta_option_index));
    AlphaM_beta = NaN(40,SimCases,length(beta_option_index));
    DM_beta = NaN(40,SimCases,length(beta_option_index));
    waves_beta = NaN(40,SimCases,length(beta_option_index));
    
    for bb = [beta_option_index]
        beta_option = bb;
        if beta_option == 42
            th_off_vec = th_off_vector42{pindex};
        elseif beta_option == 41
            th_off_vec = th_off_vector41{pindex};
        else
            th_off_vec = th_off_vector{pindex};
        end
        %th_off_vec = th_off_vector{pindex};
        TH_th = NaN(40,SimCases,length(th_off_vec));
        waves_th= NaN(40,SimCases,length(th_off_vec));
        for th_off_index = 1:1 %1:length(th_off_vec)
            close all
            if th_off_index == 1
                ERNON_vector42 = ERNON_vector42_1;
            elseif th_off_index == 2
                ERNON_vector42 = ERNON_vector42_2;
            elseif th_off_index == 3
                ERNON_vector42 = ERNON_vector42_3;
            end
            
            %====================== Model parameter values ======================%
            pref = PrefVector{pindex};        % prefecture to be analyzed
            prefGDP = GDPVector(pindex);
            % gamma = 7/5;          % recovery rate from Covid
            %         gamma = 7/8;          % recovery rate from Covid
            gamma = 7/10;          % recovery rate from Covid
            k = 2;                 % exponent of (1-h*alpha)
            hconstant = 1;         % 0 = without intercept, 1 = with intercept for h regression
            %            SimPeriod = 52;        % simulation period in weeks
            medical_start_date = datetime(2021,3,18);
            elderly_start_date = datetime(2021,5,13);
            retro2 = 13;
            
            % VacStartDate = "Apr-01";             % time until the start of vaccination process
            %
            % VacDuration = 12;       % time until vaccination pace reaches its steady pace
            RetroPeriod = 17;      % retroactive periods used to estimate gamma and delta
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
            % Columns: 1 = date, 2 = new positive, 3 = death, 4 = mobility, 5 = GDP, 6 = population
            dateD = Data(:,1) + 21916;
            N = Data(:,2);
            dD = Data(:,3);
            M = Data(:,4);
            % GDP = Data(:,5);
            POP = Data(:,6);
            GDP = Data(:,7);
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
            xaxis_vec = 0:1:length(SimDate);
            xaxis_vec2 = 1:1:length(SimDate)+tt;
            xaxis_vec3 = 1:1:Tdata+tt;
            
            
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
            RetroH = TdataGDP-4;
            %     RetroH = 15;
            if isempty(find(SimDateEN == medical_start_date)) == 0
                medical_start = find(SimDateEN == medical_start_date);
            else
                medical_start = 1;
            end
            elderly_start = find(SimDateEN == elderly_start_date);
            VacStart = find(SimDateEN == datetime(2021,4,1));
            End2020 = find(dateEN == datetime(2021,1,7));
            
            %--- Construct weekly vaccine data ---%
            if iPC == 1
                vaccine = importdata([home 'vaccine_daily.xls']);   % Import daily vaccine data
            else
                vaccine = importdata([home 'vaccine_daily.xls']); 
                %vaccine = importdata([home 'vaccine_daily.csv']);   % Import daily vaccine data
            end
            dateV = datetime(vaccine.textdata(2:end,1),'InputFormat','yyyy/MM/dd');
            %dateV = datetime(vaccine.textdata(2:end,1),'InputFormat','MM/dd/yyyy');
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
            % potentialGDP(1) = (100/(1.0122))*(1.0063^(1/12)); %GDP in 2019/12 = 100
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
            if iPC==1
                xticklabels(MonthWeekJP(xtick1))
            else
                xticklabels(MonthWeekJP(xtick1))
            end
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
            %     for i = 1:Tdata
            %         S(i+1)=S(i)-N(i);
            %         I(i+1)=I(i)+N(i)-gamma*I(i)-dD(i);
            %         R(i+1)=R(i)+gamma*I(i);
            %         D(i+1)=D(i)+dD(i);
            %         if i > TdataGDP
            %             GDP(i) = referenceGDP(i)*(1-alpha(i));
            %         end
            %     end
%             effectiveness = 0.9;
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
            % beta_tilde = -POP0*((S(2:Tdata+1)-S(1:Tdata))./(S(1:Tdata).*I(1:Tdata)));   % overall infection rate
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
            % alpha_off = mean(alpha((dateEN >= datetime(2020,9,4)) & (datetime(2020,11,26)>= dateEN )));
            % alpha_off = alpha_off_vector(pindex);
            % alpha_off = mean(alpha((dateEN >= datetime(2020,9,4)) & (datetime(2020,11,26)>= dateEN ))); % output loss without the state of emergency
            %     alpha_off = 0;
            %     if pindex == 6
            %         alpha_off = alpha_off - 0.02;
            %     end
            
            %Minimum after summer
            %alpha_off = min(alpha(dateEN >= datetime(2020,8,6)));
            
            InitialValues = [S(end),I(end),R(end),D(end)];
            
            %--- Construct time series of parameters ---%
            gammaT = gamma*ones(SimPeriod,1);
            % beta_sample = beta(end-RetroPeriod+1:end);
            beta_r = 0;
            for retrop = retro_lb:retro_ub
                beta_r = beta_r + mean(beta(end-retrop+1:end));
            end
            beta_avg = beta_r/(retro_ub-retro_lb+1);
            % beta_sample = beta(TdataGDP-RetroPeriod+1:TdataGDP);
            %     betaT = mean(beta_sample)*ones(SimPeriod,1);
            %     betaT is determined in beta_option, so this line should be removed.
            
            delta_sample = delta(end-RetroPeriod+1:end);
            % delta_average = mean(delta_sample);
            delta_average = sum(delta_sample.*(I(end-RetroPeriod+1:end)/sum(I(end-RetroPeriod+1:end))));
            
            %             pace = ps*3500000;
            pace = ps*3600000;
            %                 pace2 = pace*1.27;
            %                 pace3 = pace*1.86;
            vacpath = zeros(SimPeriod,1);
            if Vgradual == 1
                % case 1: gradually increasing Vpath
                vacpath(1:10) = (pace/10):(pace/10):pace;
                vacpath(11:end) = pace*ones(42,1);
            else
                % case 2: flat Vpath
                vacpath(5:end) = pace*ones(48,1);
            end
            %             pace_m = ps*900000;
                        medical = ps*4700000*0.8;
%                         if medical_start == 1
%                             medical = medical - sum(V1_w);
%                             if medical < 0
%                                 medical = 0;
%                             end
%                         end
            elderly = ps*36000000*0.8;
            ordinary = (125710000-36000000-4700000)*ps*0.8;
%             ordinary = (125710000-36000000)*ps*0.8;
            elderly_total = ps*36000000;
            
%             V_mat = zeros(SimPeriod,size(VP,2));
%             VT_mat = zeros(SimPeriod,4,size(VP,2));
%             delta_mat = zeros(SimPeriod,size(VP,2));
            %     [V,deltaT,VT,real_pace] = vaccine_path(pace,pace_m,SimPeriod,medical,elderly,ordinary,medical_start,elderly_start,3,0.9,delta_average,elderly_total);
            %             [V,deltaT,VT,real_pace] = vaccine_path_medical_off(pace,pace_m,SimPeriod,medical,elderly,ordinary,medical_start,elderly_start,16,effectiveness,delta_average,elderly_total,V1_w,V2_w);
            if Vsimple == 0
%                 [V,deltaT,VT] = vaccine_distribution(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2,3);
                [V,deltaT,VT] = vaccine_distribution_medical(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2,ps,POP0,3);
            else
                [V,deltaT,VT] = vaccine_distribution_simple(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2);
            end
            
            %betaT_temp_ini = 0;
            %         betaT_temp_ini = mean(beta(end-5+1:end))/beta_avg - 1; % calculate the AR(1) shock size of beta for each prefecture: initial beta = average beta in the last 5 wks
            betaT = beta_avg*ones(SimPeriod,1);
            betaT_temp = betaT_temp_ini;
            betaT(1,1) = (1+betaT_temp)*betaT(1,1);
            for i = 2:length(betaT)
                betaT_temp = betaT_temp * beta_rho;
                betaT(i) = betaT(i) * (1+betaT_temp);
            end
            if beta_option == 1
                betaT = beta_avg*ones(SimPeriod,1);
                %Normal
                % betaT = mean(beta_sample)*ones(SimPeriod,1);
                %         betaT = beta_avg*ones(SimPeriod,1);
                %         betaT(1,1) = (1+betaT_temp)*betaT(1,1);
                %         for i = 2:length(betaT)
                %             betaT_temp = betaT_temp * beta_rho;
                %             betaT(i) = betaT(i) * (1+betaT_temp);
                %         end
            elseif beta_option == 2
                %Covid-Variant
                betaT = 1.1*mean(beta_sample)*ones(SimPeriod,1);
            elseif beta_option == 3
                %Less Social Distancing during March w3,w4 and April w1,w2
                % betaT = mean(beta_sample)*ones(SimPeriod,1);
                %betaT = beta_avg*ones(SimPeriod,1);
                % betaT(1:5) =  1.2*beta_avg; % 1.3880 = mean (Dec10 - Jan7, 5 weeks)
                %1.2183 ... 2020/11/26 - 2021/1/7
                %1.1581 ... 2020/11/12 - 2020/12/24
                %(sum(ERN((dateEN >= datetime(2020,11,26)) &
                %(datetime(2020,12,31)>= dateEN ))) + 0.5*1.9750)/6.5 = 1.1601
                %(2020/11/26 - 2020/12/31 + 0.5 * 2021/1/7)
                
                % betaT(1:9) = 1.08*beta_avg;  %1.1492 ... 2020/11/05 - 2020/12/31
                % betaT(1:3) = 1.3*beta_avg; % to achieve ERN = 1.418, % we have to adjust the inside of betaT(3:5) every week so that beta increases from March 22 for 3 weeks
                %         beta_multiplier_vec = linspace(0,5,1000);
                %         [min_val, min_ind] = min(abs(max(ERN(end-15:end))*ones(1,1000)-(S(end)/POP0).*(((1+(h(2)/h(1))*alpha_off).^k).*beta_multiplier_vec.*betaT(1))./(gammaT(1)+deltaT(1))));
                %         beta_multiplier=beta_multiplier_vec(min_ind);
                %         betaT(4:6) =  beta_multiplier*betaT(4:6);
            elseif beta_option == 4
                var_initial = 0.1;     % initial variant share
                var_growth = 0.2;     % weekly growth parameter for logit model
                var_intercept = log(var_initial/(1-var_initial));
                var_share = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
                beta_bar = beta_avg/(1+var_infection*var_ss*var_initial);
                betaT = beta_bar*(1+var_infection*var_ss*var_share);
                
            elseif beta_option == 41
                var_initial = scaleA *var_initial_vector(pindex) / v_scale_vec(pindex);
                var_growth = var_growth_vec(1);
                logit_initial = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
                %         var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept
                var_intercept = logit_initial;
                var_share = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
                if var_initial > 0.7
                    beta_bar = beta_avg * 0.95; % if var_initial = 0.8 (var_share = 0.4 in the last month)
                elseif var_initial > 0.5 && var_initial < 0.7
                    beta_bar = beta_avg * 0.96; % if var_initial = 0.6 (var_share = 0.3 in the last month)
                elseif var_initial > 0.3 && var_initial < 0.5
                    beta_bar = beta_avg * 0.975; % if var_initial = 0.4 (var_share = 0.2 in the last month)
                elseif var_initial > 0.1 && var_initial < 0.3
                    beta_bar = beta_avg * 0.99; % if var_initial = 0.2 (var_share = 0.1 in the last month)
                else
                    beta_bar = beta_avg * 1;    % if var_initial = 0.0
                end
                %                 beta_bar = beta_avg/(1+var_infection*var_initial);
                if SimCases ~= 4
                    betaT = beta_bar*(1+var_infection*var_share);
                end
                
                %Assuming that current beta is higher
                %                 betaT_temp = betaT_temp_ini;
                %                 betaT2(1,1) = (1+betaT_temp)*betaT(1,1);
                %                 for i = 2:length(betaT)
                %                     betaT_temp = betaT_temp * beta_rho;
                %                     betaT2(i,1) = betaT(i) * (1+betaT_temp);
                betaT2 = betaT;
                %                 end
                
                %             betaT2(1:5) =  1.2*betaT2(1:5); % 1.3880 = mean (Dec10 - Jan7, 5 weeks)
                
                plot_var_share(1,1) = var_initial;
                plot_var_share(2:length(SimDate)+1,1) = var_share(:,1);
                % plot_betaT(1,1) = beta_avg;
                
                
                plot_betaT(1:tt,1) = beta(end-tt+1:end);
                plot_betaT(tt+1:length(SimDate)+tt,1) = betaT(:,1);
                
                % plot_betaT2(1,1) = beta_avg*(1+betaT_temp);
                plot_betaT2(1:tt,1) = beta(end-tt+1:end);
                plot_betaT2(tt+1:length(SimDate)+tt,1) = betaT2(:,1);
                
                plot_betaT3(1:tt,1) = beta(end-tt+1:end);
                plot_betaT3(tt+1:length(SimDate)+tt,1) = beta_avg*ones(length(SimDate),1);
                
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
                        figure(4101)
                    elseif l == 2
                        figure(4102)
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
                    % plot(xaxis_vec2(tt+1:end), plot_betaT2(tt+1:end,1),'-b','LineWidth',2.0)
                    plot(xaxis_vec2,plot_betaT3,'-k','LineWidth',2)
                    ax = gca;
                    ax.YAxis.FontSize = 20;
                    ax.XAxis.FontSize = 20;
                    ax.YAxis.Exponent = 0;
                    ytickformat('%,0.2f')
                    xticks(find(WeekNumber_Sim2==1))
                    xtickangle(45)
                    if l == 1
                        % legend('Raw Infection Rate with Variants', 'Raw Infection Rate w/o Variants','Past Average','FontSize',10)
                        legend('Raw Infection Rate with Variants','Past Average','FontSize',10)
                        title('Projected path of beta','FontSize',20,'FontWeight','normal')
                        xticklabels(MonthWeekEN_Sim2(WeekNumber_Sim2==1))
                        xlabel('Weeks','FontSize',20)
                        ylabel('Infection Rate','FontSize',20,'FontName',fn)
                    elseif l == 2
                        legend('感染率(変異株成長ケース)','過去の平均','FontSize',10,'FontName',fn)
                        % legend('感染率1(変異株成長ケース)','感染率2(変異株成長ケース)','過去の平均','FontSize',10,'FontName',fn)
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
                    saveas(figure(4101),[home 'Figures/' char(pref) '/beta_path' sprintf('%.0f', beta_option) '.png']);
                    saveas(figure(4102),[home 'Figures/' char(pref) '/beta_path' sprintf('%.0f', beta_option) '_jp.png']);
                end
                
                betaT = betaT2;
                
                
            elseif beta_option == 42
                var_initial = scaleB *var_initial_vector(pindex) / v_scale_vec(pindex);
                var_growth = var_growth_vec(2);
                logit_initial = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
                %         var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept
                var_intercept = logit_initial;
                var_share = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
                if var_initial > 0.7
                    beta_bar = beta_avg * 0.95; % if var_initial = 0.8 (var_share = 0.4 in the last month)
                elseif var_initial > 0.5 && var_initial < 0.7
                    beta_bar = beta_avg * 0.96; % if var_initial = 0.6 (var_share = 0.3 in the last month)
                elseif var_initial > 0.3 && var_initial < 0.5
                    beta_bar = beta_avg * 0.975; % if var_initial = 0.4 (var_share = 0.2 in the last month)
                elseif var_initial > 0.1 && var_initial < 0.3
                    beta_bar = beta_avg * 0.99; % if var_initial = 0.2 (var_share = 0.1 in the last month)
                else
                    beta_bar = beta_avg * 1;    % if var_initial = 0.0
                end
                %                 beta_bar = beta_avg/(1+var_infection*var_initial);
                if SimCases ~= 4
                    betaT = beta_bar*(1+var_infection*var_share);
                end
                %Assuming that current beta is higher
                %                 betaT_temp = betaT_temp_ini;
                %                 betaT2(1,1) = (1+betaT_temp)*betaT(1,1);
                %                 for i = 2:length(betaT)
                %                     betaT_temp = betaT_temp * beta_rho;
                %                     betaT2(i,1) = betaT(i) * (1+betaT_temp);
                betaT2 = betaT;
                %                 end
                
                %             betaT2(1:5) =  1.2*betaT2(1:5); % 1.3880 = mean (Dec10 - Jan7, 5 weeks)
                
                plot_var_share(1,1) = var_initial;
                plot_var_share(2:length(SimDate)+1,1) = var_share(:,1);
                % plot_betaT(1,1) = beta_avg;
                tt = 12; %Showing previous t periods for the plot
                
                plot_betaT(1:tt,1) = beta(end-tt+1:end);
                plot_betaT(tt+1:length(SimDate)+tt,1) = betaT(:,1);
                
                % plot_betaT2(1,1) = beta_avg*(1+betaT_temp);
                plot_betaT2(1:tt,1) = beta(end-tt+1:end);
                plot_betaT2(tt+1:length(SimDate)+tt,1) = betaT2(:,1);
                
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
                    % plot(xaxis_vec2(tt+1:end), plot_betaT2(tt+1:end,1),'-r','LineWidth',2.0)
                    plot(xaxis_vec2,plot_betaT3,'-k','LineWidth',2)
                    ax = gca;
                    ax.YAxis.FontSize = 20;
                    ax.XAxis.FontSize = 20;
                    ax.YAxis.Exponent = 0;
                    ytickformat('%,0.2f')
                    xticks(find(WeekNumber_Sim2==1))
                    xtickangle(45)
                    %             if l == 1
                    %                 legend('Raw Infection Rate with Variants', 'Raw Infection Rate w/o Variants','Past Average','FontSize',10)
                    %                 title('Projected path of beta','FontSize',20,'FontWeight','normal')
                    %                 xticklabels(MonthWeekEN_Sim2(WeekNumber_Sim2==1))
                    %                 xlabel('Weeks','FontSize',20)
                    %                 ylabel('Infection Rate','FontSize',20,'FontName',fn)
                    %             elseif l == 2
                    %                 legend('感染率1(変異株成長ケース)','感染率2(変異株成長ケース)','過去の平均','FontSize',10,'FontName',fn)
                    %                 xticklabels(MonthWeekJP_Sim2(WeekNumber_Sim2==1))
                    %                 title('βの推移','FontSize',20,'FontWeight','normal','FontName',fn)
                    %                 xlabel('週','FontSize',20)
                    %                 ylabel('感染率','FontSize',20,'FontName',fn)
                    %             end
                    if l == 1
                        % legend('Raw Infection Rate with Variants', 'Raw Infection Rate w/o Variants','Past Average','FontSize',10)
                        legend('Raw Infection Rate with Variants','Past Average','FontSize',10)
                        title('Projected path of beta','FontSize',20,'FontWeight','normal')
                        xticklabels(MonthWeekEN_Sim2(WeekNumber_Sim2==1))
                        xlabel('Weeks','FontSize',20)
                        ylabel('Infection Rate','FontSize',20,'FontName',fn)
                    elseif l == 2
                        legend('感染率(変異株成長ケース)','過去の平均','FontSize',10,'FontName',fn)
                        % legend('感染率1(変異株成長ケース)','感染率2(変異株成長ケース)','過去の平均','FontSize',10,'FontName',fn)
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
                
                betaT = betaT2;
                
            end
            
            % VacPace = (POP0/125710000)*0.9*(3000000/2);  % number of vaccinations per week
            % V = zeros(SimPeriod,1);
            % V(VacStart-1:VacStart+VacDuration-1) = 0:VacPace/VacDuration:VacPace;
            % V(VacStart+VacDuration:end) = VacPace;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Projection parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            th_on = th_on_vector(pindex)*7;         % threshold to place the state of emergency (daily new infected persons in Tokyo = 1750)
            th_on2 = th_on_vector2(pindex)*7;         % threshold to place the state of emergency (daily new infected persons in Tokyo = 1750)
            %     ERN_on = ERN_on_vector(pindex);
            %     altA = (((ERN_on*(POP0/S(end))*((gammaT(1)+deltaT(1))/betaT(1))).^(1/k))-1)*(h(1)/h(2));
            th_off = th_off_vec(th_off_index)*7;
            ERN_on = ERN_on_vector(pindex);
            if beta_option == 41
                ERN_on = ERN_on_scenario1(pindex);
            elseif beta_option == 42
                ERN_on = ERN_on_scenario2(pindex);
            end
            ERN_now = ERN_now_vector(pindex);
            altA_now = (((ERN_now*(POP0/S(end))*((gammaT(1)+deltaT(1))/beta_avg)).^(1/k))-1)*(h(1)/h(2));
            altA_on = (((ERN_on*(POP0/S(end))*((gammaT(1)+deltaT(1))/beta_avg)).^(1/k))-1)*(h(1)/h(2));
            ERNCheck = (S(end)/POP0).*(((1+(h(2)/h(1))*alpha_off).^k).*beta_avg)./(gammaT(1)+deltaT(1));
            %SimCases=4;
            DMat = nan(SimCases,40);
            AlphaMat = nan(SimCases,40);
            AlphaPath = nan(SimPeriod,40,SimCases);
            NPath = nan(SimPeriod,40,SimCases);
            SimERN = nan(SimPeriod,40,SimCases);
            BackDataN = zeros(SimPeriod+8,3,SimCases);
            BackDataAlpha = zeros(SimPeriod+8,3,SimCases);
            BackDataERN = zeros(SimPeriod+8,3,SimCases);
            BackDataDA = zeros(40,3,SimCases);
            
            %---- 1. Different thresholds to lift the state of emergency ---%
            % TL = TL_vector{pindex};
            % TL_index = TL_index_vector{pindex};
            % for i = 1:length(TL)
            %         if pindex <= 4
            %             [DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData,NPath(:,i,1),SimERN(:,i,1)] = Covid_projection_control_gradual_alt_2(InitialValues,altA_now,altA_on,alpha_off,th_on,th_off,TL(i)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,4); %1 -> 4 gradual 3/15
            %             % [DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData,NPath(:,i,1),SimERN(:,i,1)] = Covid_projection_control_gradual_alt(InitialValues,altA_now,altA_on,alpha_off,th_on,TL(i)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,4); %1 -> 4 gradual 3/15
            %         elseif pindex > 4
            %         [DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData,NPath(:,i,1),SimERN(:,i,1)] = Covid_projection_control_gradual_off(InitialValues,altA_on,alpha_off,th_on,TL(i)*7,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,4,alpha(end)); %1 -> 4 gradual 3/15
            %         end
            % end
            %---------------------------------------------------------------%
            
            %---- . Different thresholds to initiate the state of emergency ---%
            TL = TL_vector{pindex};
            TL_index = TL_index_vector{pindex};
            if beta_option == 1 || beta_option == 3
                ON = THON_vector1{pindex}*7;
                ON_index = THON_index_vector1{pindex};
            elseif beta_option == 41
                ON = THON_vector41{pindex}*7;
                ON_index = THON_index_vector41{pindex};
            elseif beta_option == 42
                ON = THON_vector42{pindex}*7;
                ON_index = THON_index_vector42{pindex};
            end
            for i = 1:length(ON)
                %[DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData,NPath(:,i,1),SimERN(:,i,1)] = Covid_projection_control_gradual_off(InitialValues,altA_on,alpha_off,ON(i),TL_index(1)*7,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+DR_index(1),alpha(end));
                [DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData(:,:,i,1),NPath(:,i,1),SimERN(:,i,1)] = Covid_projection_control_gradual_off(InitialValues,altA_on,alpha_off,ON(i),th_off,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+DR_index(1),alpha(end));
            end
            ON = ON/7;
            
            %---- 2. Different durations for economic recovery ---%
            DR = DR_vector{pindex};
            %     DR_index = [4,8]; %0->4,8->12 gradual 3/15
            %     if pindex == 3 || pindex == 4 || pindex == 7
            %         DR_index = [4,12];
            %     end
            for i = 1:length(DR)
                %         if pindex <= 4
                %             [DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData,NPath(:,i,2),SimERN(:,i,2)] = Covid_projection_control_gradual_alt_2(InitialValues,altA_now,altA_on,alpha_off,th_on,th_off,TL_index(1)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,1+DR(i));
                %             % [DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData,NPath(:,i,2),SimERN(:,i,2)] = Covid_projection_control_gradual_alt(InitialValues,altA_now,altA_on,alpha_off,th_on,TL_index(1)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,1+DR(i));
                %         elseif pindex > 4
                %             [DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData,NPath(:,i,2),SimERN(:,i,2)] = Covid_projection_control_gradual_off(InitialValues,altA_on,alpha_off,th_on,TL_index(1)*7,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+DR(i),alpha(end));
                %[DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData,NPath(:,i,2),SimERN(:,i,2)] = Covid_projection_control_gradual_off_threshold(InitialValues,altA_on,alpha_off,th_on,th_on2,TL_index(1)*7,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+DR(i),alpha(end));
                [DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData(:,:,i,2),NPath(:,i,2),SimERN(:,i,2)] = Covid_projection_control_gradual_off_threshold(InitialValues,altA_on,alpha_off,th_on,th_on2,th_off,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+DR(i),alpha(end));
                %         end
            end
            
            
            %---- 3. Different ERN on ---%
            %     TL = TL_vector{pindex};
            %     TL_index = TL_index_vector{pindex};
            %             if beta_option == 1 || beta_option == 3
            %                 ERNON = (((ERNON_vector1{pindex}.*(POP0/S(end)).*((gammaT(1)+deltaT(1))./beta_avg)).^(1/k))-1).*(h(1)/h(2));
            %                 ERNON_index = ERNON_index_vector1{pindex};
            %             elseif beta_option == 41
            %                 ERNON = (((ERNON_vector41{pindex}.*(POP0/S(end)).*((gammaT(1)+deltaT(1))./beta_avg)).^(1/k))-1).*(h(1)/h(2));
            %                 ERNON_index = ERNON_index_vector41{pindex};
            %             elseif beta_option == 42
            %                 ERNON = (((ERNON_vector42{pindex}.*(POP0/S(end)).*((gammaT(1)+deltaT(1))./beta_avg)).^(1/k))-1).*(h(1)/h(2));
            %                 ERNON_index = ERNON_index_vector42{pindex};
            %             end
            %             for i = 1:length(ERNON)
            %                 % [DMat(3,i),AlphaMat(3,i),AlphaPath(:,i,3),SimData,NPath(:,i,3),SimERN(:,i,3)] = Covid_projection_control_gradual_off(InitialValues,ERNON(i),alpha_off,th_on,TL_index(1)*7,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+DR_index(1),alpha(end));
            %                 %[DMat(3,i),AlphaMat(3,i),AlphaPath(:,i,3),SimData,NPath(:,i,3),SimERN(:,i,3)] = Covid_projection_control_gradual_off_threshold(InitialValues,ERNON(i),alpha_off,th_on,th_on2,TL_index(1)*7,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+DR_index(1),alpha(end));
            %                 [DMat(3,i),AlphaMat(3,i),AlphaPath(:,i,3),SimData(:,:,i,3),NPath(:,i,3),SimERN(:,i,3)] = Covid_projection_control_gradual_off_threshold(InitialValues,ERNON(i),alpha_off,th_on,th_on2,th_off,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+DR_index(1),alpha(end));
            %             end
            
            %-----------------------------------------------------%
            
%             if SimCases == 3
%             %------------ 3. Different vaccine pace -------------%
%             % %             VP = [25,50,75,100:50:600];
%             %             VP = [3:1:16];
%             % %             VP_index = [200,600,350];
%             %             VP_index = [6,15];
%             for i = 1:length(VP)
%                 lag = VP(i);
%                 %                 [V,deltaT,VT,real_pace] = vaccine_path(ps*VP(i)*10000,ps*900000,SimPeriod,medical,elderly,ordinary,medical_start,elderly_start,3,0.9,delta_average,elderly_total);
%                 %                 [DMat(3,i),AlphaMat(3,i),AlphaPath(:,i,3),SimData,NPath(:,i,3),SimERN(:,i,3)] = Covid_projection_control_gradual(InitialValues,altA,alpha_off,th_on,TL_index(2)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,1);
%                 
%                 if Vsimple == 0
%                     [V,deltaT,VT] = vaccine_distribution(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2,lag);
%                     [V_mat(:,i),delta_mat(:,i),VT_mat(:,:,i)] = vaccine_distribution(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2,lag);
%                 else
%                     [V,deltaT,VT] = vaccine_distribution_simple(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2);
%                     [V_mat(:,i),delta_mat(:,i),VT_mat(:,:,i)] = vaccine_distribution_simple(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2);
%                 end
%                 [DMat(3,i),AlphaMat(3,i),AlphaPath(:,i,3),SimData(:,:,i,3),NPath(:,i,3),SimERN(:,i,3)] = Covid_projection_control_gradual_off_threshold(InitialValues,altA_on,alpha_off,th_on,th_on2,th_off,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+6,alpha(end));
%                 
%             end
%             %-----------------------------------------------------%
%             end 
            
            if SimCases == 4
            %%---- 4  Different infection rate (only for variant scenario) ---%
            for i = 1:length(var_infection_vector)
                var_infection = var_infection_vector(i);
                betaT = beta_bar*(1+var_infection*var_share);
                [DMat(4,i),AlphaMat(4,i),AlphaPath(:,i,4),SimData(:,:,i,4),NPath(:,i,4),SimERN(:,i,4)] = Covid_projection_control_gradual_off_threshold(InitialValues,altA_on,alpha_off,th_on,th_on2,th_off,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+6,alpha(end));
            end
            %%---------------------------------------------------------------%
            end
            
            %%---- 4  Different thresholds to lift the state of emergency (with alpha_off = 2020 Autumn level) ---%
            % Alt_alpha_off = min(alpha(dateEN >= datetime(2020,8,6))); %Alternative value of alpha_off
            % TL = 100:50:400;
            % TL_index = [400,150];
            % for i = 1:length(TL)
            %     [DMat(4,i),AlphaMat(4,i),AlphaPath(:,i,4),SimData,NPath(:,i,4),SimERN(:,i,4)] = Covid_projection_control_gradual(InitialValues,altA,Alt_alpha_off,th_on,TL(i)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,1);
            % end
            %%---------------------------------------------------------------%
            
            
            %%---- 5  Different thresholds to lift the state of emergency (with higher beta) ---%
            % Alt_betaT = 1.7*(mean(beta_sample)*ones(SimPeriod,1)); Alternative value of betaT
            % TL = 100:50:400;
            % TL_index = [400,150];
            % for i = 1:length(TL)
            %     [DMat(5,i),AlphaMat(5,i),AlphaPath(:,i,5),SimData,NPath(:,i,5),SimERN(:,i,5)] = Covid_projection_control_gradual(InitialValues,altA,alpha_off,th_on,TL(i)*7,Alt_betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,1);
            % end
            %%---------------------------------------------------------------%
            
            
            
            %           plot_betaT(1:tt,1) = beta(end-tt+1:end);
            %         plot_betaT(tt+1:length(SimDate)+tt,1) = betaT(:,1);
            %
            %         % plot_betaT2(1,1) = beta_avg*(1+betaT_temp);
            %         plot_betaT2(1:tt,1) = beta(end-tt+1:end);
            %         plot_betaT2(tt+1:length(SimDate)+tt,1) = betaT2(:,1);
            %
            %         plot_betaT3(1:tt,1) = beta(end-tt+1:end);
            %         plot_betaT3(tt+1:length(SimDate)+tt,1) = beta_avg*ones(length(SimDate),1);
            
            
            minAlpha = min(AlphaMat,[],'all');
            
            for y = SimCases %Figure(12100)  %SimCases
                if y == 1
                    minAlpha1 = min(AlphaMat(y,:));
                    if beta_option == 1
                        save('alpha_min1.mat', 'minAlpha1')
                    end
                    minAlphaMat1 = load('alpha_min1.mat','minAlpha1');
                    minAlpha = minAlphaMat1.minAlpha1;
                elseif y == 2
                    minAlpha2 = min(AlphaMat(y,:));
                    if beta_option == 1
                        save('alpha_min2.mat', 'minAlpha2')
                    end
                    minAlphaMat2 = load('alpha_min2.mat','minAlpha2');
                    minAlpha = minAlphaMat2.minAlpha2;
                elseif y == 3
                    minAlpha3 = min(AlphaMat(y,:));
                    if beta_option == 1
                        save('alpha_min3.mat', 'minAlpha3')
                    end
                    minAlphaMat3 = load('alpha_min3.mat','minAlpha3');
                    minAlpha = minAlphaMat3.minAlpha3;
                elseif y == 4
                    minAlpha4 = min(AlphaMat(y,:));
%                     if beta_option == 1
                        save('alpha_min4.mat', 'minAlpha4')
%                     end
                    minAlphaMat4 = load('alpha_min4.mat','minAlpha4');
                    minAlpha = minAlphaMat4.minAlpha4;
                end
                if y == 1
                    %if (y==1 || y==4)
                    %TH = TL;
                    %TH_index = TL_index;
                    TH = ON;
                    TH_index = ON_index;
                elseif y == 2
                    TH = DR;
                    TH_index = DR_index;
                elseif y == 3
                    TH = VP;
                    TH_index = VP_index;
                    %                     if beta_option == 1 || beta_option == 3
                    %                         TH = ERNON_vector1{pindex};
                    %                     elseif beta_option == 41
                    %                         TH = ERNON_vector41{pindex};
                    %                     elseif beta_option == 42
                    %                         TH = ERNON_vector42{pindex};
                    %                     end
                    %                     TH_index = ERNON_index;
                elseif y == 4
                    TH = var_infection_vector;
                    TH_index = var_infection_index;
                end
                TH_th(1:length(TH),y,th_off_index) = TH';
                AlphaM = AlphaMat(y,:);
                AlphaM_th(:,th_off_index) = AlphaM;
                AlphaM = AlphaM(~isnan(AlphaM));
                DM = DMat(y,:);
                DM_th(:,th_off_index) = DM;
                DM = DM(~isnan(DM));
                AlphaM = (AlphaM - minAlpha)*prefGDP*10000;
                AlphaM_th(:,th_off_index) = (AlphaM_th(:,th_off_index) - minAlpha)*prefGDP*10000;
                
                BackDataDA(1:length(TH),:,y) = [AlphaM',DM',TH'];
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
                    ATL = AlphaM(waves==0);
                    DTL = DM(waves==0);
                    % THTL = TH(waves==0);
                    THON = TH(waves==0);
                    ATLblue = AlphaM(TH==TH_index(2));
                    DTLblue = DM(TH==TH_index(2));
                elseif y == 2
                    ADR = AlphaM(waves==0);
                    DDR = DM(waves==0);
                    THDR = TH(waves==0);
                    ADRblue = AlphaM(TH==TH_index(2));
                    DDRblue = DM(TH==TH_index(2));
                elseif y == 3
                    AERN = AlphaM(waves==0);
                    DERN = DM(waves==0);
                    % THTL = TH(waves==0);
                    THON = TH(waves==0);
                    AERNblue = AlphaM(TH==TH_index(2));
                    DERNblue = DM(TH==TH_index(2));
                elseif y == 4
                    Ainfec = AlphaM(waves==0);
                    Dinfec = DM(waves==0);
                    % THTL = TH(waves==0);
                    THON = TH(waves==0);
                    Ainfecblue = AlphaM(TH==TH_index(2));
                    Dinfecblue = DM(TH==TH_index(2));
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
                    % set(gcf,'Position',[100,100,1200,500])
                    set(gcf,'Position',[100,100,1000,800])
                    subplot(3,2,1)
                    if y == 1 || y == 2
                        for i = 1:length(TH)
                            if abs(TH(i) - TH_index(1)) < 0.0001
                                plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                                BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                                BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                                BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                            elseif abs(TH(i) - TH_index(2)) < 0.0001
                                plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                                BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                                BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                                BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                            else
                                plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                            end
                            hold on
                        end
                    elseif y == 3 || y == 4
                        for i = 1:length(TH)
                            if abs(TH(i) - TH_index(1)) < 0.0001
                                plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.2f',TH(i)))
                                BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                                BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                                BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                            elseif abs(TH(i) - TH_index(2)) < 0.0001
                                plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.2f',TH(i)))
                                BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                                BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                                BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                            else
                                plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.2f',TH(i)))
                            end
                            hold on
                        end
                    end
                    %             elseif y == 3
                    %                 if l == 1
                    %                     for i = 1:length(TH)
                    %                         if TH(i) == TH_index(1)
                    %                             plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                             BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                    %                             BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                    %                             BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                    %                         elseif TH(i) == TH_index(2)
                    %                             plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                             BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                    %                             BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                    %                             BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                    %                         elseif TH(i) == TH_index(3)
                    %                             plot([N;NPath(:,i,y)]/7,'-k','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                             BackDataN(:,3,y) = [N(end-7:end);NPath(:,i,y)];
                    %                             BackDataAlpha(:,3,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                    %                             BackDataERN(:,3,y) = [ERN(end-7:end);SimERN(:,i,y)];
                    %                         else
                    %                             plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                    %                         end
                    %                         hold on
                    %                     end
                    %                 elseif l == 2
                    %                     for i = 1:length(TH)
                    %                         if TH(i) == TH_index(1)
                    %                             plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                    %                         elseif TH(i) == TH_index(2)
                    %                             plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                    %                         elseif TH(i) == TH_index(3)
                    %                             plot([N;NPath(:,i,y)]/7,'-k','LineWidth',2,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                    %                         else
                    %                             plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                    %                         end
                    %                         hold on
                    %                     end
                    %                 end
                    %             end
                    plot(N/7,'k','LineWidth',2,'HandleVisibility','off')
                    xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
                    %             if pindex <= 4
                    %                 xline(Tdata+find(SimDateEN == datetime(2021,3,18)),'--','LineWidth',1.5,'HandleVisibility','off');
                    %             end
                    %             if y == 1
                    %                 xline(Tdata+find(SimDateEN == datetime(2021,4,1)),'--','LineWidth',1.5,'HandleVisibility','off');
                    %             end
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
                    %             if y == 3
                    %                 lgd.Location = 'north';
                    %             end
                    xtickangle(45)
                    if beta_option == 41 || beta_option == 42
                        xlim([Tdata-7 Tdata+52])
                    else
                        xlim([Tdata-7 Tdata+28])
                    end
                    %             if y == 3
                    %                 xlim([Tdata-7 Tdata+28])
                    %             end
                    
                    %--- Number of cumulative deaths ---%
                    subplot(3,2,2)
                    plot(AlphaM(waves==0),DM(waves==0),'-go','LineWidth',2,'MarkerSize',10);
                    hold on
                    plot(AlphaM(waves==1),DM(waves==1),'-mo','LineWidth',2,'MarkerSize',10);
                    hold on
                    plot(AlphaM(waves==2),DM(waves==2),'-bo','LineWidth',2,'MarkerSize',10);
                    hold on
                    %             if y ~= 3
                    text(AlphaM,DM,string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
                    hold on
                    %             end
                    %             plot(AlphaM(waves==1),DM(waves==1),'-go','LineWidth',2,'MarkerSize',10)
                    %             hold on
                    %             plot(AlphaM(waves==3),DM(waves==3),'-mo','LineWidth',2,'MarkerSize',10)
                    %             hold on
                    %             plot(AlphaM(waves==5),DM(waves==5),'-bo','LineWidth',2,'MarkerSize',10)
                    %             hold on
                    %             text(AlphaM(1),DM(1),string(TH(1)),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14)
                    %             for i = 2:length(AlphaM)
                    %                 if AlphaM(i)~=AlphaM(i-1) && DM(i)~=DM(i-1)
                    %                     text(AlphaM(i),DM(i),string(TH(i)),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14)
                    % %                 elseif AlphaM(i)==AlphaM(i-1) && DM(i)==DM(i-1) && AlphaM(i)~=AlphaM(i-2) && DM(i)~=DM(i-2)
                    % %                     text(AlphaM(i),DM(i),string(TH(i)/7),'VerticalAlignment','top','HorizontalAlignment','left','FontSize',14)
                    % %                 elseif AlphaM(i)==AlphaM(i-1) && DM(i)==DM(i-1) && AlphaM(i)==AlphaM(i-2) && DM(i)==DM(i-2) && AlphaM(i)~=AlphaM(i-3) && DM(i)~=DM(i-3)
                    % %                     text(AlphaM(i),DM(i),string(TH(i)/7),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',14)
                    % %                 elseif AlphaM(i)==AlphaM(i-1) && DM(i)==DM(i-1) && AlphaM(i)==AlphaM(i-2) && DM(i)==DM(i-2)&& AlphaM(i)==AlphaM(i-3) && DM(i)==DM(i-3)
                    % %                     text(AlphaM(i),DM(i),string(TH(i)/7),'VerticalAlignment','top','HorizontalAlignment','right','FontSize',14)
                    %                 end
                    %                 hold on
                    %             end
                    %         text(AlphaM,DM,string(TH/7),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14)
                    %         hold on
                    %             text(AlphaM(TH==TH_index(1)),DM(TH==TH_index(1)),string(TH(TH==TH_index(1))/7),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14)
                    %             hold on
                    %             text(AlphaM(TH==TH_index(2)),DM(TH==TH_index(2)),string(TH(TH==TH_index(1))/7),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14)
                    %             hold on
                    scatter(AlphaM(TH==TH_index(1)),DM(TH==TH_index(1)),250,'red','filled');
                    hold on
                    scatter(AlphaM(TH==TH_index(2)),DM(TH==TH_index(2)),250,'blue','filled');
                    %             if y == 3
                    %                 hold on
                    %                 scatter(AlphaM(TH==TH_index(3)),DM(TH==TH_index(3)),250,'black','filled');
                    %             end
                    
                    if l == 1
                        xlabel('Output Loss (hundred million yen)','FontSize',10)
                        ylabel('Cumulative Deaths','FontSize',10)
                        title('Relationship between Covid-19 and output','FontSize',10,'FontWeight','normal')
                    elseif l == 2
                        %                 xlabel('経済損失 (億円)','FontSize',20)
                        %                 ylabel('累計死亡者数','FontSize',20)
                        %                 title('コロナ感染と経済の関係','FontSize',20,'FontWeight','normal')
                        xlabel('経済損失 (億円)','FontSize',10,'FontName',fn)
                        ylabel('累計死亡者数','FontSize',10,'FontName',fn)
                        title('コロナ感染と経済の関係','FontSize',10,'FontWeight','normal','FontName',fn)
                    end
                    xlim([0,inf])
                    %             if beta_option == 1
                    %                 if pindex == 4
                    %                     ylim([1000,1300])
                    %                 elseif pindex == 6
                    %                     ylim([600,900])
                    %                 elseif pindex == 7
                    %                     ylim([800,1100])
                    %                 end
                    %             end
                    %             if beta_option == 3
                    %                 if pindex == 7
                    %                     ylim([950,1250])
                    %                 end
                    %             end
                    xtickangle(45)
                    grid on
                    ax = gca;
                    ax.YAxis.FontSize = 10;
                    ax.XAxis.FontSize = 10;
                    ax.YAxis.Exponent = 0;
                    ax.XAxis.Exponent = 0;
                    ytickformat('%,6.0f')
                    
                    %         subplot(2,2,3)
                    %         for i = 1:length(TH)
                    %             if TH(i) == TH_index(1)
                    %                 plot((1-[alpha;AlphaPath(:,i,y)])*100,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                 if y == 1
                    %                     redalpha = (1-[alpha(end-7:end);AlphaPath(:,i,y)])*100;
                    %                 end
                    %             elseif TH(i) == TH_index(2)
                    %                 plot((1-[alpha;AlphaPath(:,i,y)])*100,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                 if y == 1
                    %                     bluealpha = (1-[alpha(end-7:end);AlphaPath(:,i,y)])*100;
                    %                 end
                    %             else
                    %                 plot((1-[alpha;AlphaPath(:,i,y)])*100,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                    %             end
                    %             hold on
                    %         end
                    %         plot((1-alpha)*100,'k','LineWidth',2,'HandleVisibility','off')
                    %         xticks(find(WeekNumber==1))
                    %         if l == 1
                    %             title('GDP','FontSize',20,'FontWeight','normal')
                    %             xticklabels(MonthWeekEN(WeekNumber==1))
                    %         elseif l == 2
                    %             title('経済活動の推移','FontSize',20,'FontWeight','normal')
                    %             xticklabels(MonthWeekJP(WeekNumber==1))
                    %         end
                    %         lgd = legend;
                    %         lgd.NumColumns = 2;
                    %         xtickangle(45)
                    %         xlim([Tdata-7 Tdata+29])
                    %         ax = gca;
                    %         ax.YAxis.FontSize = 20;
                    %         ax.XAxis.FontSize = 20;
                    %         ax.YAxis.Exponent = 0;
                    %         ytickformat('%,6.0f')
                    
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
                    %             for i = 1:length(TH)
                    %                 if TH(i) == TH_index(1)
                    %                     plot([1-alpha;1-AlphaPath(:,i,2)]*100,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                 elseif TH(i) == TH_index(2)
                    %                     plot([1-alpha;1-AlphaPath(:,i,2)]*100,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                 else
                    %                     plot([1-alpha;1-AlphaPath(:,i,2)]*100,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                    %                 end
                    %                 hold on
                    %             end %99.49, 95.87
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
                        ylabel('GDP','FontSize',10)
                        title('GDP Path','FontSize',10,'FontWeight','normal')
                    elseif l == 2
                        xlabel('Weeks','FontSize',10,'FontName',fn)
                        ylabel('GDP','FontSize',10,'FontName',fn)
                        title('GDP path','FontSize',10,'FontWeight','normal','FontName',fn)
                    end
                    
                    
                    
                    if l == 1
                        if beta_option == 41
                            figure(41*1000+11+y)
                        elseif beta_option == 42
                            figure(42*1000+11+y)
                        else
                            figure(11+y)
                        end
                    elseif l == 2
                        if beta_option == 41
                            figure(41*1000+(11+y)*10+1)
                        elseif beta_option == 42
                            figure(42*1000+(11+y)*10+1)
                        else
                            figure((11+y)*10+1)
                        end
                    end
                    set(gcf,'Position',[100,100,1200,500])
                    %set(gcf,'Position',[100,100,1200,1000])
                    subplot(1,2,1)
                    if y == 1 || y == 2
                        for i = 1:length(TH)
                            if abs(TH(i) - TH_index(1)) < 0.0001
                                plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                                BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                                BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                                BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                            elseif abs(TH(i) - TH_index(2)) < 0.0001
                                plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                                BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                                BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                                BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                            else
                                plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                            end
                            hold on
                        end
                    elseif y == 3
                        for i = 1:length(TH)
                            if abs(TH(i) - TH_index(1)) < 0.0001
                                plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.2f',TH(i)))
%                                 plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                                BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                                BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                                BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                            elseif abs(TH(i) - TH_index(2)) < 0.0001
                                plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.2f',TH(i)))
%                                 plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                                BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                                BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                                BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                            else
                                plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.2f',TH(i)))
%                                 plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                            end
                            hold on
                        end
                    elseif y == 4
                        for i = 1:length(TH)
                            if abs(TH(i) - TH_index(1)) < 0.0001
                                plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.2f',TH(i)+1))
                                BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                                BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                                BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                            elseif abs(TH(i) - TH_index(2)) < 0.0001
                                plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.2f',TH(i)+1))
                                BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                                BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                                BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                            else
                                plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.2f',TH(i)+1))
                            end
                            hold on
                        end
                    end
                    %             elseif y == 3
                    %                 if l == 1
                    %                     for i = 1:length(TH)
                    %                         if TH(i) == TH_index(1)
                    %                             plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                             BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                    %                             BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                    %                             BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                    %                         elseif TH(i) == TH_index(2)
                    %                             plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                             BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                    %                             BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                    %                             BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                    %                         elseif TH(i) == TH_index(3)
                    %                             plot([N;NPath(:,i,y)]/7,'-k','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                             BackDataN(:,3,y) = [N(end-7:end);NPath(:,i,y)];
                    %                             BackDataAlpha(:,3,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                    %                             BackDataERN(:,3,y) = [ERN(end-7:end);SimERN(:,i,y)];
                    %                         else
                    %                             plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                    %                         end
                    %                         hold on
                    %                     end
                    %                 elseif l == 2
                    %                     for i = 1:length(TH)
                    %                         if TH(i) == TH_index(1)
                    %                             plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                    %                         elseif TH(i) == TH_index(2)
                    %                             plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                    %                         elseif TH(i) == TH_index(3)
                    %                             plot([N;NPath(:,i,y)]/7,'-k','LineWidth',2,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                    %                         else
                    %                             plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                    %                         end
                    %                         hold on
                    %                     end
                    %                 end
                    %             end
                    plot(N/7,'k','LineWidth',2,'HandleVisibility','off')
                    xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
                    %             if pindex <= 4
                    %                 xline(Tdata+find(SimDateEN == datetime(2021,3,18)),'--','LineWidth',1.5,'HandleVisibility','off');
                    %             end
                    %if y == 1
                    %    xline(Tdata+find(SimDateEN == datetime(2021,4,1)),'--','LineWidth',1.5,'HandleVisibility','off');
                    %end
                    ax = gca;
                    ax.YAxis.FontSize = 20;
                    ax.XAxis.FontSize = 20;
                    ax.YAxis.Exponent = 0;
                    ytickformat('%,6.0f')
                    %                     if y == 3      % uncomment if you are doing the ERN comparison.
                    %                         yline(th_on/7,'-.r',sprintf('%.0f',th_on/7),'LineWidth',0.5,'HandleVisibility','off');
                    %                         if th_on ~= th_on2 && max(waves) > 1
                    %                             yline(th_on2/7,'--r',sprintf('%.0f',th_on2/7),'LineWidth',0.5,'HandleVisibility','off');
                    %                         end
                    %                         yline(th_off/7,'-.b',sprintf('%.0f',th_off/7),'LineWidth',0.5,'HandleVisibility','off');
                    %         %                 yticks(th_off/7);
                    %         %                 ytickformat('%,6.0f')
                    %                     end
                    xticks(find(WeekNumber==1))
                    if l == 1
                        title('Projected path of new cases','FontSize',20,'FontWeight','normal')
                        xticklabels(MonthWeekEN(WeekNumber==1))
                        %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
                    elseif l == 2
                        %title('新規感染者数の推移','FontSize',20,'FontWeight','normal')
                        title('新規感染者数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
                        xticklabels(MonthWeekJP(WeekNumber==1))
                    end
                    if y ~= 3 || Vsimple == 0
                        lgd = legend;
                    end
                    lgd.NumColumns = 2;
                    %             if y == 3
                    %                 lgd.Location = 'north';
                    %             end
                    xtickangle(45)
                    if beta_option == 41 || beta_option == 42
                        xlim([Tdata-7 Tdata+52])
                    else
                        xlim([Tdata-7 Tdata+28])
                    end
                    %             if y == 3
                    %                 xlim([Tdata-7 Tdata+28])
                    %             end
                    
                    %--- Number of cumulative deaths ---%
                    subplot(1,2,2)
                    plot(AlphaM(waves==0),DM(waves==0),'-go','LineWidth',2,'MarkerSize',10);
                    hold on
                    plot(AlphaM(waves==1),DM(waves==1),'-mo','LineWidth',2,'MarkerSize',10);
                    hold on
                    plot(AlphaM(waves==2),DM(waves==2),'-bo','LineWidth',2,'MarkerSize',10);
                    hold on
                    if y == 4
                        text(AlphaM,DM,string(TH+1),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
                    elseif y ~= 3 || Vsimple == 0
                        text(AlphaM,DM,string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
                    end
                    hold on
                    %             end
                    %             plot(AlphaM(waves==1),DM(waves==1),'-go','LineWidth',2,'MarkerSize',10)
                    %             hold on
                    %             plot(AlphaM(waves==3),DM(waves==3),'-mo','LineWidth',2,'MarkerSize',10)
                    %             hold on
                    %             plot(AlphaM(waves==5),DM(waves==5),'-bo','LineWidth',2,'MarkerSize',10)
                    %             hold on
                    %             text(AlphaM(1),DM(1),string(TH(1)),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14)
                    %             for i = 2:length(AlphaM)
                    %                 if AlphaM(i)~=AlphaM(i-1) && DM(i)~=DM(i-1)
                    %                     text(AlphaM(i),DM(i),string(TH(i)),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14)
                    % %                 elseif AlphaM(i)==AlphaM(i-1) && DM(i)==DM(i-1) && AlphaM(i)~=AlphaM(i-2) && DM(i)~=DM(i-2)
                    % %                     text(AlphaM(i),DM(i),string(TH(i)/7),'VerticalAlignment','top','HorizontalAlignment','left','FontSize',14)
                    % %                 elseif AlphaM(i)==AlphaM(i-1) && DM(i)==DM(i-1) && AlphaM(i)==AlphaM(i-2) && DM(i)==DM(i-2) && AlphaM(i)~=AlphaM(i-3) && DM(i)~=DM(i-3)
                    % %                     text(AlphaM(i),DM(i),string(TH(i)/7),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',14)
                    % %                 elseif AlphaM(i)==AlphaM(i-1) && DM(i)==DM(i-1) && AlphaM(i)==AlphaM(i-2) && DM(i)==DM(i-2)&& AlphaM(i)==AlphaM(i-3) && DM(i)==DM(i-3)
                    % %                     text(AlphaM(i),DM(i),string(TH(i)/7),'VerticalAlignment','top','HorizontalAlignment','right','FontSize',14)
                    %                 end
                    %                 hold on
                    %             end
                    %         text(AlphaM,DM,string(TH/7),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14)
                    %         hold on
                    %             text(AlphaM(TH==TH_index(1)),DM(TH==TH_index(1)),string(TH(TH==TH_index(1))/7),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14)
                    %             hold on
                    %             text(AlphaM(TH==TH_index(2)),DM(TH==TH_index(2)),string(TH(TH==TH_index(1))/7),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14)
                    %             hold on
                    scatter(AlphaM(abs(TH - TH_index(1)) < 0.0001),DM(abs(TH - TH_index(1)) < 0.0001),250,'red','filled');
                    hold on
                    scatter(AlphaM(abs(TH - TH_index(2)) < 0.0001),DM(abs(TH - TH_index(2)) < 0.0001),250,'blue','filled');
                    %             if y == 3
                    %                 hold on
                    %                 scatter(AlphaM(TH==TH_index(3)),DM(TH==TH_index(3)),250,'black','filled');
                    %             end
                    
                    if l == 1
                        xlabel('Output Loss (hundred million yen)','FontSize',20)
                        ylabel('Cumulative Deaths','FontSize',20)
                        title('Relationship between Covid-19 and output','FontSize',20,'FontWeight','normal')
                    elseif l == 2
                        %                 xlabel('経済損失 (億円)','FontSize',20)
                        %                 ylabel('累計死亡者数','FontSize',20)
                        %                 title('コロナ感染と経済の関係','FontSize',20,'FontWeight','normal')
                        xlabel('経済損失 (億円)','FontSize',20,'FontName',fn)
                        ylabel('累計死亡者数','FontSize',20,'FontName',fn)
                        title('コロナ感染と経済の関係','FontSize',20,'FontWeight','normal','FontName',fn)
                    end
                    if beta_option ~= 3 || betaT_temp_ini>=0
                        xlim([0,inf])
                    end
                    %             if beta_option == 1
                    %                 if pindex == 4
                    %                     ylim([1000,1300])
                    %                 elseif pindex == 6
                    %                     ylim([600,900])
                    %                 elseif pindex == 7
                    %                     ylim([800,1100])
                    %                 end
                    %             end
                    %             if beta_option == 3
                    %                 if pindex == 7
                    %                     ylim([950,1250])
                    %                 end
                    %             end
                    xtickangle(45)
                    grid on
                    ax = gca;
                    ax.YAxis.FontSize = 20;
                    ax.XAxis.FontSize = 20;
                    ax.YAxis.Exponent = 0;
                    ax.XAxis.Exponent = 0;
                    ytickformat('%,6.0f')
                    
                    %         subplot(2,2,3)
                    %         for i = 1:length(TH)
                    %             if TH(i) == TH_index(1)
                    %                 plot((1-[alpha;AlphaPath(:,i,y)])*100,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                 if y == 1
                    %                     redalpha = (1-[alpha(end-7:end);AlphaPath(:,i,y)])*100;
                    %                 end
                    %             elseif TH(i) == TH_index(2)
                    %                 plot((1-[alpha;AlphaPath(:,i,y)])*100,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                    %                 if y == 1
                    %                     bluealpha = (1-[alpha(end-7:end);AlphaPath(:,i,y)])*100;
                    %                 end
                    %             else
                    %                 plot((1-[alpha;AlphaPath(:,i,y)])*100,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                    %             end
                    %             hold on
                    %         end
                    %         plot((1-alpha)*100,'k','LineWidth',2,'HandleVisibility','off')
                    %         xticks(find(WeekNumber==1))
                    %         if l == 1
                    %             title('GDP','FontSize',20,'FontWeight','normal')
                    %             xticklabels(MonthWeekEN(WeekNumber==1))
                    %         elseif l == 2
                    %             title('経済活動の推移','FontSize',20,'FontWeight','normal')
                    %             xticklabels(MonthWeekJP(WeekNumber==1))
                    %         end
                    %         lgd = legend;
                    %         lgd.NumColumns = 2;
                    %         xtickangle(45)
                    %         xlim([Tdata-7 Tdata+29])
                    %         ax = gca;
                    %         ax.YAxis.FontSize = 20;
                    %         ax.XAxis.FontSize = 20;
                    %         ax.YAxis.Exponent = 0;
                    %         ytickformat('%,6.0f')
                    
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
                    TAD = table([titleAD;round(BackDataDA(1:length(TH),:,y))]);
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
                    end
                end
                
                for i = 1:length(TH)
                    NPath_beta(:,i,y,find(beta_option_index == beta_option)) = NPath(:,i,y);
                    AlphaM_beta(i,y,find(beta_option_index == beta_option)) = AlphaM(i);
                    DM_beta(i,y,find(beta_option_index == beta_option)) = DM(i);
                    waves_beta(i,y,find(beta_option_index == beta_option)) = waves(i);
                end
            end
            
            %     figure(200)
            %     plot(ATL,DTL,'-go','LineWidth',2,'MarkerSize',10);
            %     hold on
            %     plot(ADR,DDR,'-go','LineWidth',2,'MarkerSize',10);
            %     hold on
            %     text(ATL,DTL,string(THTL),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
            %     hold on
            %     text(ADR,DDR,string(THDR),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
            %     hold on
            %     scatter(ATLblue,DTLblue,250,'blue','filled');
            %     hold on
            %     scatter(ADRblue,DDRblue,250,'blue','filled');
            % %     xlabel('経済損失 (億円)','FontSize',20)
            % %     ylabel('累計死亡者数','FontSize',20)
            % %     title('コロナ感染と経済の関係','FontSize',20,'FontWeight')
            %     xlabel('経済損失 (億円)','FontSize',20,'FontName',fn)
            %     ylabel('累計死亡者数','FontSize',20,'FontName',fn)
            %     title('コロナ感染と経済の関係','FontSize',20,'FontWeight','normal','FontName',fn)
            %     xtickangle(45)
            %     grid on
            %     ax = gca;
            %     ax.YAxis.FontSize = 20;
            %     ax.XAxis.FontSize = 20;
            %     ax.YAxis.Exponent = 0;
            %     ax.XAxis.Exponent = 0;
            %     ytickformat('%,6.0f')
            
%             y = 2;
%             plot_ERN(1:tt,1) = ERN(end-tt+1:end);
%             plot_ERN(tt+1:length(SimDate)+tt,1) = SimERN(:,DR_index(1),2);
%             for l = 1:2
%                 if l == 1
%                     figure(5101)
%                 elseif l == 2
%                     figure(5102)
%                 end
%                 set(gcf,'Position',[100,100,1200,500])
%                 if beta_option == 41 || beta_option == 42
%                     subplot(1,2,1)
%                     plot(xaxis_vec,plot_var_share(:,1)*100,'-r','LineWidth',2)
%                     ax = gca;
%                     ax.YAxis.FontSize = 20;
%                     ax.XAxis.FontSize = 20;
%                     ax.YAxis.Exponent = 0;
%                     ytickformat('%,3.0f')
%                     xticks(find(WeekNumber_Sim==1))
%                     if l == 1
%                         % legend('Variant Share','FontSize',10)
%                         title('Projected path of variants','FontSize',20,'FontWeight','normal')
%                         xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
%                         xlabel('Weeks','FontSize',20)
%                         ylabel('Share of Variants(%)','FontSize',20)
%                     elseif l == 2
%                         % legend('変異株割合','FontName',fn,'FontSize',10)
%                         title('新規感染者数の推移','FontSize',20,'FontWeight','normal')
%                         title('変異株シェアの推移','FontSize',20,'FontWeight','normal','FontName',fn)
%                         xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
%                         xlabel('週','FontSize',20,'FontName',fn)
%                         ylabel('変異株割合(%)','FontSize',20,'FontName',fn)
%                     end
%                     xtickangle(45)
%                     xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
%                     xlim([0,length(SimDate)]);
%                     ylim([0, 100])
%                     % lgd = legend;
%                     % lgd.Location = 'northwest';
%                     
%                     subplot(1,2,2)
%                     %             plot(xaxis_vec2(1:tt+1),plot_ERN(1:tt+1),'-k','LineWidth',2)
%                     %             hold on
%                     %             plot(xaxis_vec2(tt+1:length(SimDate)+tt),plot_ERN(tt+1:length(SimDate)+tt),'-r','LineWidth',2)
%                     %             ax = gca;
%                     %             ax.YAxis.FontSize = 20;
%                     %             ax.XAxis.FontSize = 20;
%                     %             ax.YAxis.Exponent = 0;
%                     %             ytickformat('%,0.2f')
%                     %             xticks(find(WeekNumber_Sim2==1))
%                     %             xtickangle(45)
%                     %             if l == 1
%                     %                 legend('Past ERN','Projected ERN','FontSize',10)
%                     %                 title('Projected path of ERN','FontSize',20,'FontWeight','normal')
%                     %                 xticklabels(MonthWeekEN_Sim2(WeekNumber_Sim2==1))
%                     %                 xlabel('Weeks','FontSize',20)
%                     %                 ylabel('ERN','FontSize',20,'FontName',fn)
%                     %             elseif l == 2
%                     %                 legend('過去の実効再生産数','予測された実効再生産数','FontSize',10,'FontName',fn)
%                     %                 xticklabels(MonthWeekJP_Sim2(WeekNumber_Sim2==1))
%                     %                 title('実効再生産数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
%                     %                 xlabel('週','FontSize',20)
%                     %                 ylabel('ERN','FontSize',20,'FontName',fn)
%                     %             end
%                     %             xline(tt+1,'--','LineWidth',1.5,'HandleVisibility','off');
%                     %             xlim([0,length(SimDate)+tt]);
%                     %             lgd = legend;
%                     %             lgd.Location = 'southeast';
%                     for i = 1:length(DR)
%                         if DR(i) == DR_index(1)
%                             plot([ERN;SimERN(:,i,2)],'-r','LineWidth',2,'DisplayName',sprintf('%.0f',DR(i)))
%                         elseif DR(i) == DR_index(2)
%                             plot([ERN;SimERN(:,i,2)],'-b','LineWidth',2,'DisplayName',sprintf('%.0f',DR(i)))
%                         else
%                             plot([ERN;SimERN(:,i,2)],'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',DR(i)))
%                         end
%                         hold on
%                     end
%                     if l == 1
%                         plot(xaxis_vec3(Tdata-tt+1:Tdata+1),plot_ERN(1:tt+1),'-k','LineWidth',2,'DisplayName','Actual')
%                     elseif l == 2
%                         plot(xaxis_vec3(Tdata-tt+1:Tdata+1),plot_ERN(1:tt+1),'-k','LineWidth',2,'DisplayName','これまで')
%                     end
%                     grid on
%                     xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
%                     yline(1.0,'LineWidth',1.5,'HandleVisibility','off');
%                     xticks(find(WeekNumber==1))
%                     xtickangle(45)
%                     ytickformat('%,2.4f')
%                     % xticklabels(MonthWeekJP(WeekNumber==1))
%                     xlim([Tdata-11 inf])
%                     ax = gca;
%                     ax.YAxis.FontSize = 10;
%                     ax.XAxis.FontSize = 10;
%                     if l == 1
%                         xticklabels(MonthWeekEN(WeekNumber==1))
%                         xlabel('Weeks','FontSize',10)
%                         ylabel('ERN','FontSize',10)
%                         title('ERN Path','FontSize',10,'FontWeight','normal')
%                     elseif l == 2
%                         xticklabels(MonthWeekJP(WeekNumber==1))
%                         %                 xticklabels(MonthWeekJP_Sim2(WeekNumber_Sim2==1))
%                         title('実効再生産数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
%                         xlabel('週','FontSize',20)
%                         ylabel('ERN','FontSize',20,'FontName',fn)
%                         %                 xlabel('Weeks','FontSize',10,'FontName',fn)
%                         %                 ylabel('ERN','FontSize',10,'FontName',fn)
%                         %                 title('ERN path','FontSize',10,'FontWeight','normal','FontName',fn)
%                     end
%                     legend
%                 else
%                     for i = 1:length(DR)
%                         if DR(i) == DR_index(1)
%                             plot([ERN;SimERN(:,i,2)],'-r','LineWidth',2,'DisplayName',sprintf('%.0f',DR(i)))
%                         elseif DR(i) == DR_index(2)
%                             plot([ERN;SimERN(:,i,2)],'-b','LineWidth',2,'DisplayName',sprintf('%.0f',DR(i)))
%                         else
%                             plot([ERN;SimERN(:,i,2)],'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',DR(i)))
%                         end
%                         hold on
%                     end
%                     %            plot(xaxis_vec2(Tdata-tt+1:Tdata+1),plot_ERN(1:tt+1),'-k','LineWidth',2,'DisplayName','Actual')
%                     if l == 1
%                         plot(xaxis_vec3(Tdata-tt+1:Tdata+1),plot_ERN(1:tt+1),'-k','LineWidth',2,'DisplayName','Actual')
%                     elseif l == 2
%                         plot(xaxis_vec3(Tdata-tt+1:Tdata+1),plot_ERN(1:tt+1),'-k','LineWidth',2,'DisplayName','これまで')
%                     end
%                     grid on
%                     xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
%                     yline(1.0,'LineWidth',1.5,'HandleVisibility','off');
%                     xticks(find(WeekNumber==1))
%                     xtickangle(45)
%                     ytickformat('%,2.4f')
%                     xticklabels(MonthWeekJP(WeekNumber==1))
%                     xlim([Tdata-11 inf])
%                     ax = gca;
%                     ax.YAxis.FontSize = 10;
%                     ax.XAxis.FontSize = 10;
%                     %              if l == 1
%                     %                 xlabel('Weeks)','FontSize',10)
%                     %                 ylabel('ERN','FontSize',10)
%                     %                 title('ERN Path','FontSize',10,'FontWeight','normal')
%                     %             elseif l == 2
%                     %                 xlabel('Weeks','FontSize',10,'FontName',fn)
%                     %                 ylabel('ERN','FontSize',10,'FontName',fn)
%                     %                 title('ERN path','FontSize',10,'FontWeight','normal','FontName',fn)
%                     %              end
%                     if l == 1
%                         xticklabels(MonthWeekEN(WeekNumber==1))
%                         xlabel('Weeks','FontSize',10)
%                         ylabel('ERN','FontSize',10)
%                         title('ERN Path','FontSize',10,'FontWeight','normal')
%                     elseif l == 2
%                         xticklabels(MonthWeekJP(WeekNumber==1))
%                         %                 xticklabels(MonthWeekJP_Sim2(WeekNumber_Sim2==1))
%                         title('実効再生産数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
%                         xlabel('週','FontSize',20)
%                         ylabel('ERN','FontSize',20,'FontName',fn)
%                         %                 xlabel('Weeks','FontSize',10,'FontName',fn)
%                         %                 ylabel('ERN','FontSize',10,'FontName',fn)
%                         %                 title('ERN path','FontSize',10,'FontWeight','normal','FontName',fn)
%                     end
%                     legend
%                 end
%             end
            
%             y = 3;
%             for l = 1:2
%                 if l == 1
%                     if beta_option == 41
%                         figure(41*10000+11+y)
%                     elseif beta_option == 42
%                         figure(42*10000+11+y)
%                     else
%                         figure(11+y)
%                     end
%                 elseif l == 2
%                     if beta_option == 41
%                         figure(41*10000+(11+y)*10+1)
%                     elseif beta_option == 42
%                         figure(42*10000+(11+y)*10+1)
%                     else
%                         figure((11+y)*100+1)
%                     end
%                 end
%                 set(gcf,'Position',[100,100,1200,500])
%                 subplot(1,2,1)
%                 if y == 1 || y == 2
%                     for i = 1:length(TH)
%                         if TH(i) == TH_index(1)
%                             plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
%                             BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
%                             BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
%                             BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
%                         elseif TH(i) == TH_index(2)
%                             plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
%                             BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
%                             BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
%                             BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
%                         else
%                             plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
%                         end
%                         hold on
%                     end
%                 elseif y == 3
%                     for i = 1:length(TH)
%                         if TH(i) == TH_index(1)
%                             %                                 plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.2f',TH(i)))
%                             plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
%                             BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
%                             BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
%                             BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
%                         elseif TH(i) == TH_index(2)
%                             %                                 plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.2f',TH(i)))
%                             plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
%                             BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
%                             BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
%                             BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
%                         else
%                             %                                 plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.2f',TH(i)))
%                             plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
%                         end
%                         hold on
%                     end
%                 end
%                 plot(N/7,'k','LineWidth',2,'HandleVisibility','off')
%                 xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
%                 ax = gca;
%                 ax.YAxis.FontSize = 20;
%                 ax.XAxis.FontSize = 20;
%                 ax.YAxis.Exponent = 0;
%                 ytickformat('%,6.0f')
%                 xticks(find(WeekNumber==1))
%                 if l == 1
%                     title('Projected path of new cases','FontSize',20,'FontWeight','normal')
%                     xticklabels(MonthWeekEN(WeekNumber==1))
%                     %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
%                 elseif l == 2
%                     %title('新規感染者数の推移','FontSize',20,'FontWeight','normal')
%                     title('新規感染者数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
%                     xticklabels(MonthWeekJP(WeekNumber==1))
%                 end
%                 if Vsimple == 0
%                     lgd = legend;
%                 end
%                 lgd.NumColumns = 2;
%                 xtickangle(45)
%                 if beta_option == 41 || beta_option == 42
%                     xlim([Tdata-7 Tdata+52])
%                 else
%                     xlim([Tdata-7 Tdata+28])
%                 end
%                 %--- Number of cumulative deaths ---%
%                 subplot(1,2,2)
%                 plot(VP(waves==0),DM(waves==0),'-go','LineWidth',2,'MarkerSize',10);
%                 hold on
%                 plot(VP(waves==1),DM(waves==1),'-mo','LineWidth',2,'MarkerSize',10);
%                 hold on
%                 plot(VP(waves==2),DM(waves==2),'-bo','LineWidth',2,'MarkerSize',10);
%                 hold on
%                 if Vsimple == 0
%                     text(VP,DM,string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
%                 end
%                 hold on
%                 
%                 scatter(VP(TH==TH_index(1)),DM(TH==TH_index(1)),250,'red','filled');
%                 hold on
%                 scatter(VP(TH==TH_index(2)),DM(TH==TH_index(2)),250,'blue','filled');
%                 
%                 if l == 1
%                     xlabel('Interval between the 1st and 2nd dose','FontSize',20)
%                     ylabel('Cumulative Deaths','FontSize',20)
%                     %                 title('Relationship between Covid-19 and output','FontSize',20,'FontWeight','normal')
%                 elseif l == 2
%                     xlabel('ワクチン接種の間隔','FontSize',20,'FontName',fn)
%                     ylabel('累計死亡者数','FontSize',20,'FontName',fn)
%                     %                 title('コロナ感染と経済の関係','FontSize',20,'FontWeight','normal','FontName',fn)
%                 end
%                 if beta_option ~= 3 || betaT_temp_ini>=0
%                     xlim([0,inf])
%                 end
%                 xtickangle(45)
%                 grid on
%                 ax = gca;
%                 ax.YAxis.FontSize = 20;
%                 ax.XAxis.FontSize = 20;
%                 ax.YAxis.Exponent = 0;
%                 ax.XAxis.Exponent = 0;
%                 ytickformat('%,6.0f')
%                 
%             end %End of language loop = figure loop
            
            
            %--- Save figures ---%
            if figure_save == 1
                saveas(figure(5101),[home 'Figures/' char(pref) '/ERN_path' sprintf('%.0f', beta_option) '.png']);
                saveas(figure(5102),[home 'Figures/' char(pref) '/ERN_path' sprintf('%.0f', beta_option) '_jp.png']);
                saveas(figure(2),[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
                if beta_option == 41
                    %saveas(figure(41012),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '.png']);
                    %saveas(figure(41121),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_jp.png']);
                    saveas(figure(41012),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '.png']);
                    saveas(figure(41121),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '_jp.png']);
                    saveas(figure(41013),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '.png']);
                    saveas(figure(41131),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_jp.png']);
                    if SimCases>=3
                        saveas(figure(41014),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '.png']);
                        saveas(figure(41141),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '_jp.png']);
                    end
                elseif beta_option == 42
                    %saveas(figure(42012),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '.png']);
                    %saveas(figure(42121),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_jp.png']);
                    saveas(figure(42012),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '.png']);
                    saveas(figure(42121),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '_jp.png']);
                    saveas(figure(42013),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '.png']);
                    saveas(figure(42131),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_jp.png']);
                    if SimCases>=3
                        saveas(figure(42014),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '.png']);
                        saveas(figure(42141),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '_jp.png']);
                        saveas(figure(42015),[home 'Figures/' char(pref) '/VarInfection' sprintf('%.0f', beta_option) '.png']);
                        saveas(figure(42151),[home 'Figures/' char(pref) '/VarInfection' sprintf('%.0f', beta_option) '_jp.png']);
                    end
                else
                    %saveas(figure(12),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '.png']);
                    %saveas(figure(121),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_jp.png']);
                    saveas(figure(12),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '.png']);
                    saveas(figure(121),[home 'Figures/' char(pref) '/ThresholdsON' sprintf('%.0f', beta_option) '_jp.png']);
                    saveas(figure(13),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '.png']);
                    saveas(figure(131),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_jp.png']);
                    if SimCases>=3
                        saveas(figure(14),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '.png']);
                        saveas(figure(141),[home 'Figures/' char(pref) '/EmergencyState' sprintf('%.0f', beta_option) '_jp.png']);
                    end
                    %         saveas(figure(14),[home 'Figures/' char(pref) '/VaccineVariation' sprintf('%.0f', beta_option) '.png']);
                    %         saveas(figure(141),[home 'Figures/' char(pref) '/VaccineVariation' sprintf('%.0f', beta_option) '_jp.png']);
                    %         saveas(figure(200),[home 'Figures/' char(pref) '/ThresholdsAndGradual' sprintf('%.0f', beta_option) '_jp.png']);
                end
            end
        end %end of th_off vector
        
        
        
    end %end of beta_option loop
end %end of prefecture loop

%%

if length(beta_option_index) == 3
    l = 2;
    TH = DR;
    TH_index = DR_index;
    if l == 1
        figure(5)
    elseif l == 2
        figure(6)
    end
    set(gcf,'Position',[100,100,1200,1200])
    subplot(3,2,1)
    for i = 1:length(TH)
        if TH(i) == TH_index(1)
            plot([N;NPath_beta(:,i,y,1)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
        elseif TH(i) == TH_index(2)
            plot([N;NPath_beta(:,i,y,1)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
        else
            plot([N;NPath_beta(:,i,y,1)]/7,'--','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
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
    elseif l == 2
        title('新規感染者数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP(WeekNumber==1))
    end
    lgd = legend;
    lgd.NumColumns = 2;
    xtickangle(45)
    if beta_option == 41 || beta_option == 42
        xlim([Tdata-7 Tdata+52])
    else
        xlim([Tdata-7 Tdata+28])
    end
    %--- Number of cumulative deaths ---%
    AlphaM1 = AlphaM_beta(:,2,1);
    DM1 = DM_beta(:,2,1);
    waves1 = waves_beta(:,2,1);
    subplot(3,2,2)
    plot(AlphaM1(waves1==0),DM1(waves1==0),'-go','LineWidth',2,'MarkerSize',10);
    hold on
    plot(AlphaM1(waves1==1),DM1(waves1==1),'-mo','LineWidth',2,'MarkerSize',10);
    hold on
    plot(AlphaM1(waves1==2),DM1(waves1==2),'-bo','LineWidth',2,'MarkerSize',10);
    hold on
    text(AlphaM1(1:length(TH)),DM1(1:length(TH)),string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
    hold on
    scatter(AlphaM1(TH==TH_index(1)),DM1(TH==TH_index(1)),250,'red','filled');
    hold on
    scatter(AlphaM1(TH==TH_index(2)),DM1(TH==TH_index(2)),250,'blue','filled');
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
    
    subplot(3,2,3)
    for i = 1:length(TH)
        if TH(i) == TH_index(1)
            plot([N;NPath_beta(:,i,y,2)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
        elseif TH(i) == TH_index(2)
            plot([N;NPath_beta(:,i,y,2)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
        else
            plot([N;NPath_beta(:,i,y,2)]/7,'--','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
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
    elseif l == 2
        title('新規感染者数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP(WeekNumber==1))
    end
    lgd = legend;
    lgd.NumColumns = 2;
    xtickangle(45)
    if beta_option == 41 || beta_option == 42
        xlim([Tdata-7 Tdata+52])
    else
        xlim([Tdata-7 Tdata+28])
    end
    %--- Number of cumulative deaths ---%
    AlphaM2 = AlphaM_beta(:,2,2);
    DM2 = DM_beta(:,2,2);
    waves2 = waves_beta(:,2,2);
    subplot(3,2,4)
    plot(AlphaM2(waves2==0),DM2(waves2==0),'-go','LineWidth',2,'MarkerSize',10);
    hold on
    plot(AlphaM2(waves2==1),DM2(waves2==1),'-mo','LineWidth',2,'MarkerSize',10);
    hold on
    plot(AlphaM2(waves2==2),DM2(waves2==2),'-bo','LineWidth',2,'MarkerSize',10);
    hold on
    text(AlphaM2(1:length(TH)),DM2(1:length(TH)),string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
    hold on
    scatter(AlphaM2(TH==TH_index(1)),DM2(TH==TH_index(1)),250,'red','filled');
    hold on
    scatter(AlphaM2(TH==TH_index(2)),DM2(TH==TH_index(2)),250,'blue','filled');
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
    
    
    subplot(3,2,5)
    for i = 1:length(TH)
        if TH(i) == TH_index(1)
            plot([N;NPath_beta(:,i,y,3)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
        elseif TH(i) == TH_index(2)
            plot([N;NPath_beta(:,i,y,3)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
        else
            plot([N;NPath_beta(:,i,y,3)]/7,'--','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
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
    elseif l == 2
        title('新規感染者数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP(WeekNumber==1))
    end
    lgd = legend;
    lgd.NumColumns = 2;
    xtickangle(45)
    if beta_option == 41 || beta_option == 42
        xlim([Tdata-7 Tdata+52])
    else
        xlim([Tdata-7 Tdata+28])
    end
    %--- Number of cumulative deaths ---%
    AlphaM3 = AlphaM_beta(:,2,3);
    DM3 = DM_beta(:,2,3);
    waves3 = waves_beta(:,2,3);
    subplot(3,2,6)
    plot(AlphaM3(waves3==0),DM3(waves3==0),'-go','LineWidth',2,'MarkerSize',10);
    hold on
    plot(AlphaM3(waves3==1),DM3(waves3==1),'-mo','LineWidth',2,'MarkerSize',10);
    hold on
    plot(AlphaM3(waves3==2),DM3(waves3==2),'-bo','LineWidth',2,'MarkerSize',10);
    hold on
    text(AlphaM3(1:length(TH)),DM3(1:length(TH)),string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
    hold on
    scatter(AlphaM3(TH==TH_index(1)),DM3(TH==TH_index(1)),250,'red','filled');
    hold on
    scatter(AlphaM3(TH==TH_index(2)),DM3(TH==TH_index(2)),250,'blue','filled');
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
    
    
    
    
    
end


%%
if SimCases >= 3  * 100
    for l = 1:2
        if l == 1
            if beta_option == 41
                figure(41*1000+21+y)
            elseif beta_option == 42
                figure(42*1000+21+y)
            else
                figure(11+y)
            end
        elseif l == 2
            if beta_option == 41
                figure(41*1000+(21+y)*10+1)
            elseif beta_option == 42
                figure(42*1000+(21+y)*10+1)
            else
                figure((11+y)*20+1)
            end
        end
        set(gcf,'Position',[100,100,1000,800])
        %--- Number of cumulative deaths ---%
        linSg = {'-go','--go',':go','-.go','-g*','--g*',':g*','-.g*'};
        linSm = {'-mo','--mo',':mo','-.mo','-m*','--m*',':m*','-.m*'};
        linSb = {'-bo','--bo',':bo','-.bo','-b*','--b*',':b*','-.b*'};
        linSk = {'-ko','--ko',':ko','-.ko','-k*','--k*',':k*','-.k*'};
        for i = 1:length(th_off_vec)
            plot(AlphaM_th((waves_th(:,3,i)==0),i),DM_th((waves_th(:,3,i)==0),i),linSg{i},'LineWidth',2,'MarkerSize',10,'DisplayName',sprintf('%.0f',th_off_vec(i)));
            hold on
            plot(AlphaM_th((waves_th(:,3,i)==1),i),DM_th((waves_th(:,3,i)==1),i),linSm{i},'LineWidth',2,'MarkerSize',10,'DisplayName',sprintf('%.0f',th_off_vec(i)));
            hold on
            plot(AlphaM_th((waves_th(:,3,i)==2),i),DM_th((waves_th(:,3,i)==2),i),linSb{i},'LineWidth',2,'MarkerSize',10,'DisplayName',sprintf('%.0f',th_off_vec(i)));
            plot(AlphaM_th((waves_th(:,3,i)==3),i),DM_th((waves_th(:,3,i)==3),i),linSk{i},'LineWidth',2,'MarkerSize',10,'DisplayName',sprintf('%.0f',th_off_vec(i)));
            hold on
            text(AlphaM_th(:,i),DM_th(:,i),string(TH_th(:,3,i)),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
            hold on
            scatter(AlphaM_th((TH_th(:,3,i)==TH_index(1)),i),DM_th((TH_th(:,3,i)==TH_index(1)),i),250,'red','filled','HandleVisibility','off');
            hold on
            scatter(AlphaM_th((TH_th(:,3,i)==TH_index(2)),i),DM_th((TH_th(:,3,i)==TH_index(2)),i),250,'blue','filled','HandleVisibility','off');
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
        lgd = legend;
        lgd.Location = 'northwest';
        lgd.FontSize = fs;
    end
end

%% Figures added on 2021/03/25
% cumDataEffV = effectiveness*cumsum(V2_w);
cumDataEffV = E1*cumsum(V1_w) + (E2-E1)*cumsum(V2_w);
cumSimEffV = cumsum(V)+cumDataEffV(end);
cumEffV  = [cumDataEffV; cumSimEffV];

long_delta = [delta;deltaT];
long_beta = [beta;betaT];
long_gamma = [ones(Tdata,1)*gamma;gammaT];
for i=1:length(DR)
    long_alpha(:,i) = [alpha;AlphaPath(:,i,2)];
    long_beta_tilde(:,i) = [beta_tilde;betaT.*(1+h(2)/h(1).*AlphaPath(:,i,2)).^k];
    long_ERN(:,i) = [ERN;SimERN(:,i,2)];
    long_GDP(:,i) = [GDP;(1-AlphaPath(:,i,2)).*referenceGDP(Tdata+1:Tdata+SimPeriod)];
    long_N(:,i) = [N;NPath(:,i,2)];
end
long_S = [S(1:end);SimData(2:end,1,find(DR == 4),2)];
long_I = [I(1:end);SimData(2:end,2,find(DR == 4),2)];
long_gamI = [gamma*I(2:end);gammaT.*SimData(2:end,2,find(DR == 4),2)];
long_R = [R(1:end);SimData(2:end,3,find(DR == 4),2)];
long_D = [D(1:end);SimData(2:end,4,find(DR == 4),2)];
long_dD = [D(2:end)-D(1:end-1);SimData(2:end,4,find(DR == 4),2)-SimData(1:end-1,4,find(DR == 4),2)];

figure(1111)
set(gcf,'Position',[100,100,1000,800])
subplot(3,3,1)
plot(long_S)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('Suceptibles')
subplot(3,3,2)
plot(long_N(:,find(DR == DR_index(1))))
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
% subplot(3,3,9)
% plot(cumEffV/effectiveness)
% xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% title('Cumulative Vaccinated')

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
for i=1:length(DR)
    plot(long_alpha(:,i))
    hold on
end
hold off
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('Alpha')

subplot(3,3,4)
for i = 1:length(DR)
    plot(long_GDP(:,i))
    hold on
end
hold off
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('GDP to Reference GDP')

subplot(3,3,5)
for i = 1:length(DR)
    plot(long_ERN(:,i))
    hold on
end
hold off
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
title('ERN')

subplot(3,3,6)
for i = 1:length(DR)
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
for i = 1:length(DR)
    plot(long_ERN(30:end,i))
    hold on
end
hold off
xline(Tdata-30+1,'LineWidth',1.5,'HandleVisibility','off');
title('ERN since 31st week')

% Tbound = 30;
% figure(1113)
% title('Up to 30th week')
% set(gcf,'Position',[100,100,1200,800])
% subplot(3,3,1)
% plot(long_delta(1:Tbound))
% title('Delta')
%
% subplot(3,3,2)
% plot(long_beta(1:Tbound))
% hold on
% plot(beta_tilde(1:Tbound))
% plot(long_N(1:Tbound,1)./I(1:Tbound))
% legend('beta','beta tilde','N/I')
% lgd = legend;
% lgd.Location = 'southeast';
% title('Beta')
%
% subplot(3,3,3)
% plot(long_N(1:Tbound))
% title('New Cases')
%
% subplot(3,3,4)
% plot(long_I(1:Tbound))
% title('Infected')
%
% subplot(3,3,5)
% plot(long_dD(1:Tbound))
% title('New Death')
%
% subplot(3,3,6)
% plot(long_alpha(1:Tbound))
% title('Alpha')
%
% subplot(3,3,7)
% plot(Malt(1:Tbound))
% title('Mobility')


figure(1114)
Tbound = 30;
set(gcf,'Position',[100,100,1200,800])
subplot(4,2,1)
plot(long_delta(1:Tbound))
title('Delta (1-30 week)')

subplot(4,2,2)
plot(long_beta(1:Tbound))
hold on
plot(long_beta_tilde(1:Tbound,find(DR == DR_index(1))))
plot(long_N(1:Tbound,find(DR == DR_index(1)))./long_I(1:Tbound))
plot(POP0./long_S(1:Tbound,1))
hold off
legend('beta','beta tilde','N/I','S0/St')
lgd = legend;
lgd.Location = 'southeast';
title('Beta (1-30 week)')

subplot(4,2,3)
yyaxis left
plot(long_N(1:Tbound,find(DR == DR_index(1))))
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
plot(long_alpha(1:Tbound,find(DR == DR_index(1))))
hold on
ylabel('Alpha')
yyaxis right
plot(long_ERN(1:Tbound,find(DR == DR_index(1))))
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
plot(long_beta_tilde(Tbound+1:end,find(DR == DR_index(1))))
plot(long_N(Tbound+1:Tdata,find(DR == DR_index(1)))./long_I(Tbound+1:Tdata))
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
plot(long_N(Tbound+1:end,find(DR == DR_index(1))))
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
plot(long_alpha(Tbound+1:end,find(DR == DR_index(1))))
hold on
ylabel('Alpha')
yyaxis right
plot(long_ERN(Tbound+1:end,find(DR == DR_index(1))))
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
plot(long_ERN(1:Tbound,find(DR == DR_index(1))),'-k','LineWidth',2)
legend('delta','beta','ERN')
lgd = legend;
lgd.Location = 'northwest';
title('Delta (1-30 week)')

subplot(2,3,2)
yyaxis left
plot(long_N(1:Tbound,find(DR == DR_index(1))),'LineWidth',2)
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
plot(long_alpha(1:Tbound,find(DR == DR_index(1))),'LineWidth',2)
hold on
yyaxis right
plot(long_beta(1:Tbound),'LineWidth',2)
hold on
plot(long_beta_tilde(1:Tbound,find(DR == DR_index(1))),'-k','LineWidth',2)
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
plot(long_ERN(Tbound+1:end,find(DR == DR_index(1))),'-k','LineWidth',2)
xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
legend('delta','beta','ERN')
lgd = legend;
lgd.Location = 'northeast';
title('Delta since 31st week')

subplot(2,3,5)
yyaxis left
plot(long_N(Tbound+1:end,find(DR == DR_index(1))),'LineWidth',2)
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
plot(long_alpha(Tbound+1:end,find(DR == DR_index(1))),'LineWidth',2)
hold on
yyaxis right
plot(long_beta(Tbound+1:end),'LineWidth',2)
hold on
plot(long_beta_tilde(Tbound+1:end,find(DR == DR_index(1))),'-k','LineWidth',2)
xline(Tdata-Tbound+1,'LineWidth',1.5,'HandleVisibility','off');
legend('alpha','beta','beta tilde')
title('Alpha since 31st week')




%% figure for vaccine path
% 
% figure(200)
% set(gcf,'Position',[100,100,1200,800])
% idata = 3-3+1;
% Data_vac = VT_mat(:,:,idata);
% % medStuff1_nv = [V1_w; Data_vac(:,1)]; % In future, you should change V1_w to V1past_med.
% % medStuff2_nv = [V2_w; Data_vac(:,2)]; % In future, you should change V2_w to V2past_med.
% % elderly1_nv = [zeros(Tdata,1); Data_vac(:,3)];
% % elderly2_nv = [zeros(Tdata,1); Data_vac(:,4)];
% % others1_nv = [zeros(Tdata,1); Data_vac(:,5)];
% % others2_nv = [zeros(Tdata,1); Data_vac(:,6)];
% medStuff1_nv = zeros(Tdata+size(Data_vac,1),1);
% medStuff2_nv = zeros(Tdata+size(Data_vac,1),1);
% elderly1_nv = [zeros(Tdata,1); Data_vac(:,1)];
% elderly2_nv = [zeros(Tdata,1); Data_vac(:,2)];
% others1_nv = [zeros(Tdata,1); Data_vac(:,3)];
% others2_nv = [zeros(Tdata,1); Data_vac(:,4)];
% 
% 
% % Find the cumulative number of vaccines for each agent
% medStuff1_cv = cumsum(medStuff1_nv);
% medStuff2_cv = cumsum(medStuff2_nv);
% elderly1_cv = cumsum(elderly1_nv);
% elderly2_cv = cumsum(elderly2_nv);
% others1_cv = cumsum(others1_nv);
% others2_cv = cumsum(others2_nv);
% 
% Area_nv = [medStuff1_nv,medStuff2_nv,elderly1_nv,elderly2_nv,others1_nv,others2_nv];
% Area_cv = [medStuff1_cv,medStuff2_cv,elderly1_cv,elderly2_cv,others1_cv,others2_cv];
% Area_nv = Area_nv / (ps*10000);
% Area_cv = Area_cv / (ps*100000000);
% % Area_nv = [medStuff1_nv;medStuff2_nv;elderly1_nv;elderly2_nv;others1_nv;others2_nv];
% % Area_cv = [medStuff1_cv;medStuff2_cv;elderly1_cv;elderly2_cv;others1_cv;others2_cv];
% 
% 
% for ifig = 1:1:2
%     if ifig == 1
%         subplot(3,2,1)
%         %         figure(71)
%         Areagraph = bar(Area_nv,'stacked','LineStyle','none');
%     else
%         subplot(3,2,2)
%         %         figure(72)
%         Areagraph = area(Area_cv,'LineStyle','none');
%         %         Areagraph = area(transpose(Area_cv),'LineStyle','none');
%         Areagraph = bar(Area_cv,'stacked','LineStyle','none');
%     end
%     xtickangle(45)
%     xticks(find(WeekNumber==1))
%     xlim([Tdata-7 Tdata+56])
%     xticklabels(MonthWeekEN(WeekNumber==1))
%     
%     ldfs = fs/3;
%     if ifig == 1
%         ylim([0 800])
%         ylabel('ワクチン本数（万本）','FontName',fn)
%         %             legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3),Areagraph(2),Areagraph(1)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','医療従事者２本目','医療従事者１本目','FontSize',ldfs,'Location','NorthEast','FontName',fn);
%         legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','FontSize',ldfs,'Location','NorthEast','FontName',fn);
%         title('新規ワクチン接種本数（週ごと）', 'FontSize',fs,'FontName',fn);
%     else
%         ylim([0 2.5])
%         ylabel('ワクチン本数（億本）','FontName',fn)
%         %             legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3),Areagraph(2),Areagraph(1)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','医療従事者２本目','医療従事者１本目','FontSize',ldfs,'Location','NorthWest','FontName',fn);
%         legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','FontSize',ldfs,'Location','NorthWest','FontName',fn);
%         title('累計ワクチン接種本数（週ごと）', 'FontSize',fs,'FontName',fn);
%     end
%     
%     ax = gca;
%     ax.YAxis.Color = 'k';
%     ax.YAxis.FontSize = fs;
%     ax.XAxis.FontSize = fs;
%     Areagraph(3).FaceColor = [0.2 0.6 0.8];
%     Areagraph(1).FaceColor = [0.6 0.6 0.6];
%     Areagraph(5).FaceColor = [0.4 0.4 0.8];
%     Areagraph(4).FaceColor = [0.2 0.6 0.8];
%     Areagraph(2).FaceColor = [0.6 0.6 0.6];
%     Areagraph(6).FaceColor = [0.4 0.4 0.8];
%     Areagraph(1).FaceAlpha = 1;
%     Areagraph(3).FaceAlpha = 0.7;
%     Areagraph(5).FaceAlpha = 1;
%     Areagraph(2).FaceAlpha = 0.8;
%     Areagraph(4).FaceAlpha = 0.5;
%     Areagraph(6).FaceAlpha = 0.8;
% end

%% ----------- effect of vaccines on mortality rate  ----------------------------
% subplot(3,2,3)
% % figure(73)
% delta_init = deltaT(1,1);
% plot([delta; delta_mat(:,idata)]./delta_average,'LineWidth',1.5);
% % plot([delta; deltaT(:)]./delta_init, 'k', 'LineWidth',1.5);
% 
% xtickangle(45)
% xticks(find(WeekNumber==1))
% xlim([Tdata+1 Tdata+56])
% xticklabels(MonthWeekEN(WeekNumber==1))
% 
% title('致死率（現在のレベルで標準化）','FontName',fn)
% ax = gca;
% ax.YAxis.Color = 'k';
% ax.YAxis.FontSize = fs;
% ax.XAxis.FontSize = fs;

%% ------------ Compare mean(delta) and weighted(delta) ------------- %%
% figure(171)
% delta_sample = delta(end-RetroPeriod+1:end);
% delta_path_1 = [delta; mean(delta_sample)*ones(SimPeriod,1)];
% delta_path_2 = [delta; sum(delta_sample.*(I(end-RetroPeriod+1:end)/sum(I(end-RetroPeriod+1:end))))*ones(SimPeriod,1)];
% d(1)=plot(delta_path_1,'LineWidth',1);
% hold on
% d(2)=plot(delta_path_2,'LineWidth',1);
% ylim([0 0.1])
% xlim([50 Tdata+SimPeriod])
% xtickangle(45)
% xticks(find(WeekNumber==1))
% xticklabels(MonthWeekEN(WeekNumber==1))
% legend([d(1),d(2)],["Mean Delta", "Average Delta weighted by I"])
% title('\delta')
% xline(Tdata,'LineWidth',1,'HandleVisibility','off');
% 
% subplot(3,2,4)
% % figure(172)
% sw_delta = 31;
% delta_path_1 = zeros(Tdata+SimPeriod,1);
% delta_path_2 = zeros(Tdata+SimPeriod,1);
% for i_delta = sw_delta:1:size(delta,1)
%     delta_sample = delta(i_delta-RetroPeriod+1:i_delta);
%     delta_path_1(i_delta) = mean(delta_sample);
%     delta_path_2(i_delta) = sum(delta_sample.*(I(i_delta-RetroPeriod+1:i_delta)/sum(I(i_delta-RetroPeriod+1:i_delta))));
% end
% d(1)=plot(delta_path_1,'LineWidth',1);
% hold on
% d(2)=plot(delta_path_2,'LineWidth',1);
% ylim([0 0.05])
% xlim([sw_delta Tdata])
% xtickangle(45)
% xticks(find(WeekNumber==1))
% xticklabels(MonthWeekEN(WeekNumber==1))
% legend([d(1),d(2)],["Mean Delta", "Average Delta weighted by I"])
% title('\delta')
% xline(Tdata,'LineWidth',1,'HandleVisibility','off');
% 
% subplot(3,2,5)
% % figure(173)
% sw_delta = 25;
% for i = 11:2:23
%     delta_path_1 = zeros(Tdata+SimPeriod,1);
%     for i_delta = sw_delta:1:size(delta,1)
%         delta_sample = delta(i_delta-i+1:i_delta);
%         delta_path_1(i_delta) = mean(delta_sample);
%     end
%     if i == 17
%         d(i) = plot(delta_path_1,'-r','LineWidth',1);
%     else
%         d(i) = plot(delta_path_1,'--','LineWidth',0.5);
%     end
%     hold on
% end
% % ylim([0 0.05])
% xlim([sw_delta Tdata])
% xtickangle(45)
% xticks(find(WeekNumber==1))
% xticklabels(MonthWeekEN(WeekNumber==1))
% legend([d(11:2:23)],string([11:2:23]),'FontSize',6,'Location','NorthEast')
% title('Mean \delta')
% xline(Tdata,'LineWidth',1,'HandleVisibility','off');
% 
% 
% subplot(3,2,6)
% % figure(174)
% sw_delta = 25;
% for i = 11:2:23
%     delta_path_2 = zeros(Tdata+SimPeriod,1);
%     for i_delta = sw_delta:1:size(delta,1)
%         delta_sample = delta(i_delta-i+1:i_delta);
%         delta_path_2(i_delta) = sum(delta_sample.*(I(i_delta-i+1:i_delta)/sum(I(i_delta-i+1:i_delta))));
%     end
%     if i == 17
%         d(i) = plot(delta_path_2,'-r','LineWidth',1);
%     else
%         d(i) = plot(delta_path_2,'--','LineWidth',0.5);
%     end
%     hold on
% end
% % ylim([0 0.05])
% xlim([sw_delta Tdata])
% xtickangle(45)
% xticks(find(WeekNumber==1))
% xticklabels(MonthWeekEN(WeekNumber==1))
% legend([d(11:2:23)],string([11:2:23]),'FontSize',6,'Location','NorthEast')
% title("Average \delta weighted by I")
% xline(Tdata,'LineWidth',1,'HandleVisibility','off');


%% figure for vaccine path
% cumulative_nv = zeros(Tdata+size(VT_mat,1),3); %zeros(size([V1_w; Data_vac(:,1)],1),3);  % Tdata = 52 should be adjusted.
% % data_name = ["3","6","9","12","15","18","21","52"];
% [V,deltaT,VT] = vaccine_distribution_medical(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2,ps,POP0,3);
% ldfs = fs/3;
% for idata = 1 %1:length(VP)
%     
%     % set(gcf,'Position',[100,100,1200,800])
%     Data_vac = VT;  %VT_mat(:,:,idata);
%     %     medStuff1_nv = [V1_w; Data_vac(:,1)]; % In future, you should change V1_w to V1past_med.
%     %     medStuff2_nv = [V2_w; Data_vac(:,2)]; % In future, you should change V2_w to V2past_med.
%     %     medStuff1_nv = zeros(Tdata+size(Data_vac,1),1);
%     %     medStuff2_nv = zeros(Tdata+size(Data_vac,1),1);
%     medStuff1_nv = [V1_w; Data_vac(:,5)]; % In future, you should change V1_w to V1past_med.
%     medStuff2_nv = [V2_w; Data_vac(:,6)]; % In future, you should change V2_w to V2past_med.
%     elderly1_nv = [zeros(Tdata,1); Data_vac(:,1)];
%     elderly2_nv = [zeros(Tdata,1); Data_vac(:,2)];
%     others1_nv = [zeros(Tdata,1); Data_vac(:,3)];
%     others2_nv = [zeros(Tdata,1); Data_vac(:,4)];
%     
%     % Find the cumulative number of vaccines for each agent
%     medStuff1_cv = cumsum(medStuff1_nv);
%     medStuff2_cv = cumsum(medStuff2_nv);
%     elderly1_cv = cumsum(elderly1_nv);
%     elderly2_cv = cumsum(elderly2_nv);
%     others1_cv = cumsum(others1_nv);
%     others2_cv = cumsum(others2_nv);
%     
%     Area_nv = [medStuff1_nv,medStuff2_nv,elderly1_nv,elderly2_nv,others1_nv,others2_nv];
%     Area_cv = [medStuff1_cv,medStuff2_cv,elderly1_cv,elderly2_cv,others1_cv,others2_cv];
%     Area_nv = Area_nv / (ps*10000);
%     Area_cv = Area_cv / (ps*100000000);
%     % Area_nv = [medStuff1_nv;medStuff2_nv;elderly1_nv;elderly2_nv;others1_nv;others2_nv];
%     % Area_cv = [medStuff1_cv;medStuff2_cv;elderly1_cv;elderly2_cv;others1_cv;others2_cv];
%     
%     cumulative_nv(:,idata) = medStuff1_nv + medStuff2_nv + elderly1_nv + elderly2_nv + others1_nv + others2_nv;
%     
%     for ifig = 1:1:2
%         if ifig == 1
%             %         subplot(3,2,1)
%             %         figure(71)
%             figure(idata+300)
%             %             Areagraph = area(Area_nv,'LineStyle','none');
%             Areagraph = bar(Area_nv,'stacked','LineStyle','none');
%         else
%             %         subplot(3,2,2)
%             %         figure(72)
%             figure(idata+3300)
%             %             Areagraph = area(Area_cv,'LineStyle','none');
%             Areagraph = bar(Area_cv,'stacked','LineStyle','none');
%         end
%         xtickangle(45)
%         xticks(find(WeekNumber==1))
%         xlim([Tdata-7 Tdata+44])
%         xticklabels(MonthWeekJP(WeekNumber==1))
%         
%         if ifig == 1
%             ylim([0 800])
%             ylabel('ワクチン本数（万本）','FontName',fn)
%                         legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3),Areagraph(2),Areagraph(1)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','医療従事者２本目','医療従事者１本目','FontSize',ldfs,'Location','NorthEast','FontName',fn);
% %             legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','FontSize',ldfs,'Location','NorthEast','FontName',fn);
%             title('新規ワクチン接種本数（週ごと）', 'FontSize',fs,'FontName',fn);
%         else
%             ylim([0 2.5])
%             ylabel('ワクチン本数（億本）','FontName',fn)
%                         legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3),Areagraph(2),Areagraph(1)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','医療従事者２本目','医療従事者１本目','FontSize',ldfs,'Location','NorthWest','FontName',fn);
% %             legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','FontSize',ldfs,'Location','NorthWest','FontName',fn);
%             title('累計ワクチン接種本数（週ごと）', 'FontSize',fs,'FontName',fn);
%         end
%         
%         ax = gca;
%         %         ax.YAxis.Color = 'k';
%         ax.YAxis.FontSize = fs;
%         ax.XAxis.FontSize = fs;
%         ax.YAxis.FontName = fn;
%         ax.XAxis.FontName = fn;
%         Areagraph(3).FaceColor = [0.2 0.6 0.8];
%         Areagraph(1).FaceColor = [0.6 0.6 0.6];
%         Areagraph(5).FaceColor = [0.4 0.4 0.8];
%         Areagraph(4).FaceColor = [0.2 0.6 0.8];
%         Areagraph(2).FaceColor = [0.6 0.6 0.6];
%         Areagraph(6).FaceColor = [0.4 0.4 0.8];
%         Areagraph(1).FaceAlpha = 1;
%         Areagraph(3).FaceAlpha = 0.7;
%         Areagraph(5).FaceAlpha = 1;
%         Areagraph(2).FaceAlpha = 0.8;
%         Areagraph(4).FaceAlpha = 0.5;
%         Areagraph(6).FaceAlpha = 0.8;
%     end
%     %
%     %     figure(idata+20)
%     %     if idata == 1
%     %         vplot(1) = plot(cumulative_nv(:,1)/10000, 'k', 'LineWidth',2);
%     %     elseif idata == 2
%     %         vplot(2) = plot(cumulative_nv(:,2)/10000, 'b', 'LineWidth',2);
%     %         hold on
%     %         vplot(1) = plot(cumulative_nv(:,1)/10000, 'k', 'LineWidth',2);
%     %     else
%     %         vplot(2) = plot(cumulative_nv(:,2)/10000, 'b', 'LineWidth',2);
%     %         hold on
%     %         vplot(3) = plot(cumulative_nv(:,3)/10000, 'r', 'LineWidth',2);
%     %         vplot(1) = plot(cumulative_nv(:,1)/10000, 'k', 'LineWidth',2);
%     %     end
%     %     hold off
%     %     ylim([0 800])
%     %     ylabel('ワクチン本数（万本）')
%     %
%     %     xtickangle(45)
%     %     xticks(find(WeekNumber==1))
%     %     xlim([Tdata-7 Tdata+44])
%     %     xticklabels(MonthWeekJP(WeekNumber==1))
%     %     if idata == 1
%     %         legend([vplot(1)], '基本見通し','FontSize',ldfs,'Location','NorthEast','FontName',fn);
%     %     elseif idata == 2
%     %         legend([vplot(2),vplot(1)],'１２週間間隔','基本見通し','FontSize',ldfs,'Location','NorthEast','FontName',fn);
%     %     elseif idata == 3
%     %         legend([vplot(2),vplot(1),vplot(3)],'１２週間間隔','基本見通し','一本しか打たない', 'FontSize',ldfs,'Location','NorthEast','FontName',fn);
%     %     end
%     %     title('ワクチン接種','FontName',fn);
%     %     ax = gca;
%     %     ax.YAxis.Color = 'k';
%     %     ax.YAxis.FontSize = fs;
%     %     ax.XAxis.FontSize = fs;
%     %     ax.YAxis.FontName = fn;
%     %     ax.XAxis.FontName = fn;
%     %
%     if figure_save == 1
%         if iPC == 1
%             saveas(figure(idata+300),char(strcat(home,'Figures\',num2str(VP(idata)),'_NV.png')));
%             saveas(figure(idata+3300),char(strcat(home,'Figures\',num2str(VP(idata)),'_CV.png')));
%             %             saveas(figure(idata+20),char(strcat(home,'Figures\',data_name(idata),'_VShots.png')));
%         else
%             saveas(figure(idata+300),char(strcat(home,'/shotaro/',num2str(VP(idata)),'_NV.png')));
%             saveas(figure(idata+3300),char(strcat(home,'/shotaro/',num2str(VP(idata)),'_CV.png')));
%             %             saveas(figure(idata+20),char(strcat(home,'/shotaro/',data_name(idata),'_VShots.png')));
%         end
%     end
%     
% end


%% ----------- effect of vaccines on mortality rate  ----------------------------
% figure(55)
% % for ivac = 1:length(VP)
% %     if VP(ivac) == VP_index(1)
% %         plot([delta; delta_mat(:,ivac)]./delta_average,'-r','LineWidth',1.5,'DisplayName',sprintf('%.0f',VP(ivac)));
% %     elseif VP(ivac) == VP_index(2)
% %         plot([delta; delta_mat(:,ivac)]./delta_average,'-b','LineWidth',1.5,'DisplayName',sprintf('%.0f',VP(ivac)));
% %     else
% %         plot([delta; delta_mat(:,ivac)]./delta_average,'--','LineWidth',0.5,'DisplayName',sprintf('%.0f',VP(ivac)));
% %     end
% %     hold on
% % end
% plot([delta; deltaT]./delta_average,'-r','LineWidth',1.5);
% 
% % dplot(6) = plot([delta; delta_mat(:,1)]./delta_average,'LineWidth',1.5,'DisplayName',num2str(52));
% % hold on
% % for i  = 1:1:7
% %     dplot(i) = plot([delta; delta_mat(:,i)]./delta_average,'LineWidth',1.5,'DisplayName',num2str(5*i));
% % end
% % dplot(i+1) = plot([delta; delta_mat(:,i+1)]./delta_average,'LineWidth',1.5,'DisplayName',num2str(52));
% % dplot(3) = plot([delta; delta_mat(:,2)]./delta_average, 'b', 'LineWidth',1.5); % optimistic
% % hold on
% % dplot(2) = plot([delta; delta_mat(:,3)]./delta_average, 'r', 'LineWidth',1.5); % pessimistic
% % dplot(1) = plot([delta; delta_mat(:,1)]./delta_average, 'k', 'LineWidth',1.5); % baseline
% 
% 
% xtickangle(45)
% xticks(find(WeekNumber==1))
% xlim([Tdata+1 Tdata+44])
% xticklabels(MonthWeekJP(WeekNumber==1))
% % lgd = legend;
% title('致死率（現在のレベルで標準化）','FontName',fn);
% ax = gca;
% ax.YAxis.Color = 'k';
% ax.YAxis.FontSize = fs;
% ax.XAxis.FontSize = fs;
% ax.YAxis.FontName = fn;
% ax.XAxis.FontName = fn;
% 
% if figure_save == 1
%     if iPC == 1
%         saveas(figure(55),char(strcat(home,'Figures\death_rate.png')));
%     end
% end


