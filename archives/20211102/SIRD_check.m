% Check SIRD (model and data)
%
clear variables
close all
iPC = 0; % 0 for Mac, 1 for Windows
if iPC == 1
%     home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
    home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/';
end
cd(home);
%====================== Program parameter values ======================%
vaccine_disp_switch = 0;
fs = 12; % common font size for many figures
ldfs = 12; % legend font size for vaccine path
ft = '%.1f';
if iPC == 1
    fn = 'Yu Gothic'; % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
linecolor = {'red', 'blue'};
language = {'EN','JP'};
%======================================================================%

%================== Model Fixed Parameter Values ============================%
SimPeriod = 52;        % simulation period in weeks
gamma = 7/12;          % recovery rate from Covid % Should change this to 7/12 (4/25 Kohei Machi)
k = 2;                 % exponent of (1-h*alpha)
hconstant = 1;         % 0 = without intercept, 1 = with intercept for h regression
tt = 12; % Showing previous t periods for the plot
DRi=26;%17;%26 %10; % 経済回復速度
% RetroPeriods for sample mean
retro_ub = 17; % Control the moving average of beta (beta_avg = sum_{t = lb}^{ub} (1/(ub-lb + 1) sum_{x=1}^t (1/t) beta_t)
retro_lb = 17;
RetroPeriod = 17;      % retroactive periods used to estimate gamma 
RetroPeriodDelta = 4; %17;      % retroactive periods used to estimate delta (Death rate & Tokyo standard Severe case)
RetroPeriodICU_nation = 13; %17;      % retroactive periods used to estimate delta ICU (National, Hospital)
RetroPeriodICU_pref = 13;
RetroPeriodHospital = 13;
retroH_switch = 1; %If retroH_switch == 1, retroH = TdataGDP - 4, else = retroH
RetroH = 15;
% Parameters for variants
var_infection = 0.3; %Relative increase of infectiousness (alpha varaint)
var_infection_delta = 0.4; %Relative increase of death rate (alpha varaint)
var_growth = 0.47;
% Population
POP_jp = 125710000;
medical_jp = 4700000;
elderly_jp = 36000000;
ordinary_jp = (POP_jp-elderly_jp-medical_jp);
accept_share = 0.9;
accept_share_ordinary = 0.6879; %so that an accept share of age 13-64 = 80%
% vaccine pace 
E1 = 0.45; %(0707)%(0609), %0.38
E2 = 0.815; %(0707)%(0609), %0.795
D1 = 0.865; %(0707)%(0609), %0.68
D2 = 0.96; %(0707)%(0609), %0.925
% parameters for ICU
gamma_ICU_nation = 7/28; % Recovery rate from ICU
gamma_ICU_pref = 7/28; % Recovery rate from ICU
ICU_nation_adjustment = 0.8;
ICU_pref_adjustment = 0.8;
% parameters for Hospitalizaiton
gamma_Hospital = 7/15; %7/10;
Hospital_adjustment = 0.85; 
Hospital_limit = 6406; %5882;

%================ Parameter Values (Prefecture Specific) ====================%
pref = 'Tokyo';        % prefecture to be analyzed
prefGDP = 106; % 兆円, one trillion yen (chou-yen)
covid = importdata([home 'Covid_weekly_newV.csv']);  % Import weekly Covid data by prefecture
Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
Data(isnan(Data)) = 0;
% Columns: 1 = date, 2 = new positive, 3 = death, 4 = mobility,
% 5 = GDP, 6 = population, 7 = GDO
dateD = Data(:,1) + 21916;
Tdata= size(dateD,1);    % Data period in weeks
N = Data(:,2);
dD = Data(:,3);
M = Data(:,4);
% GDP = Data(:,5);
POP = Data(:,6);
GDP = Data(:,7);
I_data = zeros(Tdata+1,1);
I_data(2:end) = Data(:,27);
R_data = zeros(Tdata+1,1);
R_data(2:end)= Data(:,28);
D_data = zeros(Tdata+1,1);
D_data(2:end) = Data(:,29);
S_data = [POP(1);POP] - (I_data + R_data + D_data);
dI_data = I_data(2:end) - I_data(1:end-1);
dR_data = R_data(2:end) - R_data(1:end-1);
dD_data = D_data(2:end) - D_data(1:end-1);
dS_data = S_data(2:end) - S_data(1:end-1);
N_pref = Data(:,26);
Severe_pref = zeros(Tdata+1,1);
Severe_pref(2:end) = Data(:,30);
hospital = zeros(Tdata+1,1);
hospital(2:end) = Data(:,25);
hospital(isnan(hospital)) = 0;
%--- Import ICU data ---%
ICU_nation = zeros(Tdata+1,1);
ICU_pref = zeros(Tdata+1,1);
BED = zeros(Tdata,1);
ICU_nation(2:Tdata+1,1) = Data(:,22);
BED(1:Tdata,1) = Data(:,23);
ICU_pref(2:Tdata+1,1) = Data(:,21);

GDP = GDP/GDP(1)*100;   % Normalization: 100 = 2020JAN

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
% elderly_start = find(SimDateEN == elderly_start_date);
% VacStart = find(SimDateEN == datetime(2021,4,1));
% End2020 = find(dateEN == datetime(2021,1,7));

%ICU
ICU_limit = 392;

%--- Construct weekly vaccine data ---%
[V1_medical_ori, V1_medical, V2_medical_ori, V2_medical,...
    V1_elderly_ori, V1_elderly, V2_elderly_ori, V2_elderly, ...
    V1_others_ori, V1_others, V2_others_ori, V2_others,...
    vs_MT,vs_ET,vs_M1,vs_E1,vs_M2,vs_E2] ...
    = ImportVaccineData(home,iPC,Data,pref,dateEN,ps,vaccine_disp_switch);

%--- Constructing the reference level of output ---%
[potentialGDP, referenceGDP, alpha] = construct_GDP(GDP,TdataGDP);

%--- Regress mobility on alpha to estimate the elasticity h ---%
[Malt,h_all,h_all_se,h,h_se] = estimate_h(M,alpha,TdataGDP,RetroH,hconstant);

%--- Compute the history of S, I, R, D in the data period ---%
[S,I,R,D]...
    = SIRD(Tdata,POP0,N,E1,E2,...
    V1_elderly,V1_medical,V1_others,V2_elderly,V2_medical,V2_others,...
    gamma,dD,TdataGDP,referenceGDP,alpha);
InitialValues = [S(end),I(end),R(end),D(end),ICU_nation(end),ICU_pref(end),hospital(end)];

gamma_data = zeros(Tdata,1);
delta_data = zeros(Tdata,1);
for i = 1:Tdata
    if I(i)>0
        gamma_data(i) = dR_data(i)/I_data(i);
        delta_data(i) = dD_data(i)/I_data(i);
    end
end
gamma_data_left = 24; %02-Jul-2020
gamma_data_right = Tdata; %04-Feb-2021
gamma_sample_data = gamma_data(gamma_data_left:gamma_data_right);
% delta_average_data = ...
%     sum(delta_sample_data.*(I_data(delta_data_left:delta_data_right)/sum(I_data(delta_data_left:delta_data_right))));
gamma_average_data = mean(gamma_sample_data);


% delta_sample_data = delta_data(end-RetroPeriodDelta+1:end);
delta_data_left = 24; %02-Jul-2020
delta_data_right = 55; %04-Feb-2021
delta_sample_data = delta_data(delta_data_left:delta_data_right);
% delta_average_data = ...
%     sum(delta_sample_data.*(I_data(delta_data_left:delta_data_right)/sum(I_data(delta_data_left:delta_data_right))));
delta_average_data = mean(delta_sample_data);

gamma_path = gamma*ones(Tdata,1);
% ICU_flow_rate_nation = zeros(Tdata,1);
% for t = 1:Tdata
%     if I(t)>0
%         ICU_flow_rate_nation(t) = (ICU_nation(t+1) - ICU_nation(t) + dD(t))/I(t);
%     end
% end

[delta,beta_tilde,ERN,beta,ICU_nation_inflow, ICU_pref_inflow, Hospital_inflow,...
    gammaT,delta_average,delta_ICU_nation_average,delta_ICU_pref_average,delta_Hospital_average,...
    ICU_nation_inflow_avg,ICU_pref_inflow_avg,Hospital_inflow_avg,beta_avg,beta_se,delta_se]...
    = Time_Series_Average(S,I,D,ICU_nation,ICU_pref,hospital,dD,N,Tdata,SimPeriod,...
        RetroPeriod,POP0,hconstant,h_all,alpha,k,...
        gamma,gamma_ICU_nation,gamma_ICU_pref,gamma_Hospital,...
        ICU_nation_adjustment,ICU_pref_adjustment,Hospital_adjustment,...
        RetroPeriodDelta,RetroPeriodICU_nation,RetroPeriodICU_pref,RetroPeriodHospital);

ICU_nation_predicted=zeros(Tdata+1,1);
ICU_pref_predicted=zeros(Tdata+1,1);
ICU_nation_predicted(2:end) = ICU_nation_predicted(1:end-1) + delta_ICU_nation_average.*I(1:end-1)*ICU_nation_inflow_avg - gamma_ICU_nation*ICU_nation_predicted(1:end-1) - delta_average.*I(1:end-1);
ICU_pref_predicted(2:end) = ICU_pref_predicted(1:end-1) + delta_ICU_pref_average.*I(1:end-1)*ICU_pref_inflow_avg - gamma_ICU_pref*ICU_pref_predicted(1:end-1) - delta_average.*I(1:end-1);

%% =========== Plot  =========== %%
format shortG
% Figure parameter
% ==========================%
lineWidth = [1.5,1.5];
lineName = {'データ', 'モデル'};
linecolor = {'black', 'red'};
% ==========================%


% Transitions of I (data and model)
figname_NI = 'I_past_transitions';
Title = 'Iの過去の推移';
yvalue = [I_data, I];
yft = '%.0f';
yft2 = '%.2f';

figure('Name',char(figname_NI));
set(gcf,'Position',[100,100,1000,600])
yyaxis left
for i = 1:2
    plot(yvalue(2:end,i),linecolor{i},'lineStyle','-','LineWidth',lineWidth(i),'DisplayName',lineName{i})
    hold on
end
ytickformat(yft)
ylabel('I: 療養者数')
yyaxis right
plot(gamma_data,'--k','LineWidth',1.0,'DisplayName','\gamma (データ: 退院者/療養者)')
hold on
plot(gamma_path,'--r','LineWidth',1.0,'DisplayName','\gamma (モデルの仮定) = 7/12')
hold on 
plot(mean(gamma_data(end-RetroPeriod+1:end))*ones(Tdata,1),'--b','LineWidth',1.0,'DisplayName',['\gamma (データ, 17週平均)=',sprintf(yft2,mean(gamma_data(end-RetroPeriod+1:end)))])
xline(Tdata-RetroPeriod+1,'LineWidth',1.5,'HandleVisibility','off');
ytickformat(yft2)
xticks(find(WeekNumber==1))
xticklabels(MonthWeekJP(WeekNumber==1))
ylabel('\gamma: 回復率')
xlabel('週')
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'b';
ax.YAxis(1).FontSize = fs; ax.YAxis(2).FontSize = fs; ax.XAxis.FontSize = fs;
ax.YAxis(1).Exponent = 0; ax.YAxis(2).Exponent = 0;
title(char(Title),'FontSize',fs,'FontWeight','normal','FontName',fn)
lgd = legend; lgd.FontSize = ldfs; lgd.Location = 'NorthWest';
xtickangle(45)

% Transitions of S (data and model)
figname = 'S_past_transitions';
Title = 'Sの過去の推移';
yvalue = [S_data, S];
yft = '%.0f';

figure('Name',char(figname));
set(gcf,'Position',[100,100,1000,600])

for i = 1:2
    plot(yvalue(2:end,i),linecolor{i},'lineStyle','-','LineWidth',lineWidth(i),'DisplayName',lineName{i})
    hold on
end
ax = gca; ax.YAxis.FontSize = ldfs; ax.XAxis.FontSize = ldfs; ax.YAxis.Exponent = 0;
ytickformat(yft)
xticks(find(WeekNumber==1))
title(char(Title),'FontSize',fs,'FontWeight','normal','FontName',fn)
xticklabels(MonthWeekJP(WeekNumber==1))
lgd = legend; lgd.FontSize = ldfs; % lgd.Location = lgdLocation;
xtickangle(45)

% Transitions of R (data and model)
figname = 'R_past_transitions';
Title = 'Rの過去の推移';
yvalue = [R_data, R];
yft = '%.0f';

figure('Name',char(figname));
set(gcf,'Position',[100,100,1000,600])

for i = 1:2
    plot(yvalue(2:end,i),linecolor{i},'lineStyle','-','LineWidth',lineWidth(i),'DisplayName',lineName{i})
    hold on
end
ax = gca; ax.YAxis.FontSize = ldfs; ax.XAxis.FontSize = ldfs; ax.YAxis.Exponent = 0;
ytickformat(yft)
xticks(find(WeekNumber==1))
title(char(Title),'FontSize',fs,'FontWeight','normal','FontName',fn)
xticklabels(MonthWeekJP(WeekNumber==1))
lgd = legend; lgd.FontSize = ldfs; % lgd.Location = lgdLocation;
xtickangle(45)

% Transitions of dD (data and model)
figname = 'dD_past_transitions';
Title = 'dDの過去の推移';
yft = '%.0f';

figure('Name',char(figname));
set(gcf,'Position',[100,100,1000,600])
yyaxis left
plot(dD_data,linecolor{1},'lineStyle','-','LineWidth',lineWidth(1),'DisplayName',lineName{1})
hold on
plot(I_data(1:end-1)*delta_average_data,'black','lineStyle',':','LineWidth',1.5,'DisplayName','Predicted (data)')
plot(I(1:end-1)*delta_average,'red','lineStyle',':','LineWidth',1.5,'DisplayName','Predicted (model)')
ytickformat(yft)
yyaxis right
plot(delta_data,'--k','LineWidth',1.0,'DisplayName','delta (data)')
hold on
plot(delta,'--r','LineWidth',1.0,'DisplayName','delta (model)')
ytickformat(yft2)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'b';
ax.YAxis(1).FontSize = fs; ax.YAxis(2).FontSize = fs; ax.XAxis.FontSize = fs;
xticks(find(WeekNumber==1))
title(char(Title),'FontSize',fs,'FontWeight','normal','FontName',fn)
xticklabels(MonthWeekJP(WeekNumber==1))
lgd = legend; lgd.FontSize = ldfs; % lgd.Location = lgdLocation;
xtickangle(45)

% Transitions of dD and I (data and model)
figname = 'dD_I_past_transitions';
Title = 'dDとIの過去の推移';
yft = '%.0f';

figure('Name',char(figname));
set(gcf,'Position',[100,100,1000,600])
yyaxis left
plot(dD_data,'black','lineStyle','-','LineWidth',lineWidth(1),'DisplayName','新規死亡者数')
hold on
plot(ICU_pref(2:end),'blue','lineStyle','-','LineWidth',1.0,'DisplayName','重症者(都基準,data)')
plot(ICU_pref_predicted(2:end),'blue','lineStyle','--','LineWidth',1.0,'DisplayName','重症者(都基準,model)')
hold on
plot(ICU_nation(2:end),'green','lineStyle','-','LineWidth',1.0,'DisplayName','重症者(国基準,data)')
plot(ICU_nation_predicted(2:end),'green','lineStyle','--','LineWidth',1.0,'DisplayName','重症者(国基準,model)')
yyaxis right
plot(I_data(2:end),'red','lineStyle','-','LineWidth',1.5,'DisplayName','療養者(data)')
plot(I(2:end),'red','lineStyle','--','LineWidth',1.5,'DisplayName','療養者(model)')
ytickformat(yft)
ytickformat(yft2)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'b';
ax.YAxis(1).FontSize = fs; ax.YAxis(2).FontSize = fs; ax.XAxis.FontSize = fs;
xticks(find(WeekNumber==1))
title(char(Title),'FontSize',fs,'FontWeight','normal','FontName',fn)
xticklabels(MonthWeekJP(WeekNumber==1))
lgd = legend; lgd.FontSize = ldfs; lgd.Location = 'NorthWest';
xtickangle(45)

% Transitions of I and ICU (data and model)
figname = 'I_ICU_past_transitions';
Title = '感染者数と重症患者数の過去の推移';
yft = '%.0f';
figure('Name',char(figname));
set(gcf,'Position',[100,100,1000,600])
yyaxis left
plot(I_data(2:end),'black','lineStyle','-','LineWidth',1.5,'DisplayName','療養者(data)')
hold on 
plot(hospital(2:end),'green','lineStyle','-','LineWidth',1.5,'DisplayName','入院患者者(data)')
hold on 
plot(ICU_pref(2:end),'red','lineStyle','-','LineWidth',1.0,'DisplayName','重症者(都基準,data)')
hold on
plot(ICU_nation(2:end),'blue','lineStyle','-','LineWidth',1.0,'DisplayName','重症者(国基準,data)')
ylim([0 10000])
yyaxis right
plot(hospital(2:end)./I_data(2:end),'green','lineStyle','--','LineWidth',1.5,'DisplayName','入院割合(data)')
hold on 
plot(ICU_pref(2:end)./I_data(2:end),'red','lineStyle','--','LineWidth',1.0,'DisplayName','重症者割合(都基準,data)')
hold on
plot(ICU_nation(2:end)./I_data(2:end),'blue','lineStyle','--','LineWidth',1.0,'DisplayName','重症者割合(国基準,data)')
ytickformat(yft)
ytickformat(yft2)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'b';
ax.YAxis(1).FontSize = fs; ax.YAxis(2).FontSize = fs; ax.XAxis.FontSize = fs;
xticks(find(WeekNumber==1))
title(char(Title),'FontSize',fs,'FontWeight','normal','FontName',fn)
xticklabels(MonthWeekJP(WeekNumber==1))
lgd = legend; lgd.FontSize = ldfs; lgd.Location = 'NorthWest';
xtickangle(45)
% plot(ICU_pref_predicted(2:end),'blue','lineStyle','--','LineWidth',1.0,'DisplayName','重症者(都基準,model)')
% plot(ICU_nation_predicted(2:end),'green','lineStyle','--','LineWidth',1.0,'DisplayName','重症者(国基準,model)')
% plot(I(2:end),'red','lineStyle','--','LineWidth',1.5,'DisplayName','療養者(model)')

