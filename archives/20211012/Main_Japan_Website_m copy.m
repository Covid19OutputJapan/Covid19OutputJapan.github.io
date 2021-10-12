% This m-file executes simulation and generates figures for the main
% analysis of "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
iPC= 0;
%home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Covid19OutputJapan.github.io/archives/20210525/';
if iPC == 1
    %home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
%     home = '\Users\kenic\Dropbox\fujii_nakata\Website\Codes\';
else
    % home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/';
%     home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/';
    home = '/Users/shotaro/Dropbox/fujii_nakata/Website/Codes/';
    %home = '/Users/koheimachi/Dropbox/fujii_nakata/Website/Codes/';
end
cd(home);

%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
mat_save = 1;    % 0 = figures won't be saved, 1 = they will be saved　in the "Figure" folder

data_switch = 1;
biweekly = 1;

fs = 16;            %common font size for many figures
xlim_tradeoff = [1,2.5];
% iDrawUB = 1;          %1 = create UB with simulations
%for iDrawUB = 0; we get error in line 924.
Nsim = 30000;         % if iDrawUB=1, this is the number of draws you use.
if iPC == 1
    fn = 'Yu Gothic';
else
    fn = 'YuGothic';
end
ICU_nation = 1; % = 1 use national definition (NHK data), = 0 use data from Tokyo Keizai
%======================================================================%

%====================== Model parameter values ======================%
pref = 'Japan';        % prefecture to be analyzed
POP0 = 125710000;      % initial population
parameter
SimPeriod = 52; %precaution
gamma = 0.770883744239911;
gamma_ICU = gamma_ICU_nation;
ICU_adjustment = ICU_nation_adjustment;
wl = [1,2];            % Results as of these weeks ago

prefecture_parameter

elderly_total = elderly_jp;
medical_total = medical_jp;
ordinary_total = ordinary_jp;
%medical = medical_total*accept_share;
medical = medical_total;
elderly = elderly_total*accept_share;
ordinary = ordinary_total*accept_share_ordinary;

%====================================================================%


%--- Import data ---%
% Covid data are recorded at weekly frequencies (Mon-Sun)
% The first week start on January 20 (Mon), 2020
covid = importdata([home 'Covid_weekly_newV.csv']);  % Import weekly Covid data
Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
date = Data(:,1) + 21916;
date1 =datetime(date,'ConvertFrom','excel','Format','MMM-yy');
date = datetime(date,'ConvertFrom','excel');

dateD = Data(:,1) + 21916;

N = Data(:,2);
dD = Data(:,3);
M = Data(:,4);
GDP = Data(:,5);
Tdata = length(N);
ICU = zeros(Tdata+1,1);
BED = zeros(Tdata,1);%

I_data = zeros(Tdata+1,1);
I_data(2:end) = Data(:,26);
R_data = zeros(Tdata+1,1);
R_data(2:end)= Data(:,27);
D_data = zeros(Tdata+1,1);
D_data(2:end) = Data(:,28);
S_data = ones(Tdata+1, 1) * POP0 - (I_data + R_data + D_data);
dI_data = I_data(2:end) - I_data(1:end-1);
dR_data = R_data(2:end) - R_data(1:end-1);
dD_data = D_data(2:end) - D_data(1:end-1);
dS_data = S_data(2:end) - S_data(1:end-1);

if data_switch == 1
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

%if ICU_nation == 1
    ICU(2:Tdata+1,1) = Data(:,22);
    BED(1:Tdata,1) = Data(:,23);
%else
%    ICU(2:Tdata+1,1) = Data(:,21);
%end
ps = POP0/125710000;
xtick1 = 1:13:Tdata;
%dateEN = datestr(date,'mmm-dd');
dateEN = datetime(date);
SimDate = date(end)+7:7:date(end)+7*SimPeriod;
%SimDateEN = datestr(SimDate);
SimDateEN = datetime(SimDate);
M = 1+0.01*M;
TdataGDP = Tdata-sum(isnan(GDP));
RetroH = TdataGDP-4;
%     RetroH = 15;
if isempty(find(SimDateEN == medical_start_date,1)) == 0
    medical_start = find(SimDateEN == medical_start_date);
else
    medical_start = 1;
end
elderly_start = find(SimDateEN == elderly_start_date);
VacStart = find(SimDateEN == datetime(2021,4,1));
End2020 = find(dateEN == datetime(2021,1,7));
Month = string(date1);

pindex = 1;
% Parameters for Vaccine Path

paces_ori = paces_ori_vec(pindex);
gradual_paces = gradual_paces_vec(pindex);
sw_vacpath = sw_vacpath_vec(pindex);
VT3share = VT3share_vec(pindex);
lag = lag_vec(pindex);
medical_duration = medical_duration_vec(pindex);
ind_date = find(date == ind_date_vec(pindex));

paces2_ori = paces2_ori_vec(pindex);
paces2 = paces2_ori*ps;
gradual_paces2 = gradual_paces2_vec(pindex);
ind_date2 = find(date == datetime(2021,6,24)); %when 職域接種 starts
date_slowdown = find(date == datetime(2021,8,26)) - Tdata; % when the vaccine paces slow down

paces3_ori = paces3_ori_vec(pindex);
paces3 = paces3_ori*ps;

%--- Construct weekly vaccine data ---%
vaccine_medical = readmatrix([home 'vaccine_daily_medical.xls']);
vaccine_elderly = readmatrix([home 'vaccine_daily_elderly.xls']);
[V1_medical, V2_medical] = vaccine_daily_to_weekly_table(vaccine_medical, ps, dateEN,iPC);
V1_medical(end) = V1_medical(end) * 7/4;
V2_medical(end) = V2_medical(end) * 7/4;
[V1_elderly, V2_elderly] = vaccine_daily_to_weekly_table(vaccine_elderly, ps, dateEN,iPC);
V1_elderly(end) = V1_elderly(end) * 7/4;
V2_elderly(end) = V2_elderly(end) * 7/4;
vaccine_others = readmatrix([home 'vaccine_daily_others.xls']);
[V1_others, V2_others] = vaccine_daily_to_weekly_table(vaccine_others, ps, dateEN,iPC);
V1_others(end) = V1_others(end) * 7/4;
V2_others(end) = V2_others(end) * 7/4;
% vaccine pace
Vsimple = 0; % 0 for vaccine_distribution; 1 for vaccine_distribution_simple
PF = 1; % 0 for AZ, 1 for PF
Vgradual = 1; % 0 for flat Vpath, 1 for gradually increasing Vpath
if Vsimple == 0
    VP = 3:1:17; %[0];
    VP_index = [3,8]; %[0,0];
else
    VP = 3; %[0];
    VP_index = [3,8]; %[0,0];
end
if PF == 0  % AZ
    E1 = 0.365; %0.615;
    E2 = 0.625; %0.64;
    D1 = 0.675; %0.8;
    D2 = 0.905; %0.85;
else   % PF
    E1 = 0.38; %0.625;
    E2 = 0.795; %0.895;
    D1 = 0.68; %0.8;
    D2 = 0.925; %0.94;
end
%--- Update forecast error data ---%
FE = load([home 'Forecast.mat']);

if biweekly == 1
   dDActual = [FE.dDActual;dD(end-1:end)];
   dNActual = [FE.dNActual;N(end-1:end)];
   ICUActual = [FE.ICUActual;floor(ICU(end-1:end))];
   if dNActual(end) ~= dNActual(end-1) && length(dNActual) == Tdata
        %     save([home 'Forecast.mat'],'dDActual','dNActual','-append');
        save([home 'Forecast.mat'],'dDActual','dNActual','ICUActual','-append');
   end
end

if biweekly == 0
    dDActual = [FE.dDActual;dD(end)];
    dNActual = [FE.dNActual;N(end)];
    ICUActual = [FE.ICUActual;floor(ICU(end))];
    if dNActual(end) ~= dNActual(end-1) && length(dNActual) == Tdata
        %     save([home 'Forecast.mat'],'dDActual','dNActual','-append');
        save([home 'Forecast.mat'],'dDActual','dNActual','ICUActual','-append');
    end
end

%--- Constructing the reference level of output ---%
potentialGDP = zeros(52*3,1);       % potential GDP for the next 3 years

potentialGDP(1) = (548182/(1.0122))*(1.0063^(1/12));

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


%--- Impute alpha (regress alpha on M)---%
Malt=M;
Malt(50)=0.5*(Malt(49)+Malt(51));
alpha = (1 - GDP(1:TdataGDP)./referenceGDP(1:TdataGDP));   % output loss in percentage
X = Malt(TdataGDP-17:TdataGDP);
Y = alpha(TdataGDP-17:TdataGDP);
XC = [ones(length(X),1), X];
s = (XC'*XC)\XC'*Y;         % OLS estimate of h with constant
reg = XC*s;
r = Y - reg;
SSE = sum(r.^2);
eps_p = zeros(Tdata-TdataGDP,1);
eps_p(1) = r(end);
for i = 1:Tdata-TdataGDP-1
    eps_p(i+1) = 1*eps_p(i);
end
alpha_pred = s(1)+s(2)*Malt(TdataGDP+1:Tdata)+eps_p;

alpha = [alpha;alpha_pred];


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
% [S,I,R,D,GDP]...
%     = SIRD(Tdata,POP0,N,E1,E2,V1_elderly,V1_medical,V1_others,V2_elderly,V2_medical,V2_others,gamma,dD,TdataGDP,referenceGDP,alpha,GDP);

if data_switch == 0
    [S, I, R, D, GDP] ...
        = SIRD(Tdata, POP0, N, E1, E2, ...
        V1_elderly, V1_medical, V1_others, V2_elderly, V2_medical, V2_others, ...
        gamma, dD, TdataGDP, referenceGDP, alpha);
else
    pastV = zeros(Tdata, 1);
    pastV(3:end) = E1 * (V1_elderly(1:end - 2) + V1_medical(1:end - 2) + V1_others(1:end - 2)) + (E2 - E1) * (V2_elderly(1:end - 2) + V2_medical(1:end - 2) + V2_others(1:end - 2));
    S = S_data;
    S(2:end) = S(2:end) - cumsum(pastV);
    R = R_data;
    R(2:end) = R(2:end) + cumsum(pastV);
    D = D_data;
end

[delta,beta_tilde,ERN,beta,ICU_inflow,...
    gammaT,delta_average,delta_ICU_average,ICU_inflow_avg,delta_sample,beta_avg]...
    = Time_Series_Average_Japan(S,I,D,ICU,dD,N,Tdata,SimPeriod,...
    RetroPeriod,POP0,gamma,hconstant,h_all,alpha,k,...
    gamma_ICU,ICU_adjustment,RetroPeriodDelta,RetroPeriodICU_nation,retro_lb,retro_ub);


%--- Eliminate the effects of vaccination from delta ---%
delta_past_avg = delta_average; %Past 17 weeks average
delta_ss = delta_average*(0.1063/1.53);
VD_elderly = D1*V1_elderly + (D2-D1)*V2_elderly;
VD_medical = D1*V1_medical + (D2-D1)*V2_medical;
VD_ordinary = (D1*V1_medical + (D2-D1)*V2_medical) + (D1*V1_others + (D2-D1)*V2_others);
share = ((1 - sum(VD_elderly(1:end-2))/elderly_total) * (delta_average - delta_ss) ...
    + (1-sum(VD_ordinary(1:end-2))/(ordinary_total))*delta_ss)/delta_average;
delta_average = delta_average/share;

%--- Eliminate the effects of vaccination from delta ---%
ICU_ss = delta_ICU_average*(0.3916/1.62); %Share of the death rate among youth to the death rate of all populaiton
share = ((1 - sum(VD_elderly(1:end-2))/elderly_total) * (delta_ICU_average - ICU_ss) ...
    + (1-sum(VD_ordinary(1:end-2))/(ordinary_total))*ICU_ss)/delta_ICU_average;
ICU_average = delta_ICU_average/share; %Eliminate the effects of vaccination in the past average value


AverageAlpha2020 = 100*mean([alpha(1);alpha(2);alpha(1:End2020)]);

%--- Plot parameters in the data period ---%
ParamList = ["alpha","ERN","beta","delta"];
ParamName = ["$\alpha$ (decline in Y)","Effective reproduction number",...
    "$\beta$ (raw infection rate)","$\delta$ (death rate)"];

%%%%%%%%%%%%%%%%% Projection starts here %%%%%%%%%%%%%%%%%
beta_sample = beta(end-RetroPeriod+1:end);
beta_average = mean(beta_sample);
betaT = mean(beta_sample)*ones(SimPeriod,1);

if Vsimple == 0

    V1_prev = E1*(V1_elderly+V1_medical+V1_others)+(E2-E1)*(V2_elderly+V2_medical+V2_others);
    V = V1_prev(end-1);
    
    V2 = sum(D1*V1_elderly(1:end-1)+(D2-D1)*V2_elderly(1:end-1));
    V_ord = sum(D1*V1_medical(1:end-1)+(D2-D1)*V2_medical(1:end-1) + D1*V1_others(1:end-1)+(D2-D1)*V2_others(1:end-1));
    delta_ss = delta_average*(0.1063/1.53); % Sheet1 (age_heterogeneity/severe_risk_2021JUN06.xlsx)
    deltaT = (1-(cumsum(V2)/elderly_total))*(delta_average-delta_ss)+(1-(cumsum(V_ord)/(POP0-elderly_total)))*delta_ss;
    delta_ICU_ss = ICU_average*(0.3916/1.62); % Sheet2 (age_heterogeneity/severe_risk_2021JUN06.xlsx)
    delta_ICU = (1-(cumsum(V2)/elderly_total))*(ICU_average-delta_ICU_ss)+(1-(cumsum(V_ord)/(POP0-elderly_total)))*delta_ICU_ss;
    
else
    [V,deltaT,VT] = vaccine_distribution_simple(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2);
end

%--- Construct time series of parameters ---%
InitialValues = [S(end),I(end),R(end),D(end),ICU(end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forecast error analysis (using next week's data) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if iPC==1
%   NextD = load('\Users\shcor\Dropbox\fujii_nakata\Website\Codes\Japan20210316.mat');
%   FE = load('\Users\shcor\Dropbox\fujii_nakata\Website\Codes\Forecast.mat');
% else
%   NextD = load('/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/Japan20210316.mat');
%   FE = load('/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/Forecast.mat');
% end
%
% AlphaNext = NextD.alpha(end);
%[CumDNext,AverageAlphaNext,SimDataNext,SimNNext,SimICUNext] = ...
%Covid_projection2(InitialValues,AlphaNext,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,ICU_inflow_avg,gamma_ICU,delta_ICU); % [CumDNext,AverageAlphaNext,SimDataNext,SimNNext]=Covid_projection(InitialValues,AlphaNext,betaT,gammaT,deltaT,V,h,k,POP0,hconstant);
% dDForecast = [FE.dDForecast;(CumDNext-D(end))];
% dNForecast = [FE.dNForecast;SimNNext];
% ICUForecast = [FE.ICUForecast;SimICUNext(2)];
% if dNForecast(end) ~= dNForecast(end-1) && length(dNForecast) == NextD.Tdata
%   if iPC==1
%       save('\Users\shcor\Dropbox\fujii_nakata\Website\Codes\Forecast.mat','dDForecast','dNForecast','ICUForecast','-append');
%   else
%       save('/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/Forecast.mat','dDForecast','dNForecast','ICUForecast','-append');
%   end
% end

%---------------%
%insert %% here
%---------------%

if mat_save==1
%     save([pref char(datetime('today','Format','yyyMMdd')) '.mat'],'alpha','Tdata','AverageAlpha','CumD','LagResults','Month')
    save([pref char(datetime('today','Format','yyyMMdd')) '.mat'],'alpha','Tdata','Month')
end
