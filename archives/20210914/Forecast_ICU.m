
% This m-file executes simulation and generates figures for the main
% analysis of "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
iPC=1;
if iPC == 1
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
%     home = '/Users/shotaro/Dropbox/fujii_nakata/Website/Codes/';
    home = '/Users/koheimachi/Dropbox/fujii_nakata/Website/Codes/';
end
cd(home);

%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
mat_save = 0;    % 0 = figures won't be saved, 1 = they will be saved　in the "Figure" folder
fs = 16;            %common font size for many figures
xlim_tradeoff = [1,2.5];
if iPC == 1
    fn = 'Yu Gothic';
else
    fn = 'YuGothic';
end
%======================================================================%

%====================== Model parameter values ======================%
pref = 'Japan';        % prefecture to be analyzed
POP0 = 125710000;      % initial population
gamma = 7/12;          % recovery rate from Covid
k = 2;                 % exponent of (1-h*alpha)
hconstant = 1;         % 0 = without intercept, 1 = with intercept for h regression
TargetAlpha = 0.5:0.1:3;      % values of alpha we simulate results
AlphaVals = [1.2,1.65,2.5];   % benchmark alpha we plot time-series figures
SimPeriod = 52;        % simulation period in weeks
medical_start_date = datetime(2021,3,18);
elderly_start_date = datetime(2021,5,13);
RetroPeriod = 17;      % retroactive periods used to estimate gamma and delta
wl = [1,2];            % Results as of these weeks ago
ICU_nation = 1; % 1 for NHK, 0 for prefecture-based criteria
%====================================================================%


%--- Import data ---%
% Covid data are recorded at weekly frequencies (Mon-Sun)
% The first week start on January 20 (Mon), 2020
covid = importdata([home 'Covid_weekly_newV.csv']);  % Import weekly Covid data
Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
date = Data(:,1) + 21916;
date1 =datetime(date,'ConvertFrom','excel','Format','MMM-yy');
date = datetime(date,'ConvertFrom','excel');
N = Data(:,2);
dD = Data(:,3);
M = Data(:,4);
GDP = Data(:,5);
Tdata= size(Data,1);    % Data period in weeks
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
if isempty(find(SimDateEN == medical_start_date)) == 0
    medical_start = find(SimDateEN == medical_start_date);
else
    medical_start = 1;
end
elderly_start = find(SimDateEN == elderly_start_date);
VacStart = find(SimDateEN == datetime(2021,4,1));
End2020 = find(dateEN == datetime(2021,1,7));
Month = string(date1);

%--- Construct weekly vaccine data ---%
if iPC == 1
    vaccine_medical = importdata([home 'vaccine_daily_medical.xls']);
    vaccine_elderly = importdata([home 'vaccine_daily_elderly.xls']);
    [V1_medical, V2_medical] = vaccine_daily_to_weekly(vaccine_medical, ps, dateEN,iPC);
    [V1_elderly, V2_elderly] = vaccine_daily_to_weekly(vaccine_elderly, ps, dateEN,iPC);
else
    vaccine_medical = readmatrix([home 'vaccine_daily_medical.xls']);
    vaccine_elderly = readmatrix([home 'vaccine_daily_elderly.xls']);
    [V1_medical, V2_medical] = vaccine_daily_to_weekly_table(vaccine_medical, ps, dateEN,iPC);
    [V1_elderly, V2_elderly] = vaccine_daily_to_weekly_table(vaccine_elderly, ps, dateEN,iPC);
    
end


% vaccine pace
Vsimple = 0; % 0 for vaccine_distribution; 1 for vaccine_distribution_simple
PF = 1; % 0 for AZ, 1 for PF
Vgradual = 1; % 0 for flat Vpath, 1 for gradually increasing Vpath
if Vsimple == 0
    VP = [3:1:17]; %[0];
    VP_index = [3,8]; %[0,0];
else
    VP = [3]; %[0];
    VP_index = [3,8]; %[0,0];
end
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
% parameters for ICU
gamma_ICU = 7/28; % Recovery rate from ICU
ICU_adjustment = 0.75;
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

% %--- Import ICU data (Option 2) ---%
% if iPC==1
%     ICUdata = importdata([home 'ICU_total.csv']);  % Import weekly Covid data by prefecture
% else
%     ICUdata = importdata([home 'ICU_total.csv']);  % Import weekly Covid data by prefecture
% end
% ICUtotal = ICUdata.data(:,:);
% ICUtotal = [zeros(4,2);ICUtotal];
% dateICU  = ICUtotal(:,1) + 21916;
% dateICU_EN = datetime(dateICU,'ConvertFrom','excel');
% ICU = ICUtotal(:,2);
%--- Import ICU data (Option 2) ---%
ICU = zeros(Tdata+1,1);
BED = zeros(Tdata,1);
if ICU_nation == 1
    ICU(2:Tdata+1,1) = Data(:,22);
    BED(1:Tdata,1) = Data(:,23);
else
    ICU(2:Tdata+1,1) = Data(:,21);
end

%%%%%%%%%%%%%%%%% Forecasting accuracy analysis %%%%%%%%%%%%%%%%%
StartDateF  = 33;
HorizonValsF=[1, 4, 8];
dNForecast=zeros(Tdata,length(HorizonValsF));
dDForecast=zeros(Tdata,length(HorizonValsF));
ICUForecast=zeros(Tdata,length(HorizonValsF));
dNActual=zeros(Tdata,length(HorizonValsF));
dDActual=zeros(Tdata,length(HorizonValsF));
ICUActual=zeros(Tdata,length(HorizonValsF));

%--- Regress mobility on alpha to estimate the elasticity h ---%

for iH = 1:length(HorizonValsF)
    HorizonF  = HorizonValsF(iH);
    EndtDateF = Tdata+1-HorizonF;
    
    for iF = StartDateF:EndtDateF
        GDP_F = GDP(1:iF-6);  % why 6????
        referenceGDP_F = referenceGDP(1:iF-1);
        
        %--- Impute alpha (regress alpha on M)---%
        Malt_F=M(1:iF-1);
        alpha_F = (1 - GDP(1:iF-6)./referenceGDP(1:iF-6));   % output loss in percentage
        X_F = Malt_F(iF-6-17:iF-6);
        Y_F = alpha_F(iF-6-17:iF-6);
        XC_F = [ones(length(X_F),1), X_F];
        s_F = (XC_F'*XC_F)\XC_F'*Y_F;         % OLS estimate of h with constant
        reg_F = XC_F*s_F;
        r_F = Y_F - reg_F;
        eps_p_F = zeros(iF-1-(iF-6),1);
        eps_p_F(1) = r_F(end);
        for i = 1:iF-1-(iF-6)-1
            eps_p_F(i+1) = 1*eps_p_F(i);
        end
        alpha_pred_F = s_F(1)+s_F(2)*Malt_F((iF-6)+1:iF-1)+eps_p_F;
        
        alpha_F = [alpha_F;alpha_pred_F];
        
        %--- Regress mobility on alpha to estimate the elasticity h ---%
        Y_F = M(4:iF-6);
        X_F = alpha_F(4:iF-6);
        if hconstant == 0
            Y_F = Y_F - 1;
            h_all_F = (X_F'*X_F)\X_F'*Y_F;              % OLS estimate of h
            reg_F = X_F*h_all_F;
            r_F = Y_F - reg_F;                % r is the residuals, which is the observed minus fitted values
            SSE_F = sum(r_F.^2);            % SSE is the sum of squared errors
            MSE_F=SSE_F/(length(Y_F)-1);      % mean squared error
            h_all_se_F=sqrt(MSE_F/sum(X_F.^2));   % standard error of h
        elseif hconstant == 1
            XC_F = [ones(length(X_F),1), X_F];
            h_all_F = (XC_F'*XC_F)\XC_F'*Y_F;         % OLS estimate of h with constant
            reg_F = XC_F*h_all_F;
            r_F = Y_F - reg_F;
            SSE_F = sum(r_F.^2);
        end
        h_F = h_all_F;
        
        %--- Compute the history of S, I, R, D in the data period ---%
        S_F = zeros(iF,1);
        I_F = zeros(iF,1);
        R_F = zeros(iF,1);
        D_F = zeros(iF,1);
        S_F(1)=POP0;
        for i = 1:iF-1
            S_F(i+1)=S_F(i)-N(i)-E1*(V1_elderly(i)+V1_medical(i))-(E2-E1)*(V2_elderly(i)+V2_medical(i));
            I_F(i+1)=I_F(i)+N(i)-gamma*I_F(i)-dD(i);
            R_F(i+1)=R_F(i)+gamma*I_F(i)+E1*(V1_elderly(i)+V1_medical(i))+(E2-E1)*(V2_elderly(i)+V2_medical(i));
            D_F(i+1)=D_F(i)+dD(i);
        end
        ICU_F = ICU(1:iF);
        
        %--- Compute the history of time-varying parameters ---%
        
        delta_F = (D_F(2:iF)-D_F(1:iF-1))./I_F(1:iF-1);
        beta_tilde_F = (POP0.*N(1:iF-1))./(S_F(1:iF-1).*I_F(1:iF-1));   % overall infection rate
        if hconstant == 0
            beta_F = beta_tilde_F./(1+h_all_F*alpha_F).^k;                                      % raw infection rate
        elseif hconstant == 1
            %     beta = beta_tilde./(h(1)+h(2)*alpha).^k;
            beta_F = beta_tilde_F./(1+(h_all_F(2)/h_all_F(1))*alpha_F).^k;
        end
        ICU_inflow_F = (ICU_F(2:iF) - ICU_F(1:iF-1) + gamma_ICU.*ICU_F(1:iF-1) + dD(1:iF-1))./(delta_F(1:iF-1).*N(1:iF-1));
        
        IV_F     = [S_F(iF),I_F(iF),R_F(iF),D_F(iF),ICU_F(iF)];
        alphaT_F = alpha(iF:iF+HorizonF-1); %Use actual values because we are interested in conditional forecasts.
        beta_average_F = mean(beta_F(iF-RetroPeriod:iF-1));
        %         delta_average_F = mean(delta_F(iF-1-RetroPeriod:iF-1));
        delta_sample = delta_F(iF-1-RetroPeriod:iF-1);
        delta_average_F = sum(delta_sample.*(I_F(iF-1-RetroPeriod:iF-1)/sum(I_F(iF-1-RetroPeriod:iF-1))));
        
        ICU_inflow_avg_F = mean(ICU_inflow_F(iF-1-RetroPeriod:iF-1))*ICU_adjustment;
        
        betaT_F = beta_average_F*ones(HorizonF,1);
        deltaT_F = delta_average_F*ones(HorizonF,1);
        V_F = zeros(HorizonF,1);
        gammaT_F=gamma*ones(HorizonF,1);
        
        
        %         [CumD_F,AverageAlpha_F,SimData_F,SimN_F]=Covid_projection(IV_F,alphaT_F,betaT_F,gammaT_F,deltaT_F,V_F,h_F,k,POP0,hconstant);
        [CumD_F,AverageAlpha_F,SimData_F,SimN_F,SimICU_F]=Covid_projection2(IV_F,alphaT_F,betaT_F,gammaT_F,deltaT_F,V_F,h_F,k,POP0,hconstant,ICU_inflow_avg_F,gamma_ICU);
        
        dNForecast(iF,iH) = sum(SimN_F(1:HorizonF));
        dDForecast(iF,iH) = CumD_F-IV_F(4);
        ICUForecast(iF,iH) = SimICU_F(HorizonF+1); %sum(SimICU_F(2:HorizonF+1));
        dNActual(iF,iH)   = sum(N(iF:iF+HorizonF-1));
        dDActual(iF,iH)   = D(iF+HorizonF-1)-D(iF-1);
        ICUActual(iF,iH)   = ICU(iF+HorizonF); %sum(ICU(iF:iF+HorizonF));
        
    end
end

xtick1_F=[StartDateF+4:8:Tdata];
xtick4_F=[StartDateF+4:8:Tdata];
xtick8_F=[StartDateF+4:8:Tdata];
fs_F=10;
fs_legend_F=10;
fs_title=12;
figure(10)
subplot(2,2,1)
plot(StartDateF:(Tdata),dNForecast(StartDateF:Tdata,1),'r',StartDateF:(Tdata),dNActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
title({'Conditional Forecast vs. Actual';'(one-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF (Tdata)])
xticks(xtick1_F)
xticklabels(Month(xtick1_F))
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ytickformat('%,6.0f')
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,2)
plot(StartDateF+3:(Tdata),dNForecast(StartDateF:Tdata-3,2),'r',StartDateF+3:(Tdata),dNActual(StartDateF:Tdata-3,2),'k','LineWidth',1.5)
title({'Conditional Forecast vs. Actual';'(four-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+3 (Tdata)])
xticks(xtick4_F)
xticklabels(Month(xtick4_F))
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,3)
plot(StartDateF:(Tdata),dNForecast(StartDateF:Tdata,1)-dNActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
title({'Conditional Forecast Errors',' (one-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF (Tdata)])
xticks(xtick1_F)
xticklabels(Month(xtick1_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
subplot(2,2,4)
plot(StartDateF+3:(Tdata),dNForecast(StartDateF:Tdata-3,2)-dNActual(StartDateF:Tdata-3,2),'k','LineWidth',1.5)
title({'Conditional Forecast Errors',' (four-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+3 (Tdata)])
xticks(xtick4_F)
xticklabels(Month(xtick4_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

figure(11)
subplot(2,2,1)
plot(StartDateF:(Tdata),dDForecast(StartDateF:Tdata,1),'r',StartDateF:(Tdata),dDActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
title({'Conditional Forecast vs. Actual';'(one-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
xlim([StartDateF (Tdata)])
xticks(xtick1_F)
xticklabels(Month(xtick1_F))
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,2)
plot(StartDateF+3:(Tdata),dDForecast(StartDateF:Tdata-3,2),'r',StartDateF+3:(Tdata),dDActual(StartDateF:Tdata-3,2),'k','LineWidth',1.5)
title({'Conditional Forecast vs. Actual';'(four-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+3 (Tdata)])
xticks(xtick4_F)
xticklabels(Month(xtick4_F))
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,3)
plot(StartDateF:(Tdata),dDForecast(StartDateF:Tdata,1)-dDActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
title({'Conditional Forecast Errors',' (one-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF (Tdata)])
xticks(xtick1_F)
xticklabels(Month(xtick1_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
subplot(2,2,4)
plot(StartDateF+3:(Tdata),dDForecast(StartDateF:Tdata-3,2)-dDActual(StartDateF:Tdata-3,2),'k','LineWidth',1.5)
title({'Conditional Forecast Errors',' (four-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+3 (Tdata)])
xticks(xtick4_F)
xticklabels(Month(xtick4_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

figure(12)
subplot(2,2,1)
plot(StartDateF:(Tdata),ICUForecast(StartDateF:Tdata,1),'r',StartDateF:(Tdata),ICUActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
title({'Conditional Forecast vs. Actual';'(one-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF (Tdata)])
xticks(xtick1_F)
xticklabels(Month(xtick1_F))
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ytickformat('%,6.0f')
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,2)
plot(StartDateF+3:(Tdata),ICUForecast(StartDateF:Tdata-3,2),'r',StartDateF+3:(Tdata),ICUActual(StartDateF:Tdata-3,2),'k','LineWidth',1.5)
title({'Conditional Forecast vs. Actual';'(four-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+3 (Tdata)])
xticks(xtick4_F)
xticklabels(Month(xtick4_F))
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,3)
plot(StartDateF:(Tdata),ICUForecast(StartDateF:Tdata,1)-ICUActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
title({'Conditional Forecast Errors',' (one-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF (Tdata)])
xticks(xtick1_F)
xticklabels(Month(xtick1_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
subplot(2,2,4)
plot(StartDateF+3:(Tdata),ICUForecast(StartDateF:Tdata-3,2)-ICUActual(StartDateF:Tdata-3,2),'k','LineWidth',1.5)
title({'Conditional Forecast Errors',' (four-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+3 (Tdata)])
xticks(xtick4_F)
xticklabels(Month(xtick4_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

figure(13)
subplot(2,2,1)
plot(StartDateF+7:(Tdata),ICUForecast(StartDateF:Tdata-7,3),'r',StartDateF+7:(Tdata),ICUActual(StartDateF:Tdata-7,3),'k','LineWidth',1.5)
title({'Conditional Forecast vs. Actual';'(eight-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+7 (Tdata)])
xticks(xtick8_F)
xticklabels(Month(xtick8_F))
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ytickformat('%,6.0f')
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,3)
plot(StartDateF+7:(Tdata),ICUForecast(StartDateF:Tdata-7,3)-ICUActual(StartDateF:Tdata-7,3),'k','LineWidth',1.5)
title({'Conditional Forecast Errors',' (eight-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+7 (Tdata)])
xticks(xtick8_F)
xticklabels(Month(xtick8_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

if figure_save == 1
    saveas(figure(12),[home 'Figures/Forecast_ICU.png']);
end


