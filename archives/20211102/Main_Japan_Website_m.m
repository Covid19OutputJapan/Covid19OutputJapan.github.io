% This m-file executes simulation and generates figures for the main
% analysis of "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
iPC= 0;
% home = '/Users/shotaro/Dropbox/fujii_nakata/Website/Covid19OutputJapan.github.io/archives/20211012/';
home = '/Users/shotaro/Dropbox/fujii_nakata/Website/Covid19OutputJapan.github.io/archives/20211012/';
home_now = '/Users/shotaro/Dropbox/fujii_nakata/Website/Codes/';
figures_dir = [home_now '/Figures/'];
if iPC == 1
    %home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
%     home = '\Users\kenic\Dropbox\fujii_nakata\Website\Codes\';
else
    % home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/';
%     home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/';
    %home = '/Users/shotaro/Dropbox/fujii_nakata/Website/Codes/';
    %home = '/Users/koheimachi/Dropbox/fujii_nakata/Website/Codes/';
end
cd(home_now);

%====================== Program parameter values ======================%
mat_save = 1; 
figure_save = 1;
data_save = 1;

week_diff = 3;

data_switch = 0;
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
SimPeriod = week_diff;
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
covid = importdata([home_now 'Covid_weekly_newV.csv']);  % Import weekly Covid data
create_japan_data
%--- Construct weekly vaccine data ---%
vaccine_medical = readmatrix([home_now 'vaccine_daily_medical.xls']);
vaccine_elderly = readmatrix([home_now 'vaccine_daily_elderly.xls']);
[V1_medical, V2_medical] = vaccine_daily_to_weekly_table(vaccine_medical, ps, date(1:Tdata+3),iPC);
V1_medical(end) = V1_medical(end);
V2_medical(end) = V2_medical(end);
[V1_elderly, V2_elderly] = vaccine_daily_to_weekly_table(vaccine_elderly, ps, date(1:Tdata+3),iPC);
V1_elderly(end) = V1_elderly(end);
V2_elderly(end) = V2_elderly(end);
vaccine_others = readmatrix([home_now 'vaccine_daily_others.xls']);
[V1_others, V2_others] = vaccine_daily_to_weekly_table(vaccine_others, ps, date(1:Tdata+3),iPC);
V1_others(end) = V1_others(end);
V2_others(end) = V2_others(end);
% vaccine pace
Vsimple = 0; % 0 for vaccine_distribution; 1 for vaccine_distribution_simple
PF = 1; % 0 for AZ, 1 for PF

%--- Update forecast error data ---%
FE = load([home 'Forecast.mat']);
dDActual = [FE.dDActual;dD(end-SimPeriod+1:end)];
dNActual = [FE.dNActual;N(end-SimPeriod+1:end)];
ICUActual = [FE.ICUActual;floor(ICU(end-SimPeriod+1:end))];

if dNActual(end) ~= dNActual(end-1) && length(dNActual) == Tdata
    %     save([home 'Forecast.mat'],'dDActual','dNActual','-append');
    save([home_now 'Forecast.mat'],'dDActual','dNActual','ICUActual');
end

% if mat_save==1
%     save([pref char(datetime('today','Format','yyyMMdd')) '.mat'],'alpha','Tdata','Month')
% end
newalpha = alpha;
newTdata = Tdata;
newMonth = Month;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Use past data to simulate for SimPeriod (week_diff)%%%%%%%%%%%%%%%%%%%%%%%
parameter
SimPeriod = week_diff;
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
covid = importdata([home 'Covid_weekly_newV.csv']);  % Import weekly Covid data
create_japan_data

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
    pastV = zeros(Tdata+3, 1);
    pastV(3:end) = E1 * (V1_elderly(1:end - 2) + V1_medical(1:end - 2) + V1_others(1:end - 2)) + (E2 - E1) * (V2_elderly(1:end - 2) + V2_medical(1:end - 2) + V2_others(1:end - 2));
    S = S_data;
    S(2:end) = S(2:end) - cumsum(pastV(1:Tdata));
    R = R_data;
    R(2:end) = R(2:end) + cumsum(pastV(1:Tdata));
    D = D_data;
end

[delta,beta_tilde,ERN,beta,ICU_inflow,...
    gammaT,delta_average,delta_ICU_average,ICU_inflow_avg,delta_sample,beta_avg]...
    = Time_Series_Average_Japan(S,I,D,ICU,dD,N,Tdata,SimPeriod,...
    RetroPeriod,POP0,gamma,hconstant,h_all,alpha,k,...
    gamma_ICU,ICU_adjustment,RetroPeriodDelta,RetroPeriodICU_nation,retro_lb,retro_ub);


%--- Eliminate the effects of vaccination from delta ---%
delta_past_avg = delta_average; %Past 17 weeks average
% delta_ss = delta_average*(0.1063/1.53);
VD_elderly = D1*V1_elderly + (D2-D1)*V2_elderly;
VD_medical = D1*V1_medical + (D2-D1)*V2_medical;
VD_ordinary = (D1*V1_medical + (D2-D1)*V2_medical) + (D1*V1_others + (D2-D1)*V2_others);

%%%%%%%%%%%%%%%%% Projection starts here %%%%%%%%%%%%%%%%%
beta_sample = beta(end-RetroPeriod+1:end);
beta_average = mean(beta_sample);
betaT = mean(beta_sample)*ones(SimPeriod,1);

if Vsimple == 0

    V1_prev = E1*(V1_elderly+V1_medical+V1_others)+(E2-E1)*(V2_elderly+V2_medical+V2_others);
    
    V = V1_prev(Tdata-1:Tdata-1+SimPeriod-1);

    V2 = sum(D1*V1_elderly(1:end-1)+(D2-D1)*V2_elderly(1:end-1));
    V_ord = sum(D1*V1_medical(1:end-1)+(D2-D1)*V2_medical(1:end-1) + D1*V1_others(1:end-1)+(D2-D1)*V2_others(1:end-1));
    
    S_elderly_prev = elderly_total - sum(E1*V1_elderly(1:Tdata-2) + (E2 - E1)*V2_elderly(1:Tdata-2));
    V1_elderly_Sim = V1_elderly(Tdata-1:Tdata+1);
    V2_elderly_Sim = V2_elderly(Tdata-1:Tdata+1);
    S_elderly_path = S_elderly_prev - cumsum(E1 * V1_elderly_Sim+ (E2 - E1) * V2_elderly_Sim);

    S_young_prev =  ordinary_total + medical_total...
        - sum(E1*(V1_others(1:Tdata-2) + V1_medical(1:Tdata-2))) ...
        - sum((E2 - E1)*(V2_others(1:Tdata-2) + V2_medical(1:Tdata-2)));
    
    V1_others_Sim = V1_others(Tdata-1:Tdata-1+SimPeriod-1);
    V2_others_Sim = V2_others(Tdata-1:Tdata-1+SimPeriod-1);
    V1_medical_Sim = V1_medical(Tdata-1:Tdata-1+SimPeriod-1);
    V2_medical_Sim = V1_medical(Tdata-1:Tdata-1+SimPeriod-1);
    V1_young_Sim = V1_others_Sim + V1_medical_Sim;
    V2_young_Sim = V2_others_Sim + V2_medical_Sim;
    S_young_path = S_young_prev - cumsum(E1 * V1_young_Sim+ (E2 - E1) * V2_young_Sim);
    elderly_share_Sim = S_elderly_path./(S_young_path+S_elderly_path);
    elderly_share_Tdata = S_elderly_prev/(S_young_prev+S_elderly_prev);

    deltaT = elderly_share_Sim./elderly_share_Tdata * delta_average .* lambda_delta(1:SimPeriod);
    delta_ICU = elderly_share_Sim./elderly_share_Tdata * delta_ICU_average .* lambda_delta_ICU_nation(1:SimPeriod);

else
    [V,deltaT,VT] = vaccine_distribution_simple(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2);
end

%--- Construct time series of parameters ---%
InitialValues = [S(Tdata+1),I(Tdata+1),R(Tdata+1),D(Tdata+1),ICU(Tdata+1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forecast error analysis (using next week's data) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NextD = load([home_now 'Japan20211102.mat']);
FE = load([home 'Forecast.mat']);

% AlphaNext = NextD.alpha(end-SimPeriod+1:end);
AlphaNext = newalpha(end-SimPeriod+1:end);
[CumDNext,AverageAlphaNext,SimDataNext,SimNNext,SimICUNext] = ...
Covid_projection2(InitialValues,AlphaNext,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,ICU_inflow_avg,gamma_ICU,delta_ICU); % [CumDNext,AverageAlphaNext,SimDataNext,SimNNext]=Covid_projection(InitialValues,AlphaNext,betaT,gammaT,deltaT,V,h,k,POP0,hconstant);
dDForecast = [FE.dDForecast;(CumDNext(1+1:1+SimPeriod)-CumDNext(1:SimPeriod))];
dNForecast = [FE.dNForecast;SimNNext];
ICUForecast = [FE.ICUForecast;SimICUNext(2:SimPeriod+1)];
save([home_now 'Forecast.mat'],'dDForecast','dNForecast','ICUForecast');

%%%%%%%%%%%%%%% make forecast error plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CurrentD = load([home_now 'Japan20211102.mat']);

Tdata =newTdata;
Month = newMonth;

markersize = 6; %marker size

StartDateF  = 36;
format short
xtick1_F=StartDateF+4:8:Tdata;
xtick4_F=StartDateF+4:8:Tdata;
fs_F=10;
fs_legend_F=10;
fs_title=12;

for l = 1:2
    if l == 1
        f=figure(10)
    elseif l == 2
        f=figure(101)
    end
    f.WindowState = 'maximized';
    subplot(2,3,1)
    plot(StartDateF:(Tdata),dNForecast(StartDateF:Tdata,1),'r',StartDateF:(Tdata),dNActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
    hold on
    plot(Tdata,dNForecast(Tdata, 1),'or','MarkerFaceColor','r', 'MarkerSize', markersize);
    hold on
    plot((Tdata-SimPeriod),dNForecast(Tdata-SimPeriod, 1),'or','MarkerFaceColor','r', 'MarkerSize', markersize);
    if l == 1
        title({'Conditional Forecast vs. Actual';'(New cases)'},'FontSize',fs_title,'FontWeight','normal')
        legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
    elseif l == 2
        title({'予測と実現値';'(新規感染者数)'},'FontSize',fs_title,'FontWeight','normal','FontName',fn)
        legend('予測','実現値','FontSize',fs_legend_F,'Location','Northwest','FontName',fn)
    end
    xlim([StartDateF (Tdata)])
    xticks(xtick1_F)
    xticklabels(Month(xtick1_F))
    ax = gca;
    ax.YAxis.Exponent = 0;
    ax.YAxis.FontSize = fs_F;
    ax.XAxis.FontSize = fs_F;
    ytickformat('%,6.0f')

    subplot(2,3,2)
    plot(StartDateF:(Tdata),dDForecast(StartDateF:Tdata,1),'r',StartDateF:(Tdata),dDActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
    hold on
    plot(Tdata,dDForecast(Tdata, 1),'or','MarkerFaceColor','r', 'MarkerSize', markersize);
    hold on
    plot((Tdata-SimPeriod),dDForecast(Tdata-SimPeriod, 1),'or','MarkerFaceColor','r', 'MarkerSize', markersize);
    if l == 1
        title({'Conditional Forecast vs. Actual';'(New deaths)'},'FontSize',fs_title,'FontWeight','normal')
        legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
    elseif l == 2
        title({'予測と実現値';'(新規死亡者数)'},'FontSize',fs_title,'FontWeight','normal','FontName',fn)
        legend('予測','実現値','FontSize',fs_legend_F,'Location','Northwest','FontName',fn)
    end
    ax = gca;
    ax.YAxis.FontSize = fs_F;
    ax.XAxis.FontSize = fs_F;
    ax.YAxis.Exponent = 0;
    ytickformat('%,6.0f')
    xlim([StartDateF (Tdata)])
    xticks(xtick1_F)
    xticklabels(Month(xtick1_F))

    subplot(2,3,3)
    plot(StartDateF:(Tdata),ICUForecast(StartDateF:Tdata,1),'r',StartDateF:(Tdata),ICUActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
    hold on
    plot(Tdata,ICUForecast(Tdata, 1),'or','MarkerFaceColor','r', 'MarkerSize', markersize);
    hold on
    plot((Tdata-SimPeriod),ICUForecast(Tdata-SimPeriod, 1),'or','MarkerFaceColor','r', 'MarkerSize', markersize);
    if l == 1
        title({'Conditional Forecast vs. Actual';'(ICU)'},'FontSize',fs_title,'FontWeight','normal')
        legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
    elseif l == 2
        title({'予測と実現値';'(重症者数)'},'FontSize',fs_title,'FontWeight','normal','FontName',fn)
        legend('予測','実現値','FontSize',fs_legend_F,'Location','Northwest','FontName',fn)
    end
    xlim([StartDateF (Tdata)])
    xticks(xtick1_F)
    xticklabels(Month(xtick1_F))
    ax = gca;
    ax.YAxis.Exponent = 0;
    ax.YAxis.FontSize = fs_F;
    ax.XAxis.FontSize = fs_F;
    ytickformat('%,6.0f')

    subplot(2,3,4)
    plot(StartDateF:(Tdata),dNForecast(StartDateF:Tdata,1)-dNActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
    hold on
    plot(Tdata,dNForecast(Tdata,1)-dNActual(Tdata,1),'ok','MarkerFaceColor','k', 'MarkerSize', markersize);
    hold on
    plot((Tdata-SimPeriod),dNForecast(Tdata-SimPeriod,1)-dNActual(Tdata-SimPeriod,1),'ok','MarkerFaceColor','k', 'MarkerSize', markersize);
    if l == 1
        title({'Conditional Forecast Errors',' (New cases)'},'FontSize',fs_title,'FontWeight','normal')
    elseif l == 2
        title({'予測誤差',' (新規感染者数)'},'FontSize',fs_title,'FontWeight','normal','FontName',fn)
    end
    xlim([StartDateF (Tdata)])
    xticks(xtick1_F)
    xticklabels(Month(xtick1_F))
    yline(0)
    ax = gca;
    ax.YAxis.FontSize = fs_F;
    ax.XAxis.FontSize = fs_F;
    ax.YAxis.Exponent = 0;
    ytickformat('%,6.0f')

    subplot(2,3,5)
    plot(StartDateF:(Tdata),dDForecast(StartDateF:Tdata,1)-dDActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
    hold on
    plot(Tdata,dDForecast(Tdata,1)-dDActual(Tdata,1),'ok','MarkerFaceColor','k', 'MarkerSize', markersize);
    hold on
    plot((Tdata-SimPeriod),dDForecast(Tdata-SimPeriod,1)-dDActual(Tdata-SimPeriod,1),'ok','MarkerFaceColor','k', 'MarkerSize', markersize);
    if l == 1
        title({'Conditional Forecast Errors',' (New deaths)'},'FontSize',fs_title,'FontWeight','normal')
    elseif l == 2
        title({'予測誤差',' (新規死亡者数)'},'FontSize',fs_title,'FontWeight','normal','FontName',fn)
    end
    xlim([StartDateF (Tdata)])
    xticks(xtick1_F)
    xticklabels(Month(xtick1_F))
    yline(0)
    ax = gca;
    ax.YAxis.FontSize = fs_F;
    ax.XAxis.FontSize = fs_F;
    ax.YAxis.Exponent = 0;
    ytickformat('%,6.0f')

    subplot(2,3,6)
    plot(StartDateF:(Tdata),ICUForecast(StartDateF:Tdata,1)-ICUActual(StartDateF:Tdata,1),'k','LineWidth',1.5)
    hold on
    plot(Tdata,ICUForecast(Tdata,1)-ICUActual(Tdata,1),'ok','MarkerFaceColor','k', 'MarkerSize', markersize);
    hold on
    plot((Tdata-SimPeriod),ICUForecast(Tdata-SimPeriod,1)-ICUActual(Tdata-SimPeriod,1),'ok','MarkerFaceColor','k', 'MarkerSize', markersize);
    if l == 1
        title({'Conditional Forecast Errors',' (ICU)'},'FontSize',fs_title,'FontWeight','normal')
    elseif l == 2
        title({'予測誤差',' (重症者数)'},'FontSize',fs_title,'FontWeight','normal','FontName',fn)
    end
    xlim([StartDateF (Tdata)])
    xticks(xtick1_F)
    xticklabels(Month(xtick1_F))
    yline(0)
    ax = gca;
    ax.YAxis.FontSize = fs_F;
    ax.XAxis.FontSize = fs_F;
    ax.YAxis.Exponent = 0;
    ytickformat('%,6.0f')
end

if figure_save == 1
    saveas(figure(10),[figures_dir 'ForecastErrors.png']);
    saveas(figure(101),[figures_dir 'ForecastErrors_jp.png']);
end

if data_save == 1
    FE_table = round([dNForecast(end),dNActual(end),dNForecast(end)-dNActual(end);
        dDForecast(end),dDActual(end),dDForecast(end)-dDActual(end);
        ICUForecast(end),ICUActual(end),ICUForecast(end)-ICUActual(end)]);
    FET = table(["Forecast","Actual","Error";FE_table]);
    writetable(FET,[figures_dir 'ForecastTable.xls'],'WriteVariableNames',false);
end
