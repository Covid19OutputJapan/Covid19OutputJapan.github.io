% This m-file executes simulation and generates figures for the main
% analysis of "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
home = '/home/takeki/Dropbox/covid19outputjapan/code/';
% home = pwd;
%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
% in the "Figure" folder
fs = 16;            % common font size for many figures
xlim_tradeoff = [2,5];
%======================================================================%

%====================== Model parameter values ======================%
POP0 = 125710000;      % initial population
gamma = 7/10;          % recovery rate from Covid
k = 2;                 % exponent of (1-h*alpha)
hconstant = 0;         % 0 = without intercept, 1 = with intercept for h regression
TargetAlpha = 0:0.1:8;      % values of alpha we simulate results
AlphaIndex = [3.5,4,5.5];   % benchmark alpha we plot time-series figures
SimPeriod = 52;        % simulation period in weeks
VacStart = 12;         % time until the start of vaccination process
VacPace = 0.8*1500000; % number of vaccinated persons per week (steady pace)
VacDuration = 6;       % time until vaccination pace reaches its steady pace
RetroPeriod = 13;      % retroactive periods used to estimate gamma and delta
wl = [1,2];            % Results as of these weeks ago
%====================================================================%


%--- Import data ---%
% Covid data are recorded at weekly frequencies (Mon-Sun)
% The first week start on January 20 (Mon), 2020
covid = importdata([home 'Data/Covid_weekly_Japan.csv']);  % Import weekly Covid data
Data = covid.data;      % Columns: 1 = month, 1 = new positive, 2 = death, 3= mobility, 4 = GDP
N = Data(:,1);
dD = Data(:,2);
M = Data(:,3);
GDP = Data(:,4);
Tdata= size(Data,1);    % Data period in weeks
Month = covid.textdata(2:end,1);
xtick1 = [1:13:Tdata, Tdata];

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
% gapGDP = zeros(length(potentialGDP),1);
% gapGDP(1) = 0.0166;
% for i = 2:length(gapGDP)
%     gapGDP(i) = gapGDP(i-1)*(0.975^(12/52));
% end
referenceGDP = potentialGDP.*(1+0.0166);
referenceGDP(1:2) = [];


%--- Plot mobility data ---%
figure(2)
yyaxis left
plot(M,'k','LineWidth',1.5)
ylabel('Mobility')
yyaxis right
plot(GDP,'b-.','LineWidth',1.5)
ylabel('GDP')
xlim([1 Tdata])
xticks(xtick1)
xticklabels(Month(xtick1))
legend('Mobility (left axis)','GDP (right axis)','FontSize',16,'Location','north');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
ax.YAxis(1).FontSize = fs;
ax.YAxis(2).FontSize = fs;
ax.XAxis.FontSize = fs;


%--- Regress mobility on alpha to estimate the elasticity h ---%
X = 100*(1 - GDP./referenceGDP(1:Tdata));   % output loss in percentage
RM = [M,X];
RM(any(isnan(RM), 2), :) = [];       % delete missing values
Y = RM(:,1);
X = RM(:,2);
if hconstant == 0
    h = -(X'*X)\X'*Y;              % OLS estimate of h
    reg = -X*h;
    r = Y - reg;                % r is the residuals, which is the observed minus fitted values
    SSE = sum(r.^2);            % SSE is the sum of squared errors
    MSE=SSE/(length(Y)-1);      % mean squared error
    h_se=sqrt(MSE/sum(X.^2));   % standard error of h
elseif hconstant == 1
    XC = [ones(length(X),1), X];
    h_all = -(XC'*XC)\XC'*Y;         % OLS estimate of h with constant
    reg = -XC*h_all;
    Y2 = Y(end-12:end);
    XC2 = XC(end-12:end,:);
    h = -(XC2'*XC2)\XC2'*Y2;
end

figure(100)
scatter(X,Y,60,'filled')
hold on
plot(X,reg,'LineWidth',1.5)
xlabel('Output loss (%)');
ylabel('Mobility');
aa = gca;
aa.YAxis.FontSize = fs;
aa.XAxis.FontSize = fs;
hold off


%--- Compute the history of S, I, R, D in the data period ---%
S = zeros(Tdata+1,1);
I = zeros(Tdata+1,1);
R = zeros(Tdata+1,1);
D = zeros(Tdata+1,1);
S(1)=POP0;
Malt=M;
Malt(50)=0.5*(Malt(49)+Malt(51));   % smooth out the new year irregularity
for i = 1:Tdata
    S(i+1)=S(i)-N(i);
    I(i+1)=I(i)+N(i)-gamma*I(i)-dD(i);
    R(i+1)=R(i)+gamma*I(i);
    D(i+1)=D(i)+dD(i); 
    if isnan(GDP(i)) == 0
        lastGDP = i;
    else
        if hconstant == 0
         GDP(i) = referenceGDP(i)*(1+0.5*mean(Malt(i-1:i))/(100*h));
        %GDP(i) = 0.5*GDP(lastGDP)+0.5*referenceGDP(i)*(1+mean(M(i-3:i))/(100*h));
        elseif hconstant == 1
            GDP(i) = referenceGDP(i)*(1+0.5*(mean(M(i-1:i))+h(1))/(100*h_all(2)));
        end
    end
end

%--- Compute the history of time-varying parameters ---%
alpha = 1 - GDP./referenceGDP(1:Tdata);                                     % output loss
delta = (D(2:Tdata+1)-D(1:Tdata))./I(1:Tdata);                              % death rate
beta_tilde = -POP0*((S(2:Tdata+1)-S(1:Tdata))./(S(1:Tdata).*I(1:Tdata)));   % overall infection rate
ERN = (S(1:end-1)/POP0).*beta_tilde./(gamma+delta);                                        % effective reproduction number
if hconstant == 0
    beta = beta_tilde./(1-h*alpha).^k;                                      % raw infection rate
elseif hconstant == 1
    beta = beta_tilde./(1-0.01*h(1)-h(2)*alpha).^k;          
end

%--- Plot varabiles in the data period ---%
figure(1)
plot(N,'k','LineWidth',2)
xlim([1 Tdata])
xticks(xtick1)
yticks([0 10000 20000 30000 40000])
xticklabels(Month(xtick1))
xtickangle(45)
ax = gca;
ax.YAxis.FontSize = 16;
ax.XAxis.FontSize = 16;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')


VarList = ["S","I","R","D","N","GDP"];
VarName = ["S","I","R","D","N","Y"];

for i = 1:length(VarList)
    figure(4)
    subplot(2,3,i)
    plot(eval(VarList(i)),'k','LineWidth',2)
    xlim([1 Tdata])
    title(VarName(i),'FontWeight','normal')
    xticks(xtick1)
    xticklabels(Month(xtick1))
    xtickangle(45)
    if i == 5
        ax = gca;
        ax.YAxis.Exponent = 0;
        ytickformat('%,6.0f')
    end
    if i == length(VarList)
        hold on 
        plot(referenceGDP,'k--','LineWidth',2)
    end
end


 %--- Plot parameters in the data period ---%
ParamList = ["alpha","ERN","beta","delta"];
ParamName = ["$\alpha$ (decline in Y)","Effective reproduction number",...
    "$\beta$ (raw infection rate)","$\delta$ (death rate)"];
for i = 1:4
    figure(3)
    subplot(2,2,i)
    plot(eval(ParamList(i)),'k','LineWidth',2)
    xlim([1 Tdata])
    title(ParamName(i),'Interpreter','latex','FontSize',14,'FontWeight','normal')
    xticks(xtick1)
    xticklabels(Month(xtick1))
    xtickangle(45)
    if i == 2
        ylim([0 3]);
        yline(1);
    end
end

%%%%%%%%%%%%%%%%% Projection starts here %%%%%%%%%%%%%%%%%

%--- Construct time series of parameters ---%
InitialValues = [S(end),I(end),R(end),D(end)];
beta_sample = beta(end-RetroPeriod+1:end);
delta_sample = delta(end-RetroPeriod+1:end);
beta_average = mean(beta_sample);
delta_average = mean(delta_sample);
betaT = beta_average*ones(SimPeriod,1);
deltaT = delta_average*ones(SimPeriod,1);
gammaT = gamma*ones(SimPeriod,1);
% V = zeros(SimPeriod,1);
% V(VacStart:end) = VacPace;
V = zeros(SimPeriod,1);
V(VacStart-1:VacStart+VacDuration-1) = 0:VacPace/VacDuration:VacPace;
V(VacStart+VacDuration:end) = VacPace;
CumD = zeros(1,length(TargetAlpha));
AverageAlpha = zeros(1,length(TargetAlpha));
LagResults = zeros(2,length(TargetAlpha),length(wl));
ParamList2 = ["alpha2","ERN2","beta2","delta2"];


%--- Sensitivity parameters ---%
BetaVals = beta_average*[1.05, 1, 0.95];
HVals = h*[0.8, 1, 1.2];
KVals = [1.5, 2, 2.5];
VacStartVals = [16, 12, 20];
VacPaceVals = [0.5*VacPace, VacPace, 2*VacPace];
FrontLoadingVals = [0.5, 1, 1.5];
SA = zeros(3,length(TargetAlpha),4);    % matrix to contain the results of sensitivity analysis
AlphaIndexVariables = zeros(SimPeriod,6,length(AlphaIndex));    % matrix to contain the projected variables for certain alpha values

alphaTT = 0:TargetAlpha(16)*2/(SimPeriod-1):(TargetAlpha(16)*2);
%--- Simulation for different alpha values ---%
for i = 1:length(TargetAlpha)
    %%% Baseline projection %%%
    alphaT = flip(0:0.01*TargetAlpha(i)*2/(SimPeriod-1):(0.01*TargetAlpha(i)*2))';
    [CumD(i),AverageAlpha(i),SimData,SimN]=Covid_projection(InitialValues,alphaT,betaT,gammaT,deltaT,V,h,k,POP0,hconstant);
    %%% Record the paths of projected variables for certian alpha values 
    if sum(TargetAlpha(i)==AlphaIndex)==1
        j = find(TargetAlpha(i)==AlphaIndex);
        SimGDP = (1-alphaT)*GDP(1);
        AlphaIndexVariables(:,:,j) = [SimData(1:end-1,:),SimN,SimGDP];
    end
    %%% Results from previous weeks %%%
    for j = 1:length(wl)
        wlag = wl(j);
        InitialValuesL = [S(end-wlag),I(end-wlag),R(end-wlag),D(end-wlag)];
        betaL = mean(beta(end-RetroPeriod+1-wlag:end-wlag))*ones(SimPeriod,1);
        deltaL = mean(delta(end-RetroPeriod+1-wlag:end-wlag))*ones(SimPeriod,1);
        [LagResults(1,i,j),LagResults(2,i,j)]=Covid_projection(InitialValuesL,alphaT,betaL,gammaT,deltaL,V,h,k,POP0,hconstant);
    end
        
    %%% Beta sensitivity %%%
    for j = 1:length(BetaVals)
        betaS = BetaVals(j)*ones(SimPeriod,1);
        [SA(j,i,1)]=Covid_projection(InitialValues,alphaT,betaS,gammaT,deltaT,V,h,k,POP0,hconstant);
    end
    %%% h sensitivity %%%
    for j = 1:length(HVals)
        HS = HVals(j);
        [SA(j,i,2)]=Covid_projection(InitialValues,alphaT,betaT,gammaT,deltaT,V,HS,k,POP0,hconstant);
    end
    %%% k sensitivity %%%
    for j = 1:length(KVals)
        kS = KVals(j);
        [SA(j,i,3)]=Covid_projection(InitialValues,alphaT,betaT,gammaT,deltaT,V,h,kS,POP0,hconstant);
    end
    %%% VacStart sensitivity %%%
    for j = 1:length(VacStartVals)
        VacStartS = VacStartVals(j);
        VS = zeros(SimPeriod,1);
        VS(VacStartS:end) = VacPace;
        [SA(j,i,4)]=Covid_projection(InitialValues,alphaT,betaT,gammaT,deltaT,VS,h,k,POP0,hconstant);
    end
end

%%% Time series figures (variables and parameters) %%%

Past = [N,GDP,ERN,zeros(Tdata,1)];
VarList2 = ["N","GDP","ERN","V"];
VarName2 = ["N","Y","Effective reproduction number","Newly vaccinated persons"];
for i = 1:length(AlphaIndex)
    alphaT = flip(0:0.01*AlphaIndex(i)*2/(SimPeriod-1):(0.01*AlphaIndex(i)*2))';
    %     alphaT = 0.01*AlphaIndex(i)*ones(SimPeriod,1);
    if hconstant == 0
        ERNT = (AlphaIndexVariables(:,1,i)/POP0).*(((1 - h*alphaT).^k).*betaT)./(gammaT+deltaT);
    elseif hconstant == 1
        ERNT = (AlphaIndexVariables(:,1,i)/POP0).*(((1 - 0.01*h(1) - h(2)*alphaT).^k).*betaT)./(gammaT+deltaT);
    end
    ProjectedData = [AlphaIndexVariables(:,5,i),AlphaIndexVariables(:,6,i),ERNT,V];
    CombinedData = [Past;ProjectedData];
    figure(5)
    for j = 1:length(VarList2)
        subplot(2,2,j)
        plot(CombinedData(:,j),'k','LineWidth',2)
        xlim([1 Tdata+SimPeriod])
        title(VarName2(j),'FontWeight','normal')
        xline(Tdata);
        xticks([1 27 53 79 98])
        xticklabels( {'Jan-20','Jul-20','Jan-21','Jul-21','Dec-21'} )
        xtickangle(45)
        if j == 3
        ylim([0 3]);
        yline(1);
        end
        hold on
    end
    
    alpha2 = [alpha;alphaT];
    beta2 = [beta;betaT];
    beta_tildeT = ((1 - h*alphaT).^k).*betaT;
    beta_tilde2 = [beta_tilde;beta_tildeT];
    ERN2 = [ERN;ERNT];
    delta2 = [delta;deltaT];
    V2 = [zeros(Tdata,1);V];
    figure(101)
    %set(gcf,'Position',[100,100,800,500])
    for j = 1:length(ParamList2)
        subplot(2,2,j)
        plot(eval(ParamList2(j)),'k','LineWidth',2)
        xlim([1 Tdata+SimPeriod])
        title(ParamName(j),'Interpreter','latex','FontSize',11,'FontWeight','normal')
        xline(Tdata);
        xticks([1 27 53 79 98])
        xticklabels( {'Jan-20','Jul-20','Jan-21','Jul-21','Dec-21'} )
        xtickangle(45)
        hold on
    end
end

%--- Trade-off figure (baseline) ---%
figure(102)
plot(100*AverageAlpha,CumD,'k','LineWidth',2)
xlabel('Output Loss (%)','FontSize',fs)
ylabel('Cumlative Deaths','FontSize',fs)
xlim(xlim_tradeoff);
yline(D(end),'r--','LineWidth',1.5);
grid on
ax = gca;
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

%--- Trade-off figure (lag) ---%
figure(8)
plot(100*AverageAlpha,CumD,'k','LineWidth',2,'DisplayName','baseline')
legend('baseline')
hold on
for i = 1:length(wl)
    wlag = wl(i);
    plot(100*LagResults(2,:,i),LagResults(1,:,i),'LineWidth',2,'DisplayName',sprintf('%.0f weeks ago',wlag))
    hold on
end
hold off
xlabel('Output Loss (%)','FontSize',fs)
ylabel('Cumlative Deaths','FontSize',fs)
xlim(xlim_tradeoff);
ax = gca;
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
lgd = legend;
lgd.NumColumns = 1;
lgd.FontSize = 15;
grid on


%--- Trade-off figures with sensitivity analysis ---%
SensitivityTitle = ["\beta","h","k","Vaccine pace"];
SensitivitySaveTitle = ["beta","h","k","VacPace"];
SensitivityLegend = {["5% higher","baseline","5% lower"], ["20% lower", "baseline", "20% higher"], ["1.5","2 (baseline)","2.5"], ...
   ["0.5*baseline","baseline","2*baseline"]};
for i = 1:length(SensitivityTitle)
    figure(7)
    set(gcf,'Position',[100,100,600,500])
    subplot(2,2,i)
    for j = 1:3
        if j == 1
            plot(TargetAlpha,SA(j,:,i),'b--','LineWidth',1.5)
        elseif j == 2
            plot(TargetAlpha,SA(j,:,i),'k-','LineWidth',1.5)
        elseif j == 3
            plot(TargetAlpha,SA(j,:,i),'r-.','LineWidth',1.5)
        end
        hold on
    end
    xlabel('Output Loss (%)','FontSize',14)
    ylabel('Cumlative Deaths','FontSize',14)
    ax = gca;
    ax.YAxis.FontSize = 14;
    ax.XAxis.FontSize = 14;
    ax.YAxis.Exponent = 0;
    ytickformat('%,6.0f')
    xlim(xlim_tradeoff);
    legend(SensitivityLegend{i},'FontSize',13)
    title(SensitivityTitle(i),'FontSize',fs,'FontWeight','normal')
    grid on
    %saveas(figure(10+i),[home 'Figures/Paper/' char(SensitivitySaveTitle(i)) 'Sensitivity.png']);
    hold off;
end

%--- Counterfactual analysis for 2020 (multiplicative deviation of alpha)---%
sw = 2;
IV_CF = [S(sw),I(sw),R(sw),D(sw)];
AlphaDeviation2 = 0.5:0.01:1.3;
CF2 = zeros(2,length(AlphaDeviation2));
for i = 1:length(AlphaDeviation2)
    alpha_CF2 = alpha(sw:end-1)*AlphaDeviation2(i);
    [CF2(1,i),CF2(2,i)]=Covid_projection(IV_CF,alpha_CF2,beta(sw:end-1),gamma*ones(Tdata-sw,1),delta(sw:end-1),zeros(Tdata-sw,1),h,k,POP0,hconstant);
end
figure(9)
plot(100*CF2(2,:),CF2(1,:),'-ok','LineWidth',2,'MarkerIndices',find(AlphaDeviation2==1),'MarkerSize',15,'MarkerFaceColor','red','MarkerEdgeColor','none')
xlabel('Output Loss (%)','FontSize',18)
ylabel('Cumlative Deaths','FontSize',18)
xlim([3.5,5]);
grid on
ax = gca;
ax.YAxis.FontSize = 18;
ax.XAxis.FontSize = 18;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

slope = mean([(CF2(1,51)-CF2(1,50))/(100*(CF2(2,51)-CF2(2,50))),(CF2(1,52)-CF2(1,51))/(100*(CF2(2,52)-CF2(2,51)))]);
SVL = (POP0/slope)/100;

%%%%%%%%%%%%%%%%% Forecasting accuracy analysis %%%%%%%%%%%%%%%%%
StartDateF  = 33;
HorizonValsF=[1, 4];
dNForecast=zeros(Tdata,length(HorizonValsF));
dDForecast=zeros(Tdata,length(HorizonValsF));
dNActual=zeros(Tdata,length(HorizonValsF));
dDActual=zeros(Tdata,length(HorizonValsF));

%--- Regress mobility on alpha to estimate the elasticity h ---%

for iH = 1:length(HorizonValsF)
HorizonF  = HorizonValsF(iH);
EndtDateF = Tdata+1-HorizonF;

for iF = StartDateF:EndtDateF

GDP_F = GDP(1:iF-6);
referenceGDP_F = referenceGDP(1:iF-1);

M_F = M(1:iF-6,1);
X_F = 100*(1 - GDP_F./referenceGDP_F(1:iF-6));

RM_F = [M_F,X_F];
RM_F(any(isnan(RM_F), 2), :) = [];  % delete missing values
Y_F = RM_F(:,1);
X_F = RM_F(:,2);
h_F = -(X_F'*X_F)\X_F'*Y_F;

%--- Compute the history of S, I, R, D in the data period ---%
S_F = zeros(iF,1);
I_F = zeros(iF,1);
R_F = zeros(iF,1);
D_F = zeros(iF,1);
S_F(1)=POP0;
for i = 1:iF-1
    S_F(i+1)=S_F(i)-N(i);
    I_F(i+1)=I_F(i)+N(i)-gamma*I(i)-dD(i);
    R_F(i+1)=R_F(i)+gamma*I_F(i);
    D_F(i+1)=D_F(i)+dD(i); 
		for i = iF-5:iF-1
			GDP_F(i) = referenceGDP(i)*(1+0.5*mean(M(i-1:i))/(100*h_F));
		end
end
alpha_F = 1 - GDP_F./referenceGDP(1:iF-1);                                 
delta_F = (D_F(2:iF)-D_F(1:iF-1))./I_F(1:iF-1);                      
beta_tilde_F = -POP0*((S_F(2:iF)-S_F(1:iF-1))./(S_F(1:iF-1).*I_F(1:iF-1))); 
beta_F = beta_tilde_F./(1-h_F*alpha_F).^k;                                  

IV_F     = [S_F(iF),I_F(iF),R_F(iF),D_F(iF)];
alphaT_F = alpha(iF:iF+HorizonF-1); %Use actual values because we are interested in conditional forecasts.
beta_average_F = mean(beta_F(iF-RetroPeriod:iF-1));
delta_average_F = mean(delta_F(iF-1-RetroPeriod:iF-1));
betaT_F = beta_average_F*ones(HorizonF,1);
deltaT_F = delta_average_F*ones(HorizonF,1);
V_F = zeros(HorizonF,1);
gammaT_F=gamma*ones(HorizonF,1);

[CumD_F,AverageAlpha_F,SimData_F,SimN_F]=Covid_projection(IV_F,alphaT_F,betaT_F,gammaT_F,deltaT_F,V_F,h_F,k,POP0,hconstant);

dNForecast(iF,iH) = sum(SimN_F(1:HorizonF));
dDForecast(iF,iH) = CumD_F-IV_F(4);
dNActual(iF,iH)   = sum(N(iF:iF+HorizonF-1));
dDActual(iF,iH)   = D(iF+HorizonF-1)-D(iF-1);

end
end

xtick1_F=[StartDateF+4:8:Tdata Tdata];
xtick4_F=[StartDateF+4:8:Tdata Tdata];
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

%%%%%%%%%%%%%%%%% Uncertainty band for the tradeoff curve %%%%%%%%%%%%%%%%%

%--- Uncertainty-band parameters ---%
beta_se  = std(beta_sample)/sqrt(length(beta_sample));
delta_se= std(delta_sample)/sqrt(length(delta_sample));
%"h_se" is computed above.

BetaValsUB = beta_average+beta_se*[-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1];
DeltaValsUB = delta_average+delta_se*[-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1];
HValsUB = h+h_se*[1,0.75,0.5,0.25, 0, -0.25, -0.5, -0.75, -1];
UB = zeros(length(BetaValsUB),length(TargetAlpha));    % matrix to contain the results of uncertainty-band analysis

for i = 1:length(TargetAlpha)
    alphaT = flip(0:0.01*TargetAlpha(i)*2/(SimPeriod-1):(0.01*TargetAlpha(i)*2))';        
    for j = 1:length(BetaValsUB)
        betaS  = BetaValsUB(j)*ones(SimPeriod,1);
        deltaS = DeltaValsUB(j)*ones(SimPeriod,1);
        HS = HValsUB(:,j);
	      [UB(j,i)]=Covid_projection(InitialValues,alphaT,betaS,gammaT,deltaS,V,HS,k,POP0,hconstant);
    end
end

Y1UB = zeros(length(AverageAlpha),9);
X1UB = 100*AverageAlpha;
Y1UB(:,1)= UB(1,:);
for iUB=2:9
	Y1UB(:,iUB)= UB(iUB,:)-UB(iUB-1,:);
end

%--- Trade-off figure with UB (baseline) ---%
figure(6)
AreaUB=area(X1UB,Y1UB,0);
set(AreaUB(1),'FaceColor',[1 1 1])
set(AreaUB(2),'FaceColor',[0.9 0.9 0.9])
set(AreaUB(3),'FaceColor',[0.8 0.8 0.8])
set(AreaUB(4),'FaceColor',[0.7 0.7 0.7])
set(AreaUB(5),'FaceColor',[0.65 0.65 0.65])
set(AreaUB(6),'FaceColor',[0.65 0.65 0.65])
set(AreaUB(7),'FaceColor',[0.7 0.7 0.7])
set(AreaUB(8),'FaceColor',[0.8 0.8 0.8])
set(AreaUB(9),'FaceColor',[0.9 0.9 0.9])
hold on
for iUB=1:9
	set(AreaUB(iUB),'LineStyle','none')
end
plot(X1UB,CumD,'k','LineWidth',1.5)
xlabel('Output Loss (%)','FontSize',16)
ylabel('Cumlative Deaths','FontSize',16)
xlim([2,5]);
yline(D(end),'r--','LineWidth',1.5);
grid on
ax = gca;
ax.YAxis.FontSize = 16;
ax.XAxis.FontSize = 16;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

%--- Save all figures ---%
if figure_save == 1
    saveas(figure(1),[home 'Figures/NewCases.png']);
    saveas(figure(2),[home 'Figures/MobilityGDPLine.png']);
    saveas(figure(3),[home 'Figures/Parameters.png']);
    saveas(figure(4),[home 'Figures/Variables.png']);
    saveas(figure(5),[home 'Figures/VariablesProjection.png']);
    saveas(figure(6),[home 'Figures/BaselineTradeoffUB.png']);
    saveas(figure(7),[home 'Figures/Sensitivity.png']);
    saveas(figure(8),[home 'Figures/LaggedTradeoff.png']);
    saveas(figure(9),[home 'Figures/CF_multiplicative.png']);
    saveas(figure(10),[home 'Figures/ForecastErrorsN.png']);
    saveas(figure(11),[home 'Figures/ForecastErrorsD.png']);
    
    saveas(figure(100),[home 'Figures/MobilityGDPScatter.png']);
    saveas(figure(101),[home 'Figures/ParametersProjection.png']);
    saveas(figure(102),[home 'Figures/BaselineTradeoff.png']);
end





