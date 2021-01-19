% This m-file executes simulation and generates figures for the 
% analysis of the state of emergency placed in Tokyo in the paper 
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
home = '/home/takeki/Dropbox/covid19outputjapan/code/';
%home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Covid19_Output/';
% home = pwd;
%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
% in the "Figure" folder
fs = 16;            % common font size for many figures
xlim_tradeoff = [2,5];
%======================================================================%

%====================== Model parameter values ======================%
pref = 'Tokyo';        % prefecture to be analyzed
gamma = 7/10;          % recovery rate from Covid
k = 2;                 % exponent of (1-h*alpha)
hconstant = 0;         % 0 = without intercept, 1 = with intercept for h regression
SimPeriod = 52;        % simulation period in weeks
VacStart = 12;         % time until the start of vaccination process
VacPace = 0.8*(36000*7)/2;  % number of vaccinations per week
RetroPeriod = 13;      % retroactive periods used to estimate gamma and delta
alpha_on = 4.2;        %3.75 = 12 weeks, 4.2 = 8 weeks, 5.5 = 4 weeks
% alpha_off = 100*mean(alpha(Month == "Sep-20" | Month == "Oct-20"  | Month == "Nov-20")); % output loss without the state of emergency
th_on = 14000;         % threshold to place the state of emergency (daily new infected persons in Tokyo = 2000)
th_off = 3500;         % threshold to remove the state of emergency (daily new infected persons in Tokyo = 500)
target_duration = 8;
%====================================================================%


%--- Import data ---%
% Covid data are recorded at weekly frequencies (Mon-Sun)
% The first week start on January 20 (Mon), 2020
covid = importdata([home 'Data/Covid_weekly_prefecture.csv']);  % Import weekly Covid data by prefecture
Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
% Columns: 1 = date, 2 = new positive, 3 = death, 4 = mobility, 5 = GDP, 6 = population
date = Data(:,1) + 21916;
N = Data(:,2);
dD = Data(:,3);
M = Data(:,4);
GDP = Data(:,5);
POP = Data(:,6);
Tdata= size(Data,1);    % Data period in weeks
POP0 = POP(1);          % initial population
xtick1 = [1:13:Tdata, Tdata];
dateJP = string(datetime(date,'ConvertFrom','excel','Format','M? dd?'));
dateEN = string(datetime(date,'ConvertFrom','excel','Format','MMM-dd'));
Month = string(datetime(date,'ConvertFrom','excel','Format','MMM-yy'));

%--- Constructing the reference level of output ---%
potentialGDP = zeros(52*3,1);
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

%--- Regress mobility on alpha to estimate the elasticity h ---%
X = 100*(1 - GDP./referenceGDP(1:Tdata));   % output loss in percentage
RM = [M,X];
RM(any(isnan(RM), 2), :) = [];  % delete missing values
Y = RM(:,1);
X = RM(:,2);
h = -(X'*X)\X'*Y;              % OLS estimate of h
reg = -h*X;

% Use data in recent five months %
RM = RM(end-21:end,:);
Y = RM(:,1);
X = RM(:,2);
h2 = -(X'*X)\X'*Y;  

%--- Compute the history of S, I, R, D in the data period ---%
S = zeros(Tdata+1,1);
I = zeros(Tdata+1,1);
R = zeros(Tdata+1,1);
D = zeros(Tdata+1,1);
S(1)=POP0;
for i = 1:Tdata
    S(i+1)=S(i)-N(i);
    I(i+1)=I(i)+N(i)-gamma*I(i)-dD(i);
    R(i+1)=R(i)+gamma*I(i);
    D(i+1)=D(i)+dD(i); 
    if isnan(GDP(i)) == 0
        lastGDP = i;
    else
        GDP(i) = 0.5*GDP(lastGDP)+0.5*referenceGDP(i)*(1+M(i)/(100*h2));
    end
end

%--- Compute the history of time-varying parameters ---%
alpha = 1 - GDP./referenceGDP(1:Tdata);                                          % output loss
delta = (D(2:Tdata+1)-D(1:Tdata))./I(1:Tdata);                          % death rate
beta_tilde = -POP0*((S(2:Tdata+1)-S(1:Tdata))./(S(1:Tdata).*I(1:Tdata)));      % overall infection rate
beta = beta_tilde./(1-h*alpha).^k;                                      % raw infection rate
BRN = beta_tilde./(gamma+delta);                                        % basic reproduction number

alpha_off = 100*mean(alpha(Month == "Sep-20" | Month == "Oct-20"  | Month == "Nov-20")); % output loss without the state of emergency


%%%%%%%%%%%%%%%%% Projection starts here %%%%%%%%%%%%%%%%%

InitialValues = [S(end),I(end),R(end),D(end)];
dateP = date(end)+7:7:date(end)+7*(SimPeriod+1);
dateP = string(datetime(dateP,'ConvertFrom','excel','Format','MMM'));

%--- Construct time series of parameters ---%
beta_sample = beta(end-RetroPeriod+1:end);
betaT = mean(beta_sample)*ones(SimPeriod,1);
delta_sample = delta(end-RetroPeriod+1:end);
deltaT = mean(delta_sample)*ones(SimPeriod,1);
gammaT = gamma*ones(SimPeriod,1);
V = zeros(SimPeriod,1);
V(VacStart-1:VacStart+3) = 0:VacPace/4:VacPace;
V(VacStart+4:end) = VacPace;

%--- Projection for different th_off values ---%
TH = (100:50:1000)*7;
TH_index = 500*7;
DM = zeros(1,length(TH));
AlphaM = zeros(1,length(TH));
AlphaPath = zeros(SimPeriod,length(TH));
NPath = zeros(SimPeriod,length(TH));

% clf
% altA = [3.65, 4, 5];
% for y = 1:length(altA)
%     alpha_on = altA(y);

figure(1)
set(gcf,'Position',[100,100,1200,500])
subplot(1,2,1)
for i = 1:length(TH)
    [DM(i),AlphaM(i),AlphaPath(:,i),SimData,NPath(:,i)] = Covid_projection_control(InitialValues,alpha_on,alpha_off,th_on,TH(i),betaT,gammaT,deltaT,V,h2,k,POP0);
    if sum(TH(i)==TH_index)==1
        plot(NPath(:,i),'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)/7))
    else
        plot(NPath(:,i),'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)/7))
    end
    hold on
end
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
ylabel('Number of new cases per week','FontSize',20)
title('Projected path of new cases','FontSize',20,'FontWeight','normal')
% lgd = legend;
% lgd.NumColumns = 2;
xticks([1 9 17 25 33 41 49])
xticklabels(dateP([1 9 17 25 33 41 49]))
xlim([1 SimPeriod])
% xticks([1 14 27 40 52])
% xticklabels(dateP([1 14 27 40 52]))
xtickangle(45)

%--- Trade-off figures ---%
waves = zeros(1,length(TH));
for i = 1:length(TH)
    svec = zeros(SimPeriod-1,1);
    for t = 1:SimPeriod-1
        svec(t) = AlphaPath(t+1,i)-AlphaPath(t,i);
    end
    waves(i) = nnz(svec);
end

subplot(1,2,2)
plot(AlphaM(waves==1),DM(waves==1),'-bo','LineWidth',1.5,'MarkerSize',10)
hold on
plot(AlphaM(waves==3),DM(waves==3),'-go','LineWidth',1.5,'MarkerSize',10)
hold on
plot(AlphaM(waves==5),DM(waves==5),'-mo','LineWidth',1.5,'MarkerSize',10)
hold on
text(AlphaM,DM,string(TH/7),'VerticalAlignment','bottom','HorizontalAlignment','left')
hold on
scatter(AlphaM(TH==500*7),DM(TH==500*7),80,'red','filled');
xlabel('Output Loss (%)','FontSize',20)
ylabel('Cumlative Death','FontSize',20)
title('Relationsihp between Covid-19 and output','FontSize',20,'FontWeight','normal')
grid on
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

% NM = sum(NPath,1)+sum(N);
% subplot(1,2,2)
% plot(AlphaM(waves==1),NM(waves==1),'-bo','LineWidth',1.5,'MarkerSize',10)
% hold on
% plot(AlphaM(waves==3),NM(waves==3),'-go','LineWidth',1.5,'MarkerSize',10)
% hold on
% plot(AlphaM(waves==5),NM(waves==5),'-mo','LineWidth',1.5,'MarkerSize',10)
% hold on
% text(AlphaM,NM,string(TH/7),'VerticalAlignment','bottom','HorizontalAlignment','left')
% hold on
% scatter(AlphaM(TH==500*7),NM(TH==500*7),80,'red','filled');
% xlabel('Output Loss (%)','FontSize',20)
% ylabel('Cumlative number of infected people','FontSize',20)
% title('Relationsihp between Covid-19 and output','FontSize',20,'FontWeight','normal')
% grid on
% ax = gca;
% ax.YAxis.FontSize = 20;
% ax.XAxis.FontSize = 20;
% ax.YAxis.Exponent = 0;
% ytickformat('%,6.0f')

% saveas(figure(1),[home 'Figures/Lockdown_Exit/Alpha' sprintf('%.2f',alpha_on) '.png']);
% saveas(figure(1),[home 'Figures/Lockdown_Exit/Alpha' sprintf('%.2f',alpha_on) 'CumN.png']);

% end
