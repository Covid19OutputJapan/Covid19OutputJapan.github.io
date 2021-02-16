
clear variables
close all
home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/';
cd(home);
%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 0;      % save back data
iPC=0;              % 1 if PC, 0 if Mac (or "home=pwd" does not work)
% in the "Figure" folder
fs = 16;            % common font size for many figures
%======================================================================%

%====================== Model parameter values ======================%
pref = 'Tokyo';        % prefecture to be analyzed
% GDP_vector = [106,40,36,24,22]*100;
prefGDP = 106;
gamma = 7/5;          % recovery rate from Covid
k = 2;                 % exponent of (1-h*alpha)
hconstant = 1;         % 0 = without intercept, 1 = with intercept for h regression
SimPeriod = 52;        % simulation period in weeks
VacStartDate = "Apr-01";             % time until the start of vaccination process

VacDuration = 12;       % time until vaccination pace reaches its steady pace
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
GDP = Data(:,5);
POP = Data(:,6);
Tdata= size(Data,1);    % Data period in weeks
POP0 = POP(1);          % initial population
xtick1 = 1:13:Tdata;
dateJP = string(datetime(dateD,'ConvertFrom','excel','Format','M?dd?'));
dateEN = string(datetime(dateD,'ConvertFrom','excel','Format','MMM-dd'));
Month = string(datetime(dateD,'ConvertFrom','excel','Format','MMM-yy'));
SimDate = dateD(end)+7:7:dateD(end)+7*SimPeriod;
SimDateEN = string(datetime(SimDate,'ConvertFrom','excel','Format','MMM-dd'));
%--- Create Month-Week labels ---%
dateP = dateD(end)+7:7:dateD(end)+7*(SimPeriod+1);
date = [dateD;dateP'];
MonthNumber = month(datetime(date,'ConvertFrom','excel'));
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
    MonthWeekJP(i) = [char(string(MonthNumber(i))) '??' char(string((WeekNumber(i)))) '?'];
    MonthWeekEN(i) = [char(string(datetime(date(i),'ConvertFrom','excel','Format','MMM'))) '-' char(string((WeekNumber(i)))) 'w'];
end

M = 1+0.01*M;
TdataGDP = Tdata-sum(isnan(GDP));
RetroH = TdataGDP-4;
VacStart = find(SimDateEN == "Apr-01");
End2020 = find(dateEN == "Jan-07");


%--- Constructing the reference level of output ---%
potentialGDP = zeros(52*3,1);
potentialGDP(1) = (100/(1.0122))*(1.0063^(1/12));
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
X = M(TdataGDP-17:TdataGDP);
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
Y = M(4:TdataGDP);
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

Y = M(TdataGDP-RetroH:TdataGDP);
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
    S(i+1)=S(i)-N(i);
    I(i+1)=I(i)+N(i)-gamma*I(i)-dD(i);
    R(i+1)=R(i)+gamma*I(i);
    D(i+1)=D(i)+dD(i);
    if i > TdataGDP
        GDP(i) = referenceGDP(i)*(1-alpha(i));
    end
end

%--- Compute the history of time-varying parameters ---%

delta = (D(2:Tdata+1)-D(1:Tdata))./I(1:Tdata);                              % death rate
beta_tilde = -POP0*((S(2:Tdata+1)-S(1:Tdata))./(S(1:Tdata).*I(1:Tdata)));   % overall infection rate
ERN = (S(1:end-1)/POP0).*beta_tilde./(gamma+delta);                                        % effective reproduction number
if hconstant == 0
    beta = beta_tilde./(1+h_all*alpha).^k;                                      % raw infection rate
elseif hconstant == 1
    %     beta = beta_tilde./(h(1)+h(2)*alpha).^k;
    beta = beta_tilde./(1+(h_all(2)/h_all(1))*alpha).^k;
end


%%%%%%%%%%%%%%%%% Projection starts here %%%%%%%%%%%%%%%%%
alpha_off = mean(alpha(Month == "Sep-20" | Month == "Oct-20"  | Month == "Nov-20")); % output loss without the state of emergency
InitialValues = [S(end),I(end),R(end),D(end)];

%--- Construct time series of parameters ---%
gammaT = gamma*ones(SimPeriod,1);
beta_sample = beta(end-RetroPeriod+1:end);
betaT = mean(beta_sample)*ones(SimPeriod,1);

delta_sample = delta(end-RetroPeriod+1:end);
delta_average = mean(delta_sample);
delta_ss = delta_average*(0.09/1.27);
deltaT = delta_average*ones(SimPeriod,1);
deltaT(VacStart+4:VacStart+4+8+18) = delta_average:(delta_ss-delta_average)/26:delta_ss;
deltaT(VacStart+4+8+18+1:end) = delta_ss;

VacPace = (POP0/125710000)*0.9*(3000000/2);  % number of vaccinations per week
V = zeros(SimPeriod,1);
V(VacStart-1:VacStart+VacDuration-1) = 0:VacPace/VacDuration:VacPace;
V(VacStart+VacDuration:end) = VacPace;


medical = 3000000*(POP0/125710000);
elderly = 30000000*(POP0/125710000);
elderly_total = 36000000*(POP0/125710000);
ordinary = 70000000*(POP0/125710000);

[V2,delta2,VT,real_pace] = vaccine_path(200000,SimPeriod,medical,elderly,ordinary,VacStart,VacStart+4,4,0.9,delta_average,elderly_total);

% 
% pace = 300000;
% VM = zeros(10000,6);        % Record vaccinations (first and second times) for each group
% pace_e = pace/2;
% pace_m = pace_e/2;
% 
% medical_vector = ones(ceil(medical/pace_m),1)*medical/ceil(medical/pace_m);
% elderly_vector = ones(ceil(elderly/pace_e),1)*elderly/ceil(elderly/pace_e);
% ordinary_vector = ones(ceil(ordinary/pace_e),1)*ordinary/ceil(ordinary/pace_e);
% ordinary_start = elderly_start + length(elderly_vector) + 1;
% real_pace = (elderly/ceil(elderly/pace_e))*2;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Projection parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th_on =1750*7;         % threshold to place the state of emergency (daily new infected persons in Tokyo = 1750)
ERN_on = 0.85;
altA = (((ERN_on*(POP0/S(end))*((gammaT(1)+deltaT(1))/betaT(1))).^(1/k))-1)*(h(1)/h(2));
% ERNCheck = (S(end)/POP0).*(((1+(h(2)/h(1))*altA).^k).*betaT(1))./(gammaT(1)+deltaT(1));
DMat = nan(3,20);
AlphaMat = nan(3,20);
AlphaPath = nan(SimPeriod,20,3);
NPath = nan(SimPeriod,20,3);

%---- 1. Different thresholds to lift the state of emergency ---%
TL = 100:50:1000;   
TL_index = [500,250];
for i = 1:length(TL)
    [DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData,NPath(:,i,1)] = Covid_projection_control_gradual(InitialValues,altA,alpha_off,th_on,TL(i)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,0);      
end
%---------------------------------------------------------------%

%---- 2. Different durations for economic recovery ---%
DR = 0:13;
DR_index = [0,6];
for i = 1:length(DR)
    [DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData,NPath(:,i,2)] = Covid_projection_control_gradual(InitialValues,altA,alpha_off,th_on,500*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,DR(i));      
end
%-----------------------------------------------------%

%------------ 3. Different vaccine pace -------------%