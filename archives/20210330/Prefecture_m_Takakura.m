% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
iPC=0;
if iPC==1
    %     home = '\Users\shcor\Dropbox\Website\Codes\';
    home = '\Users\takakurakazuma\Desktop\CovidRA\';
else
    % home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/';
    home = '/Users/takakurakazuma/Desktop/CovidRA/';
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


PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka','Aichi','Fukuoka'};
GDPVector = [106,36,23,21,40,41,20];
th_on_vector = [1250,700,400,350,500,350,350];
ERN_on_vector = [0.95,0.94,0.96,0.955,0.9,0.9,0.99];
TL_vector = {200:20:280,70:10:100,80:10:130,50:10:100,0:0.1:1,0:0.1:1,0:0.1:1};
TL_index_vector = {[280,240],[100,80],[130,110],[100,80],[100,40],[50,20],[60,10]};
DR_vector = {0:12,0:12,0:2:16,0:2:16,0:12,0:12,0:2:16};
rho_vector= [0.99,0.9,0.9,0,0.9,0.1,0];
impute_periods_vector =[45,22,25,20,45,30,30];
%rho_vector= [1.0,1.0,1.0,1.0,1.0,1.0,1.0];
%DR_index line 367
retro_ub = 17; 
retro_lb = 17;

beta_option = 3;

for pindex = 1:7 %length(PrefVector)
    close all
    
    %====================== Model parameter values ======================%
    pref = PrefVector{pindex};        % prefecture to be analyzed
    prefGDP = GDPVector(pindex);
    gamma = 7/5;          % recovery rate from Covid
    k = 2;                 % exponent of (1-h*alpha)
    hconstant = 1;         % 0 = without intercept, 1 = with intercept for h regression
    SimPeriod = 52;        % simulation period in weeks
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
        covid = importdata([home '\Covid_weekly_multi.csv']);  % Import weekly Covid data by prefecture
    else
        covid = importdata([home 'Covid_weekly_multi.csv']);  % Import weekly Covid data by prefecture
    end
    Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
    % Columns: 1 = date, 2 = new positive, 3 = death, 4 = mobility, 5=Credit, 6=POS, 7=cons_info, 8=restaurant_info, 9=joboffer 10 = FujiiGDP, 11=mixGDP 12 = population
    dateD = Data(:,1) + 21916;
    N = Data(:,2);
    dD = Data(:,3);
    M = Data(:,4);
    Credit = Data(:,5);
    POS = Data(:,6);
    consinfo = Data(:,7);
    restinfo = Data(:,8);
    job = Data(:,9);
    Fujii =Data(:,10);
    GDP = Data(:,11);
    POP = Data(:,12);
    %ERN_TK = Data(:,7);
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
    Credit = 1+0.01*Credit;
    POS = 1+0.01*POS;
    job = 1+0.01*job;
    consinfo = 1+0.01*consinfo;
    restinfo = 1+0.01*restinfo;
    Fujii = Fujii/100;
    TdataGDP = Tdata-sum(isnan(GDP));
    RetroH = TdataGDP-4;
    medical_start = find(SimDateEN == medical_start_date);
    elderly_start = find(SimDateEN == elderly_start_date);
    VacStart = find(SimDateEN == datetime(2021,4,1));
    End2020 = find(dateEN == datetime(2021,1,7));
    
    
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
    
    Mobility=M;
    Mobility(50)=0.5*(Mobility(49)+Mobility(51));
    alpha = (1 - GDP(1:TdataGDP)./referenceGDP(1:TdataGDP));   % output loss in percentage
    X = [Mobility(TdataGDP-impute_periods_vector(pindex):TdataGDP),Credit(TdataGDP-impute_periods_vector(pindex):TdataGDP),POS(TdataGDP-impute_periods_vector(pindex):TdataGDP),consinfo(TdataGDP-impute_periods_vector(pindex):TdataGDP),restinfo(TdataGDP-impute_periods_vector(pindex):TdataGDP),job(TdataGDP-impute_periods_vector(pindex):TdataGDP),Fujii(TdataGDP-impute_periods_vector(pindex):TdataGDP)];
    Y = alpha(TdataGDP-impute_periods_vector(pindex):TdataGDP);
    XC = [ones(length(X),1), X];
    s = (XC'*XC)\XC'*Y;         % OLS estimate of h with constant
    reg = XC*s;
    r_p = Y - reg; %residual
    SSE_p = sum(r_p.^2);
    eps_p = zeros(Tdata-TdataGDP,1);
    
    SSE_p
    
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
    
    alpha_pred = s(1)+s(2)*Mobility(TdataGDP+1:Tdata)+s(3)*Credit(TdataGDP+1:Tdata)+s(4)*POS(TdataGDP+1:Tdata)+s(5)*consinfo(TdataGDP+1:Tdata)+s(6)*restinfo(TdataGDP+1:Tdata)+s(7)*job(TdataGDP+1:Tdata)+s(8)*Fujii(TdataGDP+1:Tdata)+eps_p;
    
    alpha = [alpha;alpha_pred];
    
    
    %--- Plot mobility data ---%
    figure(pindex)
    title(pref)
    yyaxis left
    plot(Mobility,'k','LineWidth',1.5)
    hold on
    plot(Credit,'k','LineWidth',1.5)
    hold on
    plot(POS,'k','LineWidth',1.5)
    hold on
    plot(consinfo,'k','LineWidth',1.5)
    hold on
    plot(restinfo,'k','LineWidth',1.5)
    hold on
    plot(job,'k','LineWidth',1.5)
    hold on
    plot(Fujii,'k','LineWidth',1.5)
    hold on
    ylabel('Independent Variables')
    yyaxis right
    hold on
    plot((1-alpha)*100,'r-.','LineWidth',1.5)
    hold on
    plot((1-alpha(1:TdataGDP))*100,'b-.','LineWidth',1.5)
    ylabel('GDP')
    xlim([1 Tdata])
    xticks(xtick1)
    if iPC==1
        xticklabels(month(xtick1))
    else
        xticklabels(MonthWeekJP(xtick1))
    end
    legend('Mobility (L axis)','Credit Card (L axis)','POS (L axis)','Consumption Info (L axis)','Restaurant Info (L axis)','Job Offer (L axis)','Fujii (L axis)','Imputed GDP (R axis)', 'Data GDP (R axis)','FontSize',6,'Location','southeast');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'b';
    ax.YAxis(1).FontSize = fs;
    ax.YAxis(2).FontSize = fs;
    ax.XAxis.FontSize = fs;
    saveas(gcf,pref,'png');
   
end
