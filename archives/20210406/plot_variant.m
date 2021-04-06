% The share of variants should be correct.  
% However, this code does not consider the increase in betaT as in
% Prefecture_V_m.m, so the path of beta might be different (if
% betaT_temp_ini is not zero).

clear variables
close all
iPC=0;
if iPC==1
    %     home = '\Users\shcor\Dropbox\Website\Codes\';
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
    % home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/';
    home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/codes/';
end
cd(home);


%====================== Program parameter values ======================%
figure_save = 1;    % 0 = figures won't be saved, 1 = they will be saved
fs = 16;            % common font size for many figures
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
% medical_start_date = datetime(2021,3,18);
% elderly_start_date = datetime(2021,5,13);
%======================================================================%

PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka','Aichi','Fukuoka','Hyogo','Japan'};
GDPVector = [106,36,23,21,40,41,20,20,20];
th_on_vector = [1250,700,400,350,500,350,350,350,350];
ERN_now_vector = [0.99,0.94,0.96,0.955,0.9,0.9,0.99,0.99,0.99];
ERN_on_vector = [0.9,0.75,0.75,0.75,0.75,0.75,0.75,0.75];
ERN_on_scenario1 = 0.8;
ERN_on_scenario2 = 0.7;
TL_vector = {200:20:280,70:10:100,80:10:130,50:10:100,0:0.1:1,0:0.1:1,0:0.1:1,0:0.1:1,0:0.1:1};
TL_index_vector = {[280,240],[100,80],[130,110],[100,80],[100,40],[50,20],[60,10],[60,10],[60,10]};
DR_vector = {0:2:20,0:12,0:2:16,0:2:16,0:12,0:12,0:2:16,0:2:16,0:2:16};
rho_vector= [0.99,0.9,0.9,0,0.9,0.1,0,0,0];
impute_periods_vector =[45,22,25,20,45,30,30,30,30];
%rho_vector= [1.0,1.0,1.0,1.0,1.0,1.0,1.0];
%DR_index line 367
retro_ub = 17; 
retro_lb = 17;
%Parameters of variants
var_infection = 0.5;  % variant's infection rate (compared to nomral infection rate) [0.5, 0.7]
var_ss = 1.0;         % steady-state share of variant
var_growth_vec = [0.1695,0.5182];
% var_initial_vector = [0.000282,0.008786,0.007703,0.000395,0.000967,0,0,0.110650,0.011697]; 
% var_initial_vector = [0.000266,0.008669,0.008448,0.000419,0.000967,0,0,0.118087,0.011892]; %Up to 3/16
  var_initial_vector = [0.000685,0.009607,0.008621,0.008905,0.047151,0,0,0.133262,0.015965]; % Up to 3/23
% v_scale_vec = 0.125*ones(1,9);
  % v_scale_vec = [0.0270,0.0929,0.2053,0.0877,0.330,0.1,0.1,0.4404,0.1]; %Based on 3/1-3/7 reported values; Aichi and Fukuoka are set to 10%
  v_scale_vec = [0.03101,0.07192,0.15111,0.16345,0.19290,0.36209,0.38316,0.38657,0.21943]; %Based on 3/1-3/7 reported values; Aichi and Fukuoka are set to 10%
% Data Source: https://www.mhlw.go.jp/content/000755029.pdf
scaleA = 0.25; %Scenario A ... initial vector is multiplied by scaleA
scaleB = 0.5; %Scenario B ... initial vector is multiplied by scaleB
vec1 = [1,10,20,30,40];
vec2 = [1,2.5,4,6.5,8];
vec3 = [1,5,10,15,20];

pindex = 1;
pref = PrefVector{pindex};


beta_option = 41; %{1,3,41,42} 1 ... baseline, 3 ... 気の緩み, 
% 41 ... 変異株シナリオ1(アメリカのgrowht rate), 42 ... 変異株シナリオ2(イギリスのgrowth rate)
SimPeriod = 52;        % simulation period in weeks
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
% GDE = Data(:,7); %Gross Domestic Expenditure
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
TdataGDP = Tdata-sum(isnan(GDP));
RetroH = TdataGDP-4;
%     medical_start = find(SimDateEN == medical_start_date);
%     elderly_start = find(SimDateEN == elderly_start_date);
%     VacStart = find(SimDateEN == datetime(2021,4,1));
%     End2020 = find(dateEN == datetime(2021,1,7));
beta_avg = 1.9137;

var_initial_vec = zeros(length(vec1)+1,length(var_initial_vector));
var_initial_vec(1,:) = var_initial_vector ./ v_scale_vec;   
for i = 2:length(vec1)+1
    var_initial_vec(i,1:9) = var_initial_vector.*vec1(i-1);
    var_initial_vec(i,8) = var_initial_vector(8).*vec2(i-1);
    var_initial_vec(i,5) = var_initial_vector(5).*vec3(i-1);
end

var_share = zeros(length(SimDate),2,length(PrefVector),length(vec1)+1);
beta_bar = zeros(length(SimDate),2,length(PrefVector),length(vec1)+1);
betaT = zeros(length(SimDate),2,length(PrefVector),length(vec1)+1);
plot_var_share = zeros(length(SimDate)+1,2,length(PrefVector),length(vec1)+1);
plot_betaT = zeros(length(SimDate)+1,2,length(PrefVector),length(vec1)+1);
for i = 1:2
    for p = 1:1:length(PrefVector)
        for j = 1:length(vec1)+1
            if p == 1
            %             var_growth = 0.03;     % weekly growth parameter for logit model [obtained from Kusaka analysis]
                FirstObsDate = datetime(2020,12,31);
            elseif p == 2
                FirstObsDate = datetime(2021,2,11);
            elseif p == 3
            %             var_growth = 0.24;  % weekly growth parameter for logit model [obtained from kusaka]
                FirstObsDate = datetime(2021,2,4);
            elseif p == 4
                FirstObsDate = datetime(2021,3,11);
            elseif p == 5
                FirstObsDate = datetime(2021,2,25);
            end    
            FirstObsTime = find(dateEN == FirstObsDate);
            var_initial = var_initial_vec(j,p);
            if i == 1 && j == 1
                var_initial = scaleA*var_initial_vec(j,p);
            end
            if i == 2 && j == 1
                var_initial = scaleB*var_initial_vec(j,p);
            end
            logit_initial(i) = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
            %         var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept
            var_intercept(i) = logit_initial(i);
            var_growth = var_growth_vec(1,i);
            var_share(:,i,p,j) = exp((1:SimPeriod)'*var_growth+var_intercept(i))./(1+exp((1:SimPeriod)'*var_growth+var_intercept(i)));
            beta_bar(:,i,p,j) = beta_avg/(1+var_infection*var_initial);
            betaT(:,i,p,j) = beta_bar(:,i,p,j).*(1+var_infection.*var_share(:,i,p,j));
            plot_var_share(1,i,p,j) = var_initial;
            plot_var_share(2:length(SimDate)+1,i,p,j) = var_share(:,i,p,j);
            plot_betaT(1,i,p,j) = beta_avg;
            plot_betaT(2:length(SimDate)+1,i,p,j) = betaT(:,i,p,j);
        end
    end
%     if p == 8
%         if p == 8
%             FirstObsDate = datetime(2021,1,14);
%         end
%         FirstObsTime = find(dateEN == FirstObsDate);
%         var_initial = var_initial_vec(1,p);
%         logit_initial = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
%         %         var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept
%         var_intercept = logit_initial;
%         var_growth = var_growth_vec(1,i);
%         var_share(:,i,p) = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
%         beta_bar(:,i,p) = beta_avg/(1+var_infection*var_initial);
%         betaT(:,i,p) = beta_bar(:,i,p).*(1+var_infection.*var_share(:,i,p));
%         plot_var_share(1,i,p) = var_initial;
%         plot_var_share(2:length(SimDate)+1,i,p) = var_share(:,i,p);
%         plot_betaT(1,i,p) = beta_avg;
%         plot_betaT(2:length(SimDate)+1,i,p) = betaT(:,i,p);
%     end
end

xaxis_vec = 0:1:length(SimDate);

MonthWeekEN_Sim = MonthWeekEN(length(dateD)+1:end,1);
MonthNumber_Sim = MonthNumber(length(dateD)+1:end,1);
WeekNumber_Sim = WeekNumber(length(dateD)+1:end,1);
MonthWeekJP_Sim = MonthWeekJP(length(dateD)+1:end,1);

for l = 1:2
    if l == 1
        figure(4101)
    elseif l == 2
        figure(4102)
    end 
    set(gcf,'Position',[100,100,1400,800])
    
    subplot(2,3,1)
    plot(xaxis_vec,plot_var_share(:,1,1,1)','-b','LineWidth',2)
    hold on
%     for j = 2:length(vec1)+1
%         plot(xaxis_vec,plot_var_share(:,1,1,j)','--b','LineWidth',0.3,'DisplayName',sprintf('%.0f',vec1(j-1)))
%     end
    plot(xaxis_vec,plot_var_share(:,2,1,1)','-r','LineWidth',2)
%     for j = 2:length(vec1)+1
%         plot(xaxis_vec,plot_var_share(:,2,1,j)','--r','LineWidth',0.3,'DisplayName',sprintf('%.0f',vec1(j-1)))
%     end
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    ax.YAxis.Exponent = 0;
    ytickformat('%,0.2f')
    
    xticks(find(WeekNumber_Sim==1))
    if l == 1
        legend('ScenarioA','ScenarioB','FontSize',10)
        title('TOKYO','FontSize',20,'FontWeight','normal')
        xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
        xlabel('Weeks','FontSize',20)
        ylabel('Share of Variants','FontSize',20)
    elseif l == 2
        legend('シナリオA','シナリオB','FontName',fn,'FontSize',10)
        title('東京','FontSize',20,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
        xlabel('週','FontSize',20,'FontName',fn)
        ylabel('変異株割合','FontSize',20,'FontName',fn)
    end
    xtickangle(45)
    xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
    xlim([0,length(SimDate)]);
    ylim([0,1]);
    lgd = legend;
    lgd.Location = 'southeast';

    subplot(2,3,2)
        plot(xaxis_vec,plot_var_share(:,1,2)','-b','LineWidth',2)
    hold on
    plot(xaxis_vec,plot_var_share(:,2,2)','-r','LineWidth',2)
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    ax.YAxis.Exponent = 0;
    ytickformat('%,0.2f')
    xticks(find(WeekNumber_Sim==1))
    if l == 1
        legend('ScenarioA','ScenarioB','FontSize',10)
        title('KANAGAWA','FontSize',20,'FontWeight','normal')
        xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
        xlabel('Weeks','FontSize',20)
        ylabel('Share of Variants','FontSize',20)
    elseif l == 2
        legend('シナリオA','シナリオB','FontName',fn,'FontSize',10)
        title('神奈川','FontSize',20,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
        xlabel('週','FontSize',20,'FontName',fn)
        ylabel('変異株割合','FontSize',20,'FontName',fn)
    end
    xtickangle(45)
    xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
    xlim([0,length(SimDate)]);
    ylim([0,1]);
    lgd = legend;
    lgd.Location = 'southeast';
    
    
    subplot(2,3,3)
    plot(xaxis_vec,plot_var_share(:,1,3)','-b','LineWidth',2)
    hold on
    plot(xaxis_vec,plot_var_share(:,2,3)','-r','LineWidth',2)
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    ax.YAxis.Exponent = 0;
    ytickformat('%,0.2f')
    xticks(find(WeekNumber_Sim==1))
    if l == 1
        legend('ScenarioA','ScenarioB','FontSize',10)
        title('SAITAMA','FontSize',20,'FontWeight','normal')
        xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
        xlabel('Weeks','FontSize',20)
        ylabel('Share of Variants','FontSize',20)
    elseif l == 2
        legend('シナリオA','シナリオB','FontName',fn,'FontSize',10)
        title('埼玉','FontSize',20,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
        xlabel('週','FontSize',20,'FontName',fn)
        ylabel('変異株割合','FontSize',20,'FontName',fn)
    end
    xtickangle(45)
    xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
    xlim([0,length(SimDate)]);
    ylim([0,1]);
    lgd = legend;
    lgd.Location = 'southeast';
    
    
    
    subplot(2,3,4)
    plot(xaxis_vec,plot_var_share(:,1,4)','-b','LineWidth',2)
    hold on
    plot(xaxis_vec,plot_var_share(:,2,4)','-r','LineWidth',2)
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    ax.YAxis.Exponent = 0;
    ytickformat('%,0.2f')
    xticks(find(WeekNumber_Sim==1))
    if l == 1
        legend('ScenarioA','ScenarioB','FontSize',10)
        title('CHIBA','FontSize',20,'FontWeight','normal')
        xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
        xlabel('Weeks','FontSize',20)
        ylabel('Share of Variants','FontSize',20)
    elseif l == 2
        legend('シナリオA','シナリオB','FontName',fn,'FontSize',10)
        title('千葉','FontSize',20,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
        xlabel('週','FontSize',20,'FontName',fn)
        ylabel('変異株割合','FontSize',20,'FontName',fn)
    end
    xtickangle(45)
    xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
    xlim([0,length(SimDate)]);
    ylim([0,1]);
    lgd = legend;
    lgd.Location = 'southeast';
    
    
%     subplot(2,3,5)
%     plot(xaxis_vec,plot_var_share(:,1,8)','-b','LineWidth',2)
%     hold on
%     plot(xaxis_vec,plot_var_share(:,2,8)','-r','LineWidth',2)
%     ax = gca;
%     ax.YAxis.FontSize = 20;
%     ax.XAxis.FontSize = 20;
%     ax.YAxis.Exponent = 0;
%     ytickformat('%,0.2f')
%     xticks(find(WeekNumber_Sim==1))
%     if l == 1
%         legend('ScenarioA','ScenarioB','FontSize',10)
%         title('HYOGO','FontSize',20,'FontWeight','normal')
%         xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
%         xlabel('Weeks','FontSize',20)
%         ylabel('Share of Variants','FontSize',20)
%     elseif l == 2
%         legend('シナリオA','シナリオB','FontName',fn,'FontSize',10)
%         title('兵庫','FontSize',20,'FontWeight','normal','FontName',fn)
%         xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
%         xlabel('週','FontSize',20,'FontName',fn)
%         ylabel('変異株割合','FontSize',20,'FontName',fn)
%     end
%     xtickangle(45)
%     xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
%     xlim([0,length(SimDate)]);
%     ylim([0,1]);
%     lgd = legend;
%     lgd.Location = 'southeast';
%     if l == 1 
%         sgt = sgtitle('Projected path of variant share');
%         sgt.FontSize = fs*2;
%     elseif l == 2
%         sgt = sgtitle('変異株シェアの推移');
%         sgt.FontSize = fs*2;
%         sgt.FontName = fn;
%     end 
% 
    
    subplot(2,3,5)
    plot(xaxis_vec,plot_var_share(:,1,5)','-b','LineWidth',2)
    hold on
    plot(xaxis_vec,plot_var_share(:,2,5)','-r','LineWidth',2)
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    ax.YAxis.Exponent = 0;
    ytickformat('%,0.2f')
    xticks(find(WeekNumber_Sim==1))
    if l == 1
        legend('ScenarioA','ScenarioB','FontSize',10)
        title('OSAKA','FontSize',20,'FontWeight','normal')
        xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
        xlabel('Weeks','FontSize',20)
        ylabel('Share of Variants','FontSize',20)
    elseif l == 2
        legend('シナリオA','シナリオB','FontName',fn,'FontSize',10)
        title('大阪','FontSize',20,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
        xlabel('週','FontSize',20,'FontName',fn)
        ylabel('変異株割合','FontSize',20,'FontName',fn)
    end
    xtickangle(45)
    xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
    xlim([0,length(SimDate)]);
    ylim([0,1]);
    lgd = legend;
    lgd.Location = 'southeast';
    if l == 1 
        sgt = sgtitle('Projected path of variant share');
        sgt.FontSize = fs*2;
    elseif l == 2
        sgt = sgtitle('変異株シェアの推移');
        sgt.FontSize = fs*2;
        sgt.FontName = fn;
    end 



    subplot(2,3,6)
    plot(xaxis_vec,plot_var_share(:,1,9)','-b','LineWidth',2)
    hold on
    plot(xaxis_vec,plot_var_share(:,2,9)','-r','LineWidth',2)
    ax = gca;
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    ax.YAxis.Exponent = 0;
    ytickformat('%,0.2f')
    xticks(find(WeekNumber_Sim==1))
    if l == 1
        legend('Scenario1','Scenario2','FontSize',10)
        title('Japan','FontSize',20,'FontWeight','normal')
        xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
        xlabel('Weeks','FontSize',20)
        ylabel('Share of Variants','FontSize',20)
    elseif l == 2
        legend('シナリオA','シナリオB','FontName',fn,'FontSize',10)
        title('全国','FontSize',20,'FontWeight','normal','FontName',fn)
        xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
        xlabel('週','FontSize',20,'FontName',fn)
        ylabel('変異株割合','FontSize',20,'FontName',fn)
    end
    xtickangle(45)
    xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
    xlim([0,length(SimDate)]);
    ylim([0,1]);
    lgd = legend;
    lgd.Location = 'southeast';
    if l == 1 
        sgt = sgtitle('Projected path of variant share');
        sgt.FontSize = fs*2;
    elseif l == 2
        sgt = sgtitle('変異株シェアの推移');
        sgt.FontSize = fs*2;
        sgt.FontName = fn;
    end 

    
end 


if l == 1
    figure(4201)
elseif l == 2
    figure(4202)
end 
set(gcf,'Position',[100,100,1400,800])
    
plot(xaxis_vec,plot_var_share(:,1,pindex,1)','-b','LineWidth',2,'DisplayName',['シナリオA,scale =' sprintf('%.2f',scaleA*1/v_scale_vec(pindex))])
hold on
if pindex == 5
    plot(xaxis_vec,plot_var_share(:,1,pindex,2)','--','LineWidth',0.5,'Color','#33FFFF','DisplayName',['scale = ' sprintf('%.0f',vec3(1))])
    plot(xaxis_vec,plot_var_share(:,1,pindex,3)','--','LineWidth',0.5,'Color','#33CCFF','DisplayName',['scale = ' sprintf('%.0f',vec3(2))])
    plot(xaxis_vec,plot_var_share(:,1,pindex,4)','--','LineWidth',0.5,'Color','#33AAFF','DisplayName',['scale = ' sprintf('%.0f',vec3(3))])
    plot(xaxis_vec,plot_var_share(:,1,pindex,5)','--','LineWidth',0.5,'Color','#3366FF','DisplayName',['scale = ' sprintf('%.0f',vec3(4))])
    plot(xaxis_vec,plot_var_share(:,1,pindex,6)','--','LineWidth',0.5,'Color','#3333FF','DisplayName',['scale = ' sprintf('%.0f',vec3(5))])
    
plot(xaxis_vec,plot_var_share(:,2,pindex,1)','-r','LineWidth',2,'DisplayName',['シナリオB,scale =' sprintf('%.2f',scaleB/v_scale_vec(pindex))])

    plot(xaxis_vec,plot_var_share(:,2,pindex,2)','--','LineWidth',0.5,'Color','#FF66FF','DisplayName',['scale = ' sprintf('%.0f',vec3(1))])
    plot(xaxis_vec,plot_var_share(:,2,pindex,3)','--','LineWidth',0.5,'Color','#FF55CC','DisplayName',['scale = ' sprintf('%.0f',vec3(2))])
    plot(xaxis_vec,plot_var_share(:,2,pindex,4)','--','LineWidth',0.5,'Color','#FF44AA','DisplayName',['scale = ' sprintf('%.0f',vec3(3))])
    plot(xaxis_vec,plot_var_share(:,2,pindex,5)','--','LineWidth',0.5,'Color','#FF3366','DisplayName',['scale = ' sprintf('%.0f',vec3(4))])
    plot(xaxis_vec,plot_var_share(:,2,pindex,6)','--','LineWidth',0.5,'Color','#FF2233','DisplayName',['scale = ' sprintf('%.0f',vec3(5))])

else
    plot(xaxis_vec,plot_var_share(:,1,pindex,2)','--','LineWidth',0.5,'Color','#33FFFF','DisplayName',['scale = ' sprintf('%.0f',vec1(1))])
    plot(xaxis_vec,plot_var_share(:,1,pindex,3)','--','LineWidth',0.5,'Color','#33CCFF','DisplayName',['scale = ' sprintf('%.0f',vec1(2))])
    plot(xaxis_vec,plot_var_share(:,1,pindex,4)','--','LineWidth',0.5,'Color','#33AAFF','DisplayName',['scale = ' sprintf('%.0f',vec1(3))])
    plot(xaxis_vec,plot_var_share(:,1,pindex,5)','--','LineWidth',0.5,'Color','#3366FF','DisplayName',['scale = ' sprintf('%.0f',vec1(4))])
    plot(xaxis_vec,plot_var_share(:,1,pindex,6)','--','LineWidth',0.5,'Color','#3333FF','DisplayName',['scale = ' sprintf('%.0f',vec1(5))])
    
plot(xaxis_vec,plot_var_share(:,2,pindex,1)','-r','LineWidth',2,'DisplayName',['シナリオB,baseline,scale =' sprintf('%.2f',scaleB/v_scale_vec(pindex))])

    plot(xaxis_vec,plot_var_share(:,2,pindex,2)','--','LineWidth',0.5,'Color','#FF66FF','DisplayName',['scale = ' sprintf('%.0f',vec1(1))])
    plot(xaxis_vec,plot_var_share(:,2,pindex,3)','--','LineWidth',0.5,'Color','#FF55CC','DisplayName',['scale = ' sprintf('%.0f',vec1(2))])
    plot(xaxis_vec,plot_var_share(:,2,pindex,4)','--','LineWidth',0.5,'Color','#FF44AA','DisplayName',['scale = ' sprintf('%.0f',vec1(3))])
    plot(xaxis_vec,plot_var_share(:,2,pindex,5)','--','LineWidth',0.5,'Color','#FF3366','DisplayName',['scale = ' sprintf('%.0f',vec1(4))])
    plot(xaxis_vec,plot_var_share(:,2,pindex,6)','--','LineWidth',0.5,'Color','#FF2233','DisplayName',['scale = ' sprintf('%.0f',vec1(5))])
end

ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
ax.YAxis.Exponent = 0;
ytickformat('%,0.2f')

xticks(find(WeekNumber_Sim==1))
if l == 1
    % legend('Variant Share (Scenario1)','Variant Share (Scenario2)','FontSize',10)
    title(pref,'FontSize',20,'FontWeight','normal')
    xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
    xlabel('Weeks','FontSize',20)
    ylabel('Share of Variants','FontSize',20)
elseif l == 2
    % legend('変異株割合 (シナリオ1),baseline','変異株割合 (シナリオ2)','FontName',fn,'FontSize',10)
    title(pref,'FontSize',20,'FontWeight','normal','FontName',fn)
    xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
    xlabel('週','FontSize',20,'FontName',fn)
    ylabel('変異株割合','FontSize',20,'FontName',fn)
end
xtickangle(45)
xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
xlim([0,length(SimDate)]);
ylim([0,1]);
lgd = legend;
lgd.Location = 'southeast';
lgd.FontName = fn;
lgd.FontSize = 10;
lgd.NumColumns = 2;

if figure_save == 1
    saveas(figure(4101),[home 'Figures/' char(pref) '/VariantBetaPath.png']);
    saveas(figure(4102),[home 'Figures/' char(pref) '/VariantBetaPath_jp.png']);
    saveas(figure(4202),[home 'Figures/' char(pref) '/VariantSharePath_jp.png']);

end
    

%     subplot(1,2,2)
%     plot(xaxis_vec, plot_betaT(:,1),'-b','LineWidth',2)
%     hold on
%     plot(beta_avg*ones(length(SimDate)+1,1),'-k','LineWidth',2)
%     ax = gca;
%     ax.YAxis.FontSize = 20;
%     ax.XAxis.FontSize = 20;
%     ax.YAxis.Exponent = 0;
%     ytickformat('%,0.2f')
%     xticks(find(WeekNumber_Sim==1))
%     xtickangle(45)
%     if l == 1
%         legend('Raw Infection Rate with Variants', 'Raw Infection Rate w/o Variants','FontSize',10)
%         title('Projected path of beta','FontSize',20,'FontWeight','normal')
%         xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
%         xlabel('Weeks','FontSize',20)
%         ylabel('Infection Rate','FontSize',20,'FontName',fn)
%     elseif l == 2
%         legend('感染率(変異株成長ケース)','現在の感染率','FontSize',10,'FontName',fn)
%         xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
%         title('βの推移','FontSize',20,'FontWeight','normal','FontName',fn)
%         xlabel('週','FontSize',20)
%         ylabel('感染率','FontSize',20,'FontName',fn)
%     end
%     xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
%     xlim([0,length(SimDate)]);
%     lgd = legend;
%     lgd.Location = 'northwest';
% 
