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
    % home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/';
    home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/';
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 0;      % save back data
% in the "Figure" folder
fs = 16;            % common font size for many figures
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
%======================================================================%


PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka','Aichi','Fukuoka','Hyogo'};
GDPVector = [106,36,23,21,40,41,20,20];
th_on_vector = [1250,700,400,350,500,350,350,350];
ERN_now_vector = [0.95,0.94,0.96,0.955,0.9,0.9,0.99,0.99];
ERN_on_vector = [0.95,0.75,0.75,0.75,0.75,0.75,0.75];
ERN_on_scenario1 = 0.8; %0.8
ERN_on_scenario2 = 0.7; %0.7
TL_vector = {280:10:350,70:10:100,80:10:130,50:10:100,0:0.1:1,0:0.1:1,0:0.1:1,0:0.1:1};
TL_index_vector = {[280,350],[100,80],[130,110],[100,80],[100,40],[50,20],[60,10],[60,10]}; 
% TL_vector = {240:10:300,70:10:100,80:10:130,50:10:100,0:0.1:1,0:0.1:1,0:0.1:1,0:0.1:1};
% TL_index_vector = {[240,300],[100,80],[130,110],[100,80],[100,40],[50,20],[60,10],[60,10]};
th_off_vector = [280,80,110,80,40,20,10,10];
DR_vector = {0:2:16,0:12,0:2:16,0:2:16,0:12,0:12,0:2:16,0:2:16};
rho_vector= [0.8,0.9,0.9,0,0.9,0.1,0,0];
impute_periods_vector =[45,22,25,20,45,30,30,30];
%rho_vector= [1.0,1.0,1.0,1.0,1.0,1.0,1.0];
%DR_index line 367
retro_ub = 17; 
retro_lb = 17;
betaT_temp_ini = 0.0;
% 0.2 ... 35人, 0.3 ... 60人, 0.25 ... 50人, 0.185 & 0.5 ... 30人
beta_rho = 0.6;
% 0.3, 0.5 のときTL_vector_index [280,420]
% 0.25, 0.5 のときTL_vector_index [280,380]
% 0.2, 0.6 のときTL_vector_index [280,360]

%Parameters of variants
var_infection = 0.5;  % variant's infection rate (compared to nomral infection rate) [0.5, 0.7]
var_ss = 1.0;         % steady-state share of variant
% v_scale = 10;
var_growth_vec = [0.1695,0.5182];
var_initial_vector = [0.000282,0.008786,0.007703,0.000395,0.000967,0,0,0.110650,0.011697]; % in percentage
% v_scale_vec = 0.125*ones(1,9);
v_scale_vec = [0.0270,0.0929,0.2053,0.0877,0.330,0.1,0.1,0.4404,0.1]; %Based on 3/1-3/7 reported values; Aichi and Fukuoka are set to 10%
vec1 = [1,5,10,20,50];
vec2 = [1,2.5,4,6.5,8];



beta_option = 1; %{1,3,41,42} 1 ... baseline, 3 ... 気の緩み, 
% 41 ... 変異株シナリオ1(アメリカのgrowth rate), 42 ... 変異株シナリオ2(イギリスのgrowth rate)

for pindex = 7:7  %:length(PrefVector)
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
%     RetroH = 15;
    if isempty(find(SimDateEN == medical_start_date)) == 0
        medical_start = find(SimDateEN == medical_start_date);
    else
        medical_start = 1;
    end
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
    %Normal
    alpha_off = mean(alpha((dateEN >= datetime(2020,10,1)) & (datetime(2020,11,26)>= dateEN ))); % output loss without the state of emergency
    mean(alpha((dateEN >= datetime(2020,9,4)) & (datetime(2020,11,26)>= dateEN )))

%     alpha_off = mean(alpha((dateEN >= datetime(2020,9,4)) & (datetime(2020,11,26)>= dateEN ))); % output loss without the state of emergency
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
    delta_average = mean(delta_sample);
    
    pace = ps*3500000;
    pace_m = ps*900000; %520000?
    medical = ps*4700000*0.8;
    if medical_start == 1
        past_m = length(dateEN) - find(dateEN == medical_start_date) + 1;
        medical = medical - pace_m*past_m;
        if medical < 0
            medical = 0;
        end
    end        
    elderly = ps*36000000*0.8;
    ordinary = (125710000-36000000-3700000)*ps*0.8;
    elderly_total = ps*36000000;
    
    [V,deltaT,VT,real_pace] = vaccine_path(pace,pace_m,SimPeriod,medical,elderly,ordinary,medical_start,elderly_start,3,0.9,delta_average,elderly_total);
%     [V,deltaT,VT,real_pace] = vaccine_path(pace,pace_m,SimPeriod+past_m,medical,elderly,ordinary,medical_start,elderly_start,3,0.9,delta_average,elderly_total);
%     V(1:past_m) = 
%     initS - cumsum(V(1:past_m)) 
%     initR +  cumsum(V(1:past_m)) 
    
    
    betaT = beta_avg*ones(SimPeriod,1);
    betaT_temp = betaT_temp_ini;
    betaT(1,1) = (1+betaT_temp)*betaT(1,1);
    for i = 2:length(betaT)
        betaT_temp = betaT_temp * beta_rho;
        betaT(i) = betaT(i) * (1+betaT_temp);
    end
    if beta_option == 1
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
        % betaT(1:7) =  1.1*beta_avg;  %1.2183 ... 2020/11/26 - 2021/1/7
        betaT(1:9) =  1.08*beta_avg;  %1.1492 ... 2020/11/05 - 2020/12/31
        % betaT(1:3) =  1.3*beta_avg; % to achieve ERN = 1.418, % we have to adjust the inside of betaT(3:5) every week so that beta increases from March 22 for 3 weeks
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
        
        var_initial = var_initial_vector(pindex) / v_scale_vec(pindex);
        var_growth = var_growth_vec(1);
        logit_initial = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
%         var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept
        var_intercept = logit_initial;
        var_share = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
        beta_bar = beta_avg/(1+var_infection*var_initial);
        betaT = beta_bar*(1+var_infection*var_share);
        
        %Assuming that current beta is higher
        betaT_temp = betaT_temp_ini;
        betaT2(1,1) = (1+betaT_temp)*betaT(1,1);
        for i = 2:length(betaT)
            betaT_temp = betaT_temp * beta_rho;
            betaT2(i,1) = betaT(i) * (1+betaT_temp);
        end
        
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
            plot(xaxis_vec,plot_var_share(:,1),'-r','LineWidth',2)
            ax = gca;
            ax.YAxis.FontSize = 20;
            ax.XAxis.FontSize = 20;
            ax.YAxis.Exponent = 0;
            ytickformat('%,0.2f')
            xticks(find(WeekNumber_Sim==1))
            if l == 1
                legend('Variant Share','FontSize',10)
                title('Projected path of variants','FontSize',20,'FontWeight','normal')
                xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
                xlabel('Weeks','FontSize',20)
                ylabel('Share of Variants','FontSize',20)
            elseif l == 2
                legend('変異株割合','FontName',fn,'FontSize',10)
                title('新規感染者数の推移','FontSize',20,'FontWeight','normal')
                title('変異株シェアの推移','FontSize',20,'FontWeight','normal','FontName',fn)
                xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
                xlabel('週','FontSize',20,'FontName',fn)
                ylabel('変異株割合','FontSize',20,'FontName',fn)
            end
            xtickangle(45)
            xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
            xlim([0,length(SimDate)]);
            lgd = legend;
            lgd.Location = 'northwest';

            subplot(1,2,2)
            plot(xaxis_vec2(tt+1:end), plot_betaT(tt+1:end,1),'-b','LineWidth',2)
            hold on
            plot(xaxis_vec2(tt+1:end), plot_betaT2(tt+1:end,1),'-r','LineWidth',2.0)
            plot(xaxis_vec2,plot_betaT3,'-k','LineWidth',2)
            ax = gca;
            ax.YAxis.FontSize = 20;
            ax.XAxis.FontSize = 20;
            ax.YAxis.Exponent = 0;
            ytickformat('%,0.2f')
            xticks(find(WeekNumber_Sim2==1))
            xtickangle(45)
            if l == 1
                legend('Raw Infection Rate with Variants', 'Raw Infection Rate w/o Variants','Past Average','FontSize',10)
                title('Projected path of beta','FontSize',20,'FontWeight','normal')
                xticklabels(MonthWeekEN_Sim2(WeekNumber_Sim2==1))
                xlabel('Weeks','FontSize',20)
                ylabel('Infection Rate','FontSize',20,'FontName',fn)
            elseif l == 2
                legend('感染率1(変異株成長ケース)','感染率2(変異株成長ケース)','過去の平均','FontSize',10,'FontName',fn)
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
        saveas(figure(4101),[home 'Figures/' char(pref) '/beta_path' sprintf('%.0f', beta_option) '.png']);
        saveas(figure(4102),[home 'Figures/' char(pref) '/beta_path' sprintf('%.0f', beta_option) '_jp.png']);
        
        betaT = betaT2;


    elseif beta_option == 42
        var_initial = var_initial_vector(pindex) / v_scale_vec(pindex);
        var_growth = var_growth_vec(2);
        logit_initial = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
%         var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept
        var_intercept = logit_initial;
        var_share = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
        beta_bar = beta_avg/(1+var_infection*var_initial);
        betaT = beta_bar*(1+var_infection*var_share);
        %Assuming that current beta is higher
        betaT_temp = betaT_temp_ini;
        betaT2(1,1) = (1+betaT_temp)*betaT(1,1);
        for i = 2:length(betaT)
            betaT_temp = betaT_temp * beta_rho;
            betaT2(i,1) = betaT(i) * (1+betaT_temp);
        end

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
            plot(xaxis_vec,plot_var_share(:,1),'-r','LineWidth',2)
            ax = gca;
            ax.YAxis.FontSize = 20;
            ax.XAxis.FontSize = 20;
            ax.YAxis.Exponent = 0;
            ytickformat('%,0.2f')
            xticks(find(WeekNumber_Sim==1))
            if l == 1
                legend('Variant Share','FontSize',10)
                title('Projected path of variants','FontSize',20,'FontWeight','normal')
                xticklabels(MonthWeekEN_Sim(WeekNumber_Sim==1))
                xlabel('Weeks','FontSize',20)
                ylabel('Share of Variants','FontSize',20)
            elseif l == 2
                legend('変異株割合','FontName',fn,'FontSize',10)
                title('新規感染者数の推移','FontSize',20,'FontWeight','normal')
                title('変異株シェアの推移','FontSize',20,'FontWeight','normal','FontName',fn)
                xticklabels(MonthWeekJP_Sim(WeekNumber_Sim==1))
                xlabel('週','FontSize',20,'FontName',fn)
                ylabel('変異株割合','FontSize',20,'FontName',fn)
            end
            xtickangle(45)
            xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
            xlim([0,length(SimDate)]);
            lgd = legend;
            lgd.Location = 'northwest';

            subplot(1,2,2)
            plot(xaxis_vec2(tt+1:end), plot_betaT(tt+1:end,1),'-b','LineWidth',2)
            hold on
            plot(xaxis_vec2(tt+1:end), plot_betaT2(tt+1:end,1),'-r','LineWidth',2.0)
            plot(xaxis_vec2,plot_betaT3,'-k','LineWidth',2)
            ax = gca;
            ax.YAxis.FontSize = 20;
            ax.XAxis.FontSize = 20;
            ax.YAxis.Exponent = 0;
            ytickformat('%,0.2f')
            xticks(find(WeekNumber_Sim2==1))
            xtickangle(45)
            if l == 1
                legend('Raw Infection Rate with Variants', 'Raw Infection Rate w/o Variants','Past Average','FontSize',10)
                title('Projected path of beta','FontSize',20,'FontWeight','normal')
                xticklabels(MonthWeekEN_Sim2(WeekNumber_Sim2==1))
                xlabel('Weeks','FontSize',20)
                ylabel('Infection Rate','FontSize',20,'FontName',fn)
            elseif l == 2
                legend('感染率1(変異株成長ケース)','感染率2(変異株成長ケース)','過去の平均','FontSize',10,'FontName',fn)
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
        
        saveas(figure(4201),[home 'Figures/' char(pref) '/beta_path' sprintf('%.0f', beta_option) '.png']);
        saveas(figure(4202),[home 'Figures/' char(pref) '/beta_path' sprintf('%.0f', beta_option) '_jp.png']);

        
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
%     ERN_on = ERN_on_vector(pindex);
%     altA = (((ERN_on*(POP0/S(end))*((gammaT(1)+deltaT(1))/betaT(1))).^(1/k))-1)*(h(1)/h(2));
    th_off = th_off_vector(pindex)*7;
    ERN_on = ERN_on_vector(pindex);
    if beta_option == 41
        ERN_on = ERN_on_scenario1;
    elseif beta_option == 42
        ERN_on = ERN_on_scenario2;
    end
    ERN_now = ERN_now_vector(pindex);
    altA_now = (((ERN_now*(POP0/S(end))*((gammaT(1)+deltaT(1))/beta_avg)).^(1/k))-1)*(h(1)/h(2));
    altA_on = (((ERN_on*(POP0/S(end))*((gammaT(1)+deltaT(1))/beta_avg)).^(1/k))-1)*(h(1)/h(2));
    ERNCheck = (S(end)/POP0).*(((1+(h(2)/h(1))*alpha_off).^k).*beta_avg)./(gammaT(1)+deltaT(1));
    SimCases=3;
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
    TL = TL_vector{pindex};
    TL_index = TL_index_vector{pindex};
    for i = 1:length(TL)
%         if pindex <= 4
%             [DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData,NPath(:,i,1),SimERN(:,i,1)] = Covid_projection_control_gradual_alt_2(InitialValues,altA_now,altA_on,alpha_off,th_on,th_off,TL(i)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,4); %1 -> 4 gradual 3/15
%             % [DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData,NPath(:,i,1),SimERN(:,i,1)] = Covid_projection_control_gradual_alt(InitialValues,altA_now,altA_on,alpha_off,th_on,TL(i)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,4); %1 -> 4 gradual 3/15
%         elseif pindex > 4
            [DMat(1,i),AlphaMat(1,i),AlphaPath(:,i,1),SimData,NPath(:,i,1),SimERN(:,i,1)] = Covid_projection_control_gradual_off(InitialValues,altA_on,alpha_off,th_on,TL(i)*7,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,4,alpha(end)); %1 -> 4 gradual 3/15
%         end
    end
    %---------------------------------------------------------------%
    
    %---- 2. Different durations for economic recovery ---%
    DR = DR_vector{pindex};
    DR_index = [4,8]; %0->4,8->12 gradual 3/15
    if pindex == 3 || pindex == 4 || pindex == 7
        DR_index = [4,12];
    end 
    for i = 1:length(DR)
%         if pindex <= 4
%             [DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData,NPath(:,i,2),SimERN(:,i,2)] = Covid_projection_control_gradual_alt_2(InitialValues,altA_now,altA_on,alpha_off,th_on,th_off,TL_index(1)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,1+DR(i));
%             % [DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData,NPath(:,i,2),SimERN(:,i,2)] = Covid_projection_control_gradual_alt(InitialValues,altA_now,altA_on,alpha_off,th_on,TL_index(1)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,1+DR(i));
%         elseif pindex > 4
            [DMat(2,i),AlphaMat(2,i),AlphaPath(:,i,2),SimData,NPath(:,i,2),SimERN(:,i,2)] = Covid_projection_control_gradual_off(InitialValues,altA_on,alpha_off,th_on,TL_index(1)*7,betaT,gammaT,deltaT,V,h,k,POP0,hconstant,1+DR(i),alpha(end));
%         end
    end
    %-----------------------------------------------------%
    %
    %     %------------ 3. Different vaccine pace -------------%
    %     VP = [25,50,75,100:50:600];
    %     VP_index = [200,600,350];
    %     for i = 1:length(VP)
    %         [V,deltaT,VT,real_pace] = vaccine_path(ps*VP(i)*10000,ps*900000,SimPeriod,medical,elderly,ordinary,medical_start,elderly_start,3,0.9,delta_average,elderly_total);
    %         [DMat(3,i),AlphaMat(3,i),AlphaPath(:,i,3),SimData,NPath(:,i,3),SimERN(:,i,3)] = Covid_projection_control_gradual(InitialValues,altA,alpha_off,th_on,TL_index(2)*7,betaT,gammaT,deltaT,V,h_all,k,POP0,hconstant,1);
    %     end
    %     %-----------------------------------------------------%
    
    
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
    
    minAlpha = min(AlphaMat,[],'all');
    
    for y = 2:2 %1:2  %SimCases
        if y == 1
            %if (y==1 || y==4)
            TH = TL;
            TH_index = TL_index;
        elseif y == 2
            TH = DR;
            TH_index = DR_index;
        elseif y == 3
            TH = VP;
            TH_index = VP_index;
        end
        AlphaM = AlphaMat(y,:);
        AlphaM = AlphaM(~isnan(AlphaM));
        DM = DMat(y,:);
        DM = DM(~isnan(DM));
        AlphaM = (AlphaM - minAlpha)*prefGDP*10000;
        
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
        
        if y == 1
            ATL = AlphaM(waves==0);
            DTL = DM(waves==0);
            THTL = TH(waves==0);
            ATLblue = AlphaM(TH==TH_index(2));
            DTLblue = DM(TH==TH_index(2));
        elseif y == 2
            ADR = AlphaM(waves==0);
            DDR = DM(waves==0);
            THDR = TH(waves==0);
            ADRblue = AlphaM(TH==TH_index(2));
            DDRblue = DM(TH==TH_index(2));
        end
        
        for l = 1:2
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
                    if TH(i) == TH_index(1)
                        plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                        BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                        BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                        BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                    elseif TH(i) == TH_index(2)
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
                if l == 1
                    for i = 1:length(TH)
                        if TH(i) == TH_index(1)
                            plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                            BackDataN(:,1,y) = [N(end-7:end);NPath(:,i,y)];
                            BackDataAlpha(:,1,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                            BackDataERN(:,1,y) = [ERN(end-7:end);SimERN(:,i,y)];
                        elseif TH(i) == TH_index(2)
                            plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                            BackDataN(:,2,y) = [N(end-7:end);NPath(:,i,y)];
                            BackDataAlpha(:,2,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                            BackDataERN(:,2,y) = [ERN(end-7:end);SimERN(:,i,y)];
                        elseif TH(i) == TH_index(3)
                            plot([N;NPath(:,i,y)]/7,'-k','LineWidth',2,'DisplayName',sprintf('%.0f',TH(i)))
                            BackDataN(:,3,y) = [N(end-7:end);NPath(:,i,y)];
                            BackDataAlpha(:,3,y) = [alpha(end-7:end);AlphaPath(:,i,y)];
                            BackDataERN(:,3,y) = [ERN(end-7:end);SimERN(:,i,y)];
                        else
                            plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',sprintf('%.0f',TH(i)))
                        end
                        hold on
                    end
                elseif l == 2
                    for i = 1:length(TH)
                        if TH(i) == TH_index(1)
                            plot([N;NPath(:,i,y)]/7,'-r','LineWidth',2,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                        elseif TH(i) == TH_index(2)
                            plot([N;NPath(:,i,y)]/7,'-b','LineWidth',2,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                        elseif TH(i) == TH_index(3)
                            plot([N;NPath(:,i,y)]/7,'-k','LineWidth',2,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                        else
                            plot([N;NPath(:,i,y)]/7,'--','LineWidth',0.3,'DisplayName',[sprintf('%.0f',TH(i)),'万'])
                        end
                        hold on
                    end
                end
            end
            plot(N/7,'k','LineWidth',2,'HandleVisibility','off')
            xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
%             if pindex <= 4
%                 xline(Tdata+find(SimDateEN == datetime(2021,3,18)),'--','LineWidth',1.5,'HandleVisibility','off');
%             end
            if y == 1
                xline(Tdata+find(SimDateEN == datetime(2021,4,1)),'--','LineWidth',1.5,'HandleVisibility','off');
            end
            ax = gca;
            ax.YAxis.FontSize = 20;
            ax.XAxis.FontSize = 20;
            ax.YAxis.Exponent = 0;
            ytickformat('%,6.0f')
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
            lgd = legend;
            lgd.NumColumns = 2;
            if y == 3
                lgd.Location = 'north';
            end
            xtickangle(45)
            if beta_option == 41 || beta_option == 42
                xlim([Tdata-7 Tdata+52])
            else 
                xlim([Tdata-7 Tdata+28])
            end 
            if y == 3
                xlim([Tdata-7 Tdata+28])
            end
            
            %--- Number of cumulative deaths ---%
            subplot(1,2,2)
            plot(AlphaM(waves==0),DM(waves==0),'-go','LineWidth',2,'MarkerSize',10);
            hold on
            plot(AlphaM(waves==1),DM(waves==1),'-mo','LineWidth',2,'MarkerSize',10);
            hold on
            plot(AlphaM(waves==2),DM(waves==2),'-bo','LineWidth',2,'MarkerSize',10);
            hold on
            if y ~= 3
                text(AlphaM,DM,string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
                hold on
            end
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
            if y == 3
                hold on
                scatter(AlphaM(TH==TH_index(3)),DM(TH==TH_index(3)),250,'black','filled');
            end
            
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
        end
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
                writetable(TN,[home 'Figures/' char(pref) '/BackData_Thresholds' char(pref) '_' sprintf('%.0f', beta_option) '_v.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                writetable(TAD,[home 'Figures/' char(pref) '/BackData_Thresholds' char(pref) '_' sprintf('%.0f', beta_option) '_v.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
            elseif y == 2
                writetable(TN,[home 'Figures/' char(pref) '/BackData_GradualRecovery' char(pref) '_' sprintf('%.0f', beta_option) '_v.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                writetable(TAD,[home 'Figures/' char(pref) '/BackData_GradualRecovery' char(pref) '_' sprintf('%.0f', beta_option) '_v.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
            elseif y == 3
                writetable(TN,[home 'Figures/' char(pref) '/BackData_VaccineVariation' char(pref) '_' sprintf('%.0f', beta_option) '_v.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
                writetable(TAD,[home 'Figures/' char(pref) '/BackData_VaccineVariation' char(pref) '_' sprintf('%.0f', beta_option) '_v.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
            end
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
    
    %--- Save figures ---%
    if figure_save == 1
        
        saveas(figure(2),[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
        if beta_option == 41
            saveas(figure(41012),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '.png']);
            saveas(figure(41121),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_jp.png']);
            saveas(figure(41013),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '.png']);
            saveas(figure(41131),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_jp.png']);
        elseif beta_option == 42
            saveas(figure(42012),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '.png']);
            saveas(figure(42121),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_jp.png']);
            saveas(figure(42013),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '.png']);
            saveas(figure(42131),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_jp.png']);
        else
            saveas(figure(12),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_v.png']);
            saveas(figure(121),[home 'Figures/' char(pref) '/Thresholds' sprintf('%.0f', beta_option) '_v_jp.png']);
            saveas(figure(13),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_v.png']);
            saveas(figure(131),[home 'Figures/' char(pref) '/GradualRecovery' sprintf('%.0f', beta_option) '_v_jp.png']);
        %         saveas(figure(14),[home 'Figures/' char(pref) '/VaccineVariation' sprintf('%.0f', beta_option) '.png']);
        %         saveas(figure(141),[home 'Figures/' char(pref) '/VaccineVariation' sprintf('%.0f', beta_option) '_jp.png']);
        %         saveas(figure(200),[home 'Figures/' char(pref) '/ThresholdsAndGradual' sprintf('%.0f', beta_option) '_jp.png']);
        end
    end
    
    
    
    %--- Update forecast error data ---%
    % if iPC==1
    %     FE = load('\Users\shcor\Dropbox\fujii_nakata\Website\Codes\Forecast_Tokyo.mat');
    % else
    %     FE = load('/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Forecast_Tokyo.mat');
    % end
    %
    % dDActual = [FE.dDActual;dD(end)];
    % dNActual = [FE.dNActual;N(end)];
    %
    % if iPC==1
    %     save('\Users\shcor\Dropbox\fujii_nakata\Website\Codes\Forecast_Tokyo.mat','dDActual','dNActual','-append');
    % else
    %     save('/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Forecast_Tokyo.mat','dDActual','dNActual','-append');
    % end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Forecast error analysis (using next week's data) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if iPC==1
    %   NextD = load('\Users\shcor\Dropbox\fujii_nakata\Website\Codes\Tokyo20210221.mat');
    % else
    %   NextD = load('/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/Tokyo20210221.mat');
    % end
    
    % AlphaNext = NextD.alpha(end);
    % [CumDNext,AverageAlphaNext,SimDataNext,SimNNext]=Covid_projection(InitialValues,AlphaNext,betaT,gammaT,deltaT,V,h,k,POP0,hconstant);
    % dDForecast = [FE.dDForecast;(CumDNext-D(end))];
    % dNForecast = [FE.dNForecast;SimNNext];
    % if dNForecast(end) ~= dNForecast(end-1) && length(dNForecast) == NextD.Tdata
    %   if iPC==1
    %       save('\Users\shcor\Dropbox\fujii_nakata\Website\Codes\Forecast_Tokyo.mat','dDForecast','dNForecast','-append');
    %   else
    %       save('/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes/Forecast_Tokyo.mat','dDForecast','dNForecast','-append');
    %   end
    % end
    %
    %     if mat_save==1
    %         save([pref char(datetime('today','Format','yyyMMdd')) '.mat'])
    %     end
    
end %end of prefecture loop


% figure(9991)
% plot(M)
% hold on
% plot(Malt)
% legend('M','Malt')
% 
% figure(9992)
% plot(beta_sample)
% legend('beta sample')

figure(9993)
plot(beta(end-17:end))
legend('beta path in the past 17 weeks')

figure(9994)
plot(betaT(1:10))
legend('Simulated Beta')

cumV = zeros(length(V),1);
for i = 1:length(V)
    cumV(i,1) = sum(V(1:i));
end

% figure(9995)
% plot(SimData(:,1))
% hold on 
% plot(SimData(:,3))
% plot(cumV)
% legend('S','R','V')
% 
% 
% figure(9996)
% plot(SimData(:,2))
% hold on 
% plot(SimData(:,4))
% legend('I','D')

