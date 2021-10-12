% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
iPC=0; % 0 for Mac, 1 for Windows
if iPC==1
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\Cabinet_2021APR28_DR\';
    %     home = 'C:\Users\tak10\Dropbox\fujii_nakata\Website\Codes\Cabinet_2021APR28_DR\';
else
    %home ='/Users/machikohei/Dropbox/fujii_nakata/Website/Codes/Cabinet_2021APR28_DR/';
    home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/Cabinet_2021APR28_DR/';
    %home = '/Users/ymaeda/Documents/ym_doc/Tokyo_Univrsity_MA/Reserach Assistant/fujii_nakata/Codes/';
    %home = '/Users/shotaro/Dropbox/fujii_nakata/Website/Codes/Cabinet_2021APR28_DR/';
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 0;      % save back data in the "Figure" folder
vaccine_figure_loop = 0; % =0 appear only once; =1 appear every loop; 
beta_figure_loop = 0; % =0 appear only once; =1 appear every loop; 
%======================================================================%

%================== Parameter Values ============================%
PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka','Aichi','Fukuoka','Hyogo'};
GDPVector = [106,36,23,21,40,41,20,20]; % 兆円 per year, one trillion yen (chou-yen)

parameter % m file for setting parameters

for pindex = 1:1 %:length(PrefVector) %change this parameter for prefecture
    %---------------------------------------------------%
    %--------- Import Prefecture-specifc data ----------%
    %---------------------------------------------------%
    % Covid data are recorded at weekly frequencies (Mon-Sun)
    % The first week start on January 20 (Mon), 2020
    import_prefecture
    
    %--- Construct weekly vaccine data ---%
    [V1_w, V2_w] = import_vaccine(home, iPC,dateEN,ps);
    
    %--- Constructing the reference level of output ---%
    [potentialGDP, referenceGDP, alpha] = construct_GDP(GDP,TdataGDP);
    
    %--- Regress mobility on alpha to estimate the elasticity h ---%　関数化？
    [Malt,h_all,h_all_se,h,h_se] = estimate_h(M,alpha,TdataGDP,RetroH,hconstant);
    
    %--- Plot mobility data ---%
    figure(2)
    plot_mobility(Malt,alpha,Tdata,TdataGDP,MonthWeekJP,xtick1,fs,16)
    if figure_save == 1
        saveas(figure(2),[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Main analysis starts here %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% CHANGE HERE %%%%%%%%%%%%%%%%%%%%%%%%%
    TH = 1.0:0.1:1.7; %0.3:0.05:0.7;
    TH_index = [1.0,1.3,1.5];
    % Change both paces and infection rate at the same time
    %     TH = [3600000,7000000,3600000,7000000];
    %     TH_index = [3600000,7000000,3600000,7000000];
    %     TH2 = [0.3, 0.3 ,0.5, 0.5];
    %     TH2_index = [0.3, 0.5];
    if max(TH)<3
            ft = '%.2f';
    else
            ft = '%.0f';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Construct Empty Matrices
    DMat = nan(1,length(TH));
    AlphaMat = nan(1,length(TH));
    SimData = nan(SimPeriod+1,4,length(TH));
    AlphaPath = nan(SimPeriod,length(TH));
    NPath = nan(SimPeriod,length(TH));
    SimERN = nan(SimPeriod,length(TH));
    BackDataN = zeros(SimPeriod+8,length(TH_index));
    BackDataAlpha = zeros(SimPeriod+8,length(TH_index));
    BackDataERN = zeros(SimPeriod+8,length(TH_index));
    BackDataDA = zeros(length(TH),3);
    
    for iTH = 1:length(TH)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% CHANGE HERE %%%%%%%%%%%%%%%%%%%%%%%%%
        var_infection = TH(iTH) - 1
        %         paces_ori = TH(iTH);
        %         var_infection = TH2(iTH);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %--- Compute the history of S, I, R, D in the data period ---%
        S = zeros(Tdata+1,1);
        I = zeros(Tdata+1,1);
        R = zeros(Tdata+1,1);
        D = zeros(Tdata+1,1);
        S(1)=POP0;
        for i = 1:Tdata
            S(i+1)=S(i)-N(i)-E1*V1_w(i)-(E2-E1)*V2_w(i);
            I(i+1)=I(i)+N(i)-gamma*I(i)-dD(i);
            R(i+1)=R(i)+gamma*I(i)+E1*V1_w(i)+(E2-E1)*V2_w(i);
            D(i+1)=D(i)+dD(i);
            if i > TdataGDP
                GDP(i) = referenceGDP(i)*(1-alpha(i));
            end
        end
        
        %--- Compute the history of time-varying parameters ---%
        delta = (D(2:Tdata+1)-D(1:Tdata))./I(1:Tdata);                              % death rate
        beta_tilde = (POP0.*N(1:Tdata))./((S(1:Tdata).*I(1:Tdata)));   % overall infection rate
        ERN = (S(1:end-1)/POP0).*beta_tilde./(gamma+delta);                                        % effective reproduction number
        if hconstant == 0
            beta = beta_tilde./(1+h_all*alpha).^k;                                      % raw infection rate
        elseif hconstant == 1
            %     beta = beta_tilde./(h(1)+h(2)*alpha).^k;
            beta = beta_tilde./(1+(h_all(2)/h_all(1))*alpha).^k;
        end
        
        %--- Construct time series of parameters ---%
        gammaT = gamma*ones(SimPeriod,1);
        delta_sample = delta(end-RetroPeriod+1:end);
        delta_average = sum(delta_sample.*(I(end-RetroPeriod+1:end)/sum(I(end-RetroPeriod+1:end))));
        
        %--- Construct vaccine dstribution ---%
        paces = ps*paces_ori; %3600000;
        vacpath = zeros(SimPeriod,1);
        vacpath(1+sw_vacpath:gradual_paces) = (paces/(gradual_paces-sw_vacpath)):(paces/(gradual_paces-sw_vacpath)):paces;
        vacpath(gradual_paces+1:end) = paces*ones(SimPeriod-gradual_paces,1);
        elderly_total = ps*elderly_jp;
        medical_total = ps*medical_jp;
        ordinary_total = ps*ordinary_jp;
        medical = medical_total*accept_share;
        elderly = elderly_total*accept_share;
        ordinary = ordinary_total*accept_share;
        % Construct Vaccine Path
        [V,deltaT,VT] = vaccine_distribution_medical(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2,ps,POP0,3);
        %[V,deltaT,VT] = vaccine_distribution_simple(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2);
        delta_ss = delta_average*(0.09/1.28);
        delta_ss = (delta_average - delta_ss)*(1-accept_share*D2)+delta_ss; %death rate after full vaccination of elderly = 2021/09/30
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        beta_r = 0;
        for retrop = retro_lb:retro_ub
            beta_r = beta_r + mean(beta(end-retrop+1:end));
        end
        beta_avg = beta_r/(retro_ub-retro_lb+1);
        betaT = beta_avg*ones(SimPeriod,1);
        
        var_intercept = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
        % var_intercept = logit_initial - var_growth * (length(dateEN)-FirstObsTime); %Counterfactual intercept
        
        % step 1: extrapolation in the variant share in the last 17 weeks     -mean(retro_lb,retro_ub)+1 = -17 + 1
        var_share = exp((1:SimPeriod)'*var_growth+var_intercept)./(1+exp((1:SimPeriod)'*var_growth+var_intercept));
        var_share_prev = exp((-mean(retro_lb,retro_ub)+1:0)'*var_growth+var_intercept)./(1+exp((-mean(retro_lb,retro_ub)+1:0)'*var_growth+var_intercept));
        % step 2: relative infection rate in the last 17 weeks
        relative_var_infection_prev = 1 + var_share_prev * var_infection;
        % step 3:
        no_var_beta = beta(end-mean(retro_lb,retro_ub)+1:end) ./ relative_var_infection_prev;
        beta_bar = mean(no_var_beta);
        betaT = beta_bar*(1+var_infection*var_share);
        betaT_woAR1 = betaT;
        
        betaT_temp = betaT_temp_ini;
        betaT(1,1) = (1+betaT_temp)*betaT(1,1);
        for i = 2:length(betaT)
            betaT_temp = betaT_temp * beta_rho;
            betaT(i) = betaT(i) * (1+betaT_temp);
        end
        
        alpha_off = mean(alpha((dateEN >= datetime(2020,10,1)) & (datetime(2020,11,26)>= dateEN ))); % output loss without the state of emergency
        InitialValues = [S(end),I(end),R(end),D(end)];
        altA_on = (((ERN_on*(POP0/S(end))*((gammaT(1)+delta_average)/beta_avg)).^(1/k))-1)*(h(1)/h(2));
        ERNCheck = (S(end)/POP0).*(((1+(h(2)/h(1))*alpha_off).^k).*beta_avg)./(gammaT(1)+delta_average);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if pindex < 5
            [DMat(iTH),AlphaMat(iTH),AlphaPath(:,iTH),SimData(:,:,iTH),NPath(:,iTH),SimERN(:,iTH)] ...
                = Covid_projection_control_gradual_off_threshold_delta2...
                (InitialValues,altA_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,betaT,gammaT,deltaT,delta_ss,V,h,k,POP0,hconstant,DRi,alpha(end));
%             [~,minAlphaMat(iTH)] ...
%                 = Covid_projection_control_gradual_off_threshold_delta2...
%                 (InitialValues,altA_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,beta_avg*ones(SimPeriod,1),gammaT,deltaT,delta_ss,V,h,k,POP0,hconstant,DRi,alpha(end));
        else
            [DMat(iTH),AlphaMat(iTH),AlphaPath(:,iTH),SimData(:,:,iTH),NPath(:,iTH),SimERN(:,iTH)] ...
                = Covid_projection_control_gradual_threshold_delta2...
                (InitialValues,altA_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,betaT,gammaT,deltaT,delta_ss,V,h,k,POP0,hconstant,DRi);
%             [~,minAlphaMat(iTH)] ...
%                 = Covid_projection_control_gradual_threshold_delta2...
%                 (InitialValues,altA_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,beta_avg*ones(SimPeriod,1),gammaT,deltaT,delta_ss,V,h,k,POP0,hconstant,DRi);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Construct Backdata at Index Values
        if abs(TH(iTH) - TH_index(th_wave)) < 0.0001
            BackDataN(:,th_wave) = [N(end-7:end);NPath(:,iTH)];
            BackDataAlpha(:,th_wave) = [alpha(end-7:end);AlphaPath(:,iTH)];
            BackDataERN(:,th_wave) = [ERN(end-7:end);SimERN(:,iTH)];
            th_wave = th_wave + 1;
            if th_wave > length(TH_index)
                th_wave = length(TH_index);
            end
        end
        
        % Plot betaT 
        if beta_figure_loop == 0
            if iTH == 1
                figure(400) % Figure for BetaT with Variant Share
                set(gcf,'Position',[100,100,1200,500])
                plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaT,betaT_woAR1,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
                if figure_save == 1
                    saveas(figure(400),[home 'Figures/' char(pref) '/beta_path' '.png']);
                end
                
                betaT_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT;
                betaT_woAR1_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT_woAR1;
                beta_tilde_avg = beta_avg*((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k);
                
                figure(401)
                set(gcf,'Position',[100,100,1200,500])
                plot_beta(var_initial,var_share,beta_tilde,beta_tilde_avg,betaT_tilde,betaT_woAR1_tilde,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
                title('β tildeの推移（',sprintf(ft,TH(iTH)),'）','FontSize',20,'FontWeight','normal','FontName',fn)
            end
        elseif beta_figure_loop == 1
            figure(400+iTH) % Figure for BetaT with Variant Share
            set(gcf,'Position',[100,100,1200,500])
            plot_beta(var_initial,var_share,beta,beta_avg,betaT,betaT_woAR1,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
            title(string(['βの推移（',sprintf(ft,TH(iTH)),'）']),'FontSize',20,'FontWeight','normal','FontName',fn)
            
            betaT_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT;
            betaT_woAR1_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT_woAR1;
            beta_tilde_avg = beta_avg*((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k);
            
            figure(440 + iTH)
            set(gcf,'Position',[100,100,1200,500])
            plot_beta(var_initial,var_share,beta_tilde,beta_tilde_avg,betaT_tilde,betaT_woAR1_tilde,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
            title(string(['β tildeの推移（',sprintf(ft,TH(iTH)),'）']),'FontSize',20,'FontWeight','normal','FontName',fn)
        end
        
        %figure for vaccine path
        if vaccine_figure_loop == 0
            if iTH == 1
                plot_vaccinepath(200,VT,V1_w,V2_w,SimPeriod,ps,MonthWeekJP,WeekNumber,Tdata,fs,fn);
                plot_deltapath(201,delta,deltaT,delta_average,MonthWeekJP,WeekNumber,Tdata,fs,fn);
            end
        elseif vaccine_figure_loop == 1
            plot_vaccinepath(200+iTH,VT,V1_w,V2_w,SimPeriod,ps,MonthWeekJP,WeekNumber,Tdata,fs,fn);
            subplot(1,2,1)
            title(string(['新規ワクチン接種本数（週ごと）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
            subplot(1,2,2)
            title(string(['累計ワクチン接種本数（週ごと）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
            plot_deltapath(240+iTH,delta,deltaT,delta_average,MonthWeekJP,WeekNumber,Tdata,fs,fn);
            title(string(['致死率（現在のレベルで標準化）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
        end
        
    end % End of iTH loop
        
    %minAlpha = min(minAlphaMat); %関数化, minimum alpha when variants have no effects.
    minAlpha = alpha_off;
    
    AlphaM = AlphaMat(~isnan(AlphaMat));
    AlphaM = (AlphaM - minAlpha)*prefGDP*10000;
    DM = DMat(~isnan(DMat));
    BackDataDA(1:length(TH),:) = [round(AlphaM'),round(DM'),TH'];
    
    %--- Record how many times on and off are triggered ---%
    waves = zeros(1,length(TH));
    for i = 1:length(TH)
        svec = zeros(SimPeriod-1,1);
        for t = 1:SimPeriod-1
            svec(t) = AlphaPath(t+1,i)-AlphaPath(t,i);
        end
        waves(i) = sum(svec>0);
    end
    %waves_th(1:length(waves),y,th_off_index) = waves;
    
    for l = 1:2 %1:2 when english version needed
        % Generate graphs for the website
        lng = language{l};
        figname = 100 + l;
        figure(figname)
        set(gcf,'Position',[100,100,1200,500])
        subplot(1,2,1)
        plot_SimN2(TH,TH_index,N,NPath,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,fs,fn,ft,l);
        
        %--- Number of cumulative deaths ---%
        subplot(1,2,2)
        plot_Tradeoff(AlphaM,DM,waves,TH,TH_index,l,linecolor,fs,fn)
        if figure_save == 1
            saveas(figure(figname),[home 'Figures/' char(pref) '/MainResult_' char(lng) '.png']);
        end
    end %End of language loop = figure loop
    
    if data_save == 1
        titleN = strings(1,1+length(TH_index)*3);
        titleN(1) = "週";
        for ti = 1:length(TH_index)
            titleN(1,1+ti) = string(['新規感染者数（',sprintf(ft,TH_index(ti)),'）']);
            titleN(1,1+length(TH_index)+ti) = string(['経済活動（',sprintf(ft,TH_index(ti)),'）']);
            titleN(1,1+length(TH_index)*2+ti) = string(['実効再生産数（',sprintf(ft,TH_index(ti)),'）']);
        end
        TN = table([titleN;MonthWeekJP(Tdata-7:end-1),round(BackDataN(:,1:length(TH_index))/7),round(100*(1-BackDataAlpha(:,1:length(TH_index))),1),round(BackDataERN(:,1:length(TH_index)),2)]);
        titleAD = ["経済損失（億円）","死亡者数","ケース"];
        TAD = table([titleAD;BackDataDA(1:length(TH),:)]);
        
        writetable(TN,[home 'Figures/' char(pref) '/BackData_' char(pref)  '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
        writetable(TAD,[home 'Figures/' char(pref) '/BackData_' char(pref) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
    end
    
    figname = 140 + pindex;
    figure(figname)
    set(gcf,'Position',[100,100,1400,600])
    subplot(1,3,1)
    plot_SimN2(TH,TH_index,N,NPath,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,fs,fn,ft,l)
    subplot(1,3,2)
    plot_Tradeoff(AlphaM,DM,waves,TH,TH_index,2,linecolor,fs,fn)
    subplot(1,3,3)
    plot_Alpha(alpha,AlphaPath,TH,TH_index,MonthWeekEN,WeekNumber,Tdata,linecolor,ft,fs,fn,2)
    
end %end of prefecture loop

