% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata
%
% To investigate the effects of Indian Variant, make "variant_India" = 1.
%


clear variables
close all
iPC=0; % 0 for Mac, 1 for Windows
if iPC==1
     home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
    home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/';
    home = '/Users/koheimachi/Dropbox/fujii_nakata/Website/Codes/';
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 0;      % save back data
vaccine_figure_loop = 0; % =0 appear only once; =1 appear every loop;
beta_figure_loop = 0; % =0 appear only once; =1 appear every loop;
vaccine_disp_switch = 1; % =0 not display the summary of # of the vaccinated ; =1 display
variant_India = 1; % = 0 ... without Indian variant, = 1 ... with Indian variant
ICU_nation = 1; % = 1 use national definition (NHK data), = 0 use data from Tokyo Keizai
% in the "Figure" folder
fs = 20;            % common font size for many figures
ldfs = 8;           % legend font size for vaccine path
ft = '%.0f';
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
linecolor = {'black','blue','red'};
language = {'EN','JP'};
%======================================================================%

%================== Model Fixed Parameter Values ============================%
parameter

%================== Parameter Values (Prefecture Specific) ============================%
prefecture_parameter


for pindex = [1] %:length(PrefVector) %change this parameter for prefecture
    %--- Import data ---%
    % Covid data are recorded at weekly frequencies (Mon-Sun)
    % The first week start on January 20 (Mon), 2020
    import_prefecture

    %--- Construct weekly vaccine data ---%　
    [V1_medical_ori, V1_medical, V2_medical_ori, V2_medical,...
        V1_elderly_ori, V1_elderly, V2_elderly_ori, V2_elderly, ...
        vs_MT,vs_ET,vs_M1,vs_E1,vs_M2,vs_E2] ...
        = ImportVaccineData(home,iPC,Data,pref,dateEN,ps,vaccine_disp_switch);

     %--- Constructing the reference level of output ---%
    [potentialGDP, referenceGDP, alpha] = construct_GDP(GDP,TdataGDP);

    %--- Regress mobility on alpha to estimate the elasticity h ---%
    [Malt,h_all,h_all_se,h_ori,h_se] = estimate_h(M,alpha,TdataGDP,RetroH,hconstant);


    %--- Plot mobility data ---%
    figname = string(['Mobility_GDP_' char(pref)]);
    f = figure('Name',figname);
    plot_mobility(Malt,alpha,Tdata,TdataGDP,MonthWeekJP,xtick1,fs,16)
    if figure_save == 1
        %saveas(figure(2),[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
        saveas(f,[home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
    end

    %--- Import ICU data (Option 2) ---%
    ICU = zeros(Tdata+1,1);
    BED = zeros(Tdata,1);
    if ICU_nation == 1
        ICU(2:Tdata+1,1) = Data(:,22);
        BED(1:Tdata,1) = Data(:,23);
    else
        ICU(2:Tdata+1,1) = Data(:,21);
    end
    %--- Plot ICU data ---%
    figname = string(['ICU_transition_' char(pref)]);
    f = figure('Name',figname);
    plot(ICU, 'LineWidth', 1.5)
    title('Transition of  ICU')
    ytickformat('%,6.0f')
    xticks(find(WeekNumber==1))
    xticklabels(MonthWeekJP(WeekNumber==1))
    lgd.NumColumns = 2;
    xtickangle(45)


    %--- Compute the history of S, I, R, D in the data period ---%
    [S,I,R,D]...
        = SIRD(Tdata,POP0,N,E1,E2,...
                V1_elderly,V1_medical,V2_elderly,V2_medical,...
                gamma,dD,TdataGDP,referenceGDP,alpha);

    %--- Compute the history of time-varying parameters ---%
    [delta,beta_tilde,ERN,beta,ICU_inflow,...
        gammaT,delta_average,ICU_inflow_avg]...
            = Time_Series_Average(S,I,D,ICU,dD,N,Tdata,SimPeriod,...
                RetroPeriod,POP0,gamma,hconstant,h_all,alpha,k,...
                gamma_ICU,ICU_adjustment);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%% Main analysis starts here %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% CHANGE HERE %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     TH = {100:100:800,70:10:100,80:10:130,50:10:100,[100:100:1000, 1060],20:10:100,10:10:100,50:10:100}; % 解除基準分析をコントロールしている Cell array
    % Different threshold for lifting the state of emergency
    fig_title = "Baseline";
    data_title = "Baseline";
    if pindex == 1 || pindex == 5
        state = 1;
    else
        state = 0;
    end

    alpha_May = mean(alpha((dateEN >= datetime(2020,5,07)) & (datetime(2020,5,28)>= dateEN )));
    alpha_Jan = mean(alpha((dateEN >= datetime(2021,1,07)) & (datetime(2021,1,28)>= dateEN )));
    if pindex == 1
        alpha_scale = 0.94;
        alpha_on = alpha_scale*(0.5*alpha_May+0.5*alpha_Jan);
        TH = [200:50:500];
        TH_index = [300,400,500];
        h_scale = 1;
        beta_shock_after_emergency = 0.15;
        rho_after_emergency = 0.8;
    elseif pindex == 5
        alpha_on = 1*alpha_May;
        TH = [50:50:200];
        TH_index = [50,100,200];
        h_scale = 1;
        beta_shock_after_emergency = 0.7;
        rho_after_emergency = 0.85;
        linecolor = {'black','red','blue','green'};
    end

    h = h_ori;

    %h(2) = h_scale*h_ori(2);

%     if variant_India == 1
% %         DRi = 10;
%         TH_index = [300,500,700];
%     end

    if max(TH)<3
        ft = '%.2f';
    else
        ft = '%.0f';
    end


% for iAlpha = 1:length(alpha_on_vector_sim)
%     alpha_on = alpha_on_vector_sim(iAlpha)


    DMat = nan(1,length(TH));
    AlphaMat = nan(1,length(TH));
    SimData = nan(SimPeriod+1,5,length(TH));
    AlphaPath = nan(SimPeriod,length(TH));
    NPath = nan(SimPeriod,length(TH));
    SimERN = nan(SimPeriod,length(TH));
    THonPath = nan(SimPeriod,length(TH));
    BackDataN = zeros(SimPeriod+8,length(TH_index));
    BackDataAlpha = zeros(SimPeriod+8,length(TH_index));
    BackDataERN = zeros(SimPeriod+8,length(TH_index));
    BackDataDA = zeros(length(TH),3);

    for iTH = 1:length(TH)

        %%%%%%%%%%% Change here %%%%%%%%%%%%%%%%%%%
        th_off1 = TH(iTH)*7
        th_off2 = TH(iTH)*7;
        th_off3 = TH(iTH)*7;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %--- Construct vaccine distribution and delta path---%
        [V,deltaT,VT] = ...
            vaccine_distribution_medical(V1_medical,V2_medical,...
            V1_elderly,V2_elderly,...
            elderly_jp,medical_jp,ordinary_jp,accept_share,...
            delta_average,lag,medical_duration,...
            paces_ori,sw_vacpath,gradual_paces,...
            E1,E2,D1,D2,ps,POP0,SimPeriod);

        [var_share,var_prev,var_initial] = var_share_prev(Data,SimPeriod,var_ss,var_growth);
        deltaT = construct_delta_variant(RetroPeriod,...
                    retro_lb,retro_ub,deltaT,delta,delta_average,var_prev,var_share,var_infection_delta,I);


        %% figure for vaccine path
        if vaccine_figure_loop == 0
            if iTH == 1
                plot_vaccinepath(200,VT,V1_medical,V2_medical,V1_elderly,V2_elderly,SimPeriod,ps,MonthWeekJP,WeekNumber,Tdata,fs,ldfs,fn);
                plot_deltapath(201,delta,deltaT,deltaT(1),MonthWeekJP,WeekNumber,Tdata,fs,fn,iTH);
            end
        elseif vaccine_figure_loop == 1
            plot_vaccinepath(200+iTH,VT,V1_medical,V2_medical,V1_elderly,V2_elderly,SimPeriod,ps,MonthWeekJP,WeekNumber,Tdata,fs,ldfs,fn);
            subplot(1,2,1)
            title(string(['新規ワクチン接種本数（週ごと）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
            subplot(1,2,2)
            title(string(['累計ワクチン接種本数（週ごと）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
            plot_deltapath(240,delta,deltaT,deltaT(1),MonthWeekJP,WeekNumber,Tdata,fs,fn,iTH);
            title(string(['致死率（現在のレベルで標準化）']), 'FontSize',fs,'FontName',fn);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Beta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        beta_r = 0;
        for retrop = retro_lb:retro_ub
            beta_r = beta_r + mean(beta(end-retrop+1:end));
        end
        beta_avg = beta_r/(retro_ub-retro_lb+1);

        [betaT,betaT_woAR1,beta_bar] = construct_beta(SimPeriod,retro_lb,retro_ub,...
            beta,beta_avg,var_prev,var_share,var_infection);

        %%%%  インド株　　　%%%%
        %%%%%%%%%%%%%%%%%%%%%​
        if variant_India == 1
            var_intercept2 = log(var_initial2/(1-var_initial2));
            var_share2 = exp((1:SimPeriod)'*var_growth2+var_intercept2)./(1+exp((1:SimPeriod)'*var_growth2+var_intercept2));
            var_power = 1+var_share2*var_growth2;
            var_boost = ones(SimPeriod,1);
            var_boost(var_start:end) = var_power(1:SimPeriod-var_start+1);
            var_shareInd = ones(SimPeriod,1)*var_initial2;
            var_shareInd(var_start:end) = var_share2(1:SimPeriod-var_start+1);
            betaT_Eng = betaT;
            betaT_Eng = beta_AR1(betaT_temp_ini, beta_rho, betaT_Eng, start_beta);
            betaT = betaT.*var_boost;
        end


        betaT = beta_AR1(betaT_temp_ini, beta_rho, betaT, start_beta);

        alpha_off = mean(alpha((dateEN >= datetime(2020,2,7)) & (datetime(2020,2,28)>= dateEN ))); % output loss without the state of emergency
        InitialValues = [S(end),I(end),R(end),D(end),ICU(end)];


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [DMat(iTH),AlphaMat(iTH),AlphaPath(:,iTH),SimData(:,:,iTH),NPath(:,iTH),SimERN(:,iTH),THonPath(:,iTH),SimICU(:,iTH),betaShock] ...
                = Covid_projection_ICU4...
                (InitialValues,alpha_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,betaT,gammaT,deltaT,deltaT(1),V,h,k,POP0,hconstant,DRi,state,ICU_inflow_avg,gamma_ICU,beta_shock_after_emergency,rho_after_emergency);

        % Plot betaT
        if beta_figure_loop == 0
            if iTH == 1
                %figure(400) % Figure for BetaT with Variant Share
                f=figure('Name','BetaPath');
                set(gcf,'Position',[100,100,1200,500])
                if variant_India == 1
                    %plot_beta_India(var_initial,var_initial2,var_share,var_shareInd,beta,beta_avg,betaT,betaT_Eng,betaT_woAR1,dateD,SimDate,Tdata,MonthWeekJP,MonthNumber,WeekNumber,fn)
                    plot_beta_India( var_initial,var_initial2,var_share,var_shareInd,beta,beta_avg,betaShock,betaT_Eng,betaT,  dateD,SimDate,Tdata,MonthWeekJP,MonthNumber,WeekNumber,fn)
                else
                    %plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaT,betaT_woAR1,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
                    plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaT,betaShock,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
                end
                if figure_save == 1
                    saveas(f,[home 'Figures/' char(pref) '/beta_path' '.png']);
                end
                betaShock_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaShock;
                betaT_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT;
                betaT_woAR1_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT_woAR1;
                beta_tilde_avg = beta_avg*((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k);

                %figure(401)
                figure('Name','BetaTildePath')
                set(gcf,'Position',[100,100,1200,500])
                plot_beta(var_initial,var_share,beta_tilde,beta_tilde_avg,betaT_tilde,betaT_woAR1_tilde,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
                title('β tildeの推移','FontSize',20,'FontWeight','normal','FontName',fn)
            end
        elseif beta_figure_loop == 1
            %figure(400+iTH) % Figure for BetaT with Variant Share
            figure('Name',string(['BetaPath_', sprintf(ft,TH(iTH))]))
            set(gcf,'Position',[100,100,1200,500])
            plot_beta_India( var_initial,var_initial2,var_share,var_shareInd,beta,beta_avg,betaShock,betaT_Eng,betaT,  dateD,SimDate,Tdata,MonthWeekJP,MonthNumber,WeekNumber,fn)
            %plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaT,betaT_woAR1,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
            %plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaShock,betaT,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
            title(string(['βの推移（',sprintf(ft,TH(iTH)),'）']),'FontSize',20,'FontWeight','normal','FontName',fn)

            betaT_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT;
            betaShock_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaShock;
            betaT_woAR1_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT_woAR1;
            beta_tilde_avg = beta_avg*((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k);

            %figure(440 + iTH)
            figure('Name',string(['BetaTildePath_', sprintf(ft,TH(iTH))]))
            set(gcf,'Position',[100,100,1200,500])
            %plot_beta(var_initial,var_share,beta_tilde,beta_tilde_avg,betaT_tilde,betaT_woAR1_tilde,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
            plot_beta(var_initial,var_share,beta_tilde,beta_tilde_avg,betaShock_tilde,betaT_tilde,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
            title(string(['β tildeの推移（',sprintf(ft,TH(iTH)),'）']),'FontSize',20,'FontWeight','normal','FontName',fn)
        end


    end

    %minAlpha = min(minAlphaMat); %minimum alpha when variants have no effects.
    minAlpha = alpha_off; % 経済損失0 = 2020年10-11月のGDP level

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

    for l = 1:2 %1:2 when english version needed
        % Generate graphs for the website
        lng = language{l};
        %figname = 100 + l;
        %figure(figname);
        figname = string(['MainResults_' char(lng)]);
        f = figure('Name',figname);
        set(gcf,'Position',[100,100,1200,800])
        subplot(2,2,1)
        [BackDataN,BackDataAlpha,BackDataERN] = plot_SimN(TH,TH_index,N,NPath,alpha,AlphaPath,ERN,SimERN,THonPath,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,20,fn,ft,l);
        %--- Number of cumulative deaths ---%
        subplot(2,2,2)
        plot_Tradeoff(AlphaM,DM,waves,TH,TH_index,l,linecolor,fs,fn)
        %--- Number of people who are in ICU ---%
        subplot(2,2,3)
        if ICU_nation == 1
            BackDataICU = ...
                plot_ICU_GOV(TH,TH_index,ICU,SimICU,BED(Tdata),...
                MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,...
                linecolor,fs,fn,ft,l)
        else
            BackDataICU = plot_ICU(TH,TH_index,ICU,SimICU,ICU_limit,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,20,fn,ft,l,th_off1);
        end
        %xlim([Tdata-7 Tdata+33])
%         subplot(2,2,4)
%         plot_ICU_N(TH,TH_index,N,NPath,ICU,SimICU,ICU_limit,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,20,fn,ft,l,th_off1/7,th_off2/7,th_off3/7)

%         %--- Plot ICU inflow ---%
%         figname = 'ICU_inflow';
%         f = figure('Name',figname);
%         set(gcf,'Position',[100,100,1200,500])
%         subplot(1,2,1)
%         plot(ICU_inflow, 'LineWidth', 1.5)
%         title('ICU inflow')
%         ytickformat('%,6.0f')
%         xticks(find(WeekNumber==1))
%         xticklabels(MonthWeekJP(WeekNumber==1))
%         lgd.NumColumns = 2;
%         xtickangle(45)
%         subplot(1,2,2)
%         plot(ICU_inflow.*delta, 'LineWidth', 1.5)
%         title('ICU inflow * \delta')
%         xticks(find(WeekNumber==1))
%         xticklabels(MonthWeekJP(WeekNumber==1))
%         lgd.NumColumns = 2;
%         xtickangle(45)
%         xlim([20 Tdata])

        if figure_save == 1
            %saveas(f,[home 'Figures/' char(pref) '/MainResult_' char(lng) '.png']);
            if variant_India == 1
                saveas(f,[home 'Figures/' char(pref) '/' char(fig_title) '_India_' char(lng) '.png']);
            else
                saveas(f,[home 'Figures/' char(pref) '/' char(fig_title) '_' char(lng) '.png']);
            end
            %saveas(figure(figname),[home 'Figures/' char(pref) '/MainResult_' char(lng) '.png']);
        end
    end %End of language loop = figure loop

    if data_save == 1
        titleN = strings(1,1+length(TH_index)*3);
        titleN(1) = "週";
        for ti = 1:length(TH_index)
            titleN(1,1+ti) = string(['新規感染者数（',sprintf('%.0f',TH_index(ti)),'）']);
            titleN(1,1+length(TH_index)+ti) = string(['経済活動（',sprintf('%.0f',TH_index(ti)),'）']);
            titleN(1,1+length(TH_index)*2+ti) = string(['実効再生産数（',sprintf('%.0f',TH_index(ti)),'）']);
            titleN(1,1+length(TH_index)*3+ti) = string(['重症者数（',sprintf('%.0f',TH_index(ti)),'）']);
        end
        TN = table([titleN;MonthWeekJP(Tdata-7:end-1),round(BackDataN(:,1:length(TH_index))/7),round(100*(1-BackDataAlpha(:,1:length(TH_index))),1),round(BackDataERN(:,1:length(TH_index)),2),round(BackDataICU(2:end,1:length(TH_index)))]);
        titleAD = ["経済損失（億円）","死亡者数","ケース"];
        TAD = table([titleAD;BackDataDA(1:length(TH),:)]);
        if variant_India == 1
            writetable(TN,[home 'Figures/' char(pref) '/BackData_' char(data_title) '_India_' char(pref)  '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
            writetable(TAD,[home 'Figures/' char(pref) '/BackData_' char(data_title) '_India_' char(pref) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
        else
            writetable(TN,[home 'Figures/' char(pref) '/BackData_' char(data_title) '_' char(pref)  '.xls'],'Sheet','新規感染者数（1日平均）','WriteVariableNames',false);
            writetable(TAD,[home 'Figures/' char(pref) '/BackData_' char(data_title) '_' char(pref) '.xls'],'Sheet','経済損失と死亡者数','WriteVariableNames',false);
        end
    end

    %figname = 140 + pindex; %Plotting New Cases + Trade Off + Alpha Path
    %figure(figname)
    figname = string(['MainResult+Alpha_' char(pref)]);
    figure('Name',figname)
    set(gcf,'Position',[100,100,1400,600])
    subplot(1,3,1)
    plot_SimN(TH,TH_index,N,NPath,alpha,AlphaPath,ERN,SimERN,THonPath,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,20,fn,ft,l);
    subplot(1,3,2)
    plot_Tradeoff(AlphaM,DM,waves,TH,TH_index,2,linecolor,fs,fn)
    subplot(1,3,3)
    plot_Alpha(alpha,AlphaPath,TH,TH_index,MonthWeekEN,WeekNumber,Tdata,linecolor,ft,fs,fn,2)


%     end
end %end of prefecture loop
