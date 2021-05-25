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
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 1;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 1;      % save back data
vaccine_figure_loop = 0; % =0 appear only once; =1 appear every loop;
beta_figure_loop = 0; % =0 appear only once; =1 appear every loop;
vaccine_disp_switch = 1; % =0 not display the summry of # of the vaccinated ; =1 display
variant_India = 1; % = 0 ... without Indian variant, = 1 ... with Indian variant
% in the "Figure" folder
fs = 16;            % common font size for many figures
ldfs = 8;           % legend font size for vaccine path
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
    %%
    %====================== Model parameter values ======================%
    pref = PrefVector{pindex};        % prefecture to be analyzed
    prefGDP = GDPVector(pindex);
    %====================================================================%
    betaT_temp_ini = betaT_temp_ini_vec(pindex);
    start_beta = start_beta_vec(pindex);
    beta_rho = beta_rho_vec(pindex);
    var_ss = var_ss_vector(pindex);

    th_on1 = th_on_vector(pindex)*7;         % threshold to place the state of emergency
    th_on2 = th_on_vector2(pindex)*7;         % threshold to place the state of emergency
    th_off1 = th_off_vector(pindex)*7;
    th_off2 = th_off_vector2(pindex)*7;
    th_off3 = th_off_vector3(pindex)*7;
    alpha_on = alpha_on_vector(pindex); %(((ERN_on*(POP0/S(end))*((gammaT(1)+delta_average)/beta_avg)).^(1/k))-1)*(h(1)/h(2));
    ICU_limit = ICU_limit_vec(pindex);

    %--- Import data ---%
    % Covid data are recorded at weekly frequencies (Mon-Sun)
    % The first week start on January 20 (Mon), 2020
    import_prefecture

    %--- Construct weekly vaccine data ---%　
    [V1_medical_ori, V1_medical, V2_medical_ori, V2_medical,...
        V1_elderly_ori, V1_elderly, V2_elderly_ori, V2_elderly, ...
        vs_MT,vs_ET,vs_M1,vs_E1,vs_M2,vs_E2] ...
        = ImportVaccineData(home,iPC,Data,pref,dateEN,ps,vaccine_disp_switch);
    %------------------------------------------------------------------------%

     %--- Constructing the reference level of output ---%
    [potentialGDP, referenceGDP, alpha] = construct_GDP(GDP,TdataGDP);

    %--- Regress mobility on alpha to estimate the elasticity h ---%　関数化？
    [Malt,h_all,h_all_se,h,h_se] = estimate_h(M,alpha,TdataGDP,RetroH,hconstant);

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
    ICU(2:Tdata+1,1) = Data(:,21);
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
        alpha_on = 0.5*alpha_May+0.5*alpha_Jan;
        TH = [100:100:700];
        TH_index = [400,500,700];
    elseif pindex == 5
        alpha_on = 1*alpha_May;
        TH = [200:100:800];
        TH_index = [500,600,700];
    end

    if variant_India == 1
%         DRi = 10;
        TH_index = [300,500,700];
    end
    
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

        %--- Compute the history of time-varying parameters ---%
        delta = (D(2:Tdata+1)-D(1:Tdata))./I(1:Tdata);                              % death rate
        beta_tilde = (POP0.*N(1:Tdata))./((S(1:Tdata).*I(1:Tdata)));   % overall infection rate
        ERN = (S(1:end-1)/POP0).*beta_tilde./(gamma+delta);                                        % effective reproduction number
        if hconstant == 0
            beta = beta_tilde./(1+h_all*alpha).^k;                                      % raw infection rate
        elseif hconstant == 1
            beta = beta_tilde./(1+(h_all(2)/h_all(1))*alpha).^k;
        end
        ICU_inflow = (ICU(2:Tdata+1) - ICU(1:Tdata) + gamma_ICU.*ICU(1:Tdata) + dD(1:Tdata))./(delta(1:Tdata).*N(1:Tdata));

        %--- Construct time series of parameters ---%
        gammaT = gamma*ones(SimPeriod,1);
        delta_sample = delta(end-RetroPeriod+1:end);
        delta_average = sum(delta_sample.*(I(end-RetroPeriod+1:end)/sum(I(end-RetroPeriod+1:end))));
        ICU_inflow_avg = mean(ICU_inflow(end-RetroPeriod+1:end))*ICU_adjustment;

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
        elderly = elderly - (sum(V1_elderly));
        [V,deltaT,VT] = vaccine_distribution_medical(vacpath,medical,V1_medical,V2_medical,elderly,V1_elderly,V2_elderly, ordinary,elderly_total,delta_average,E1,E2,D1,D2,ps,POP0,lag,medical_duration);
        %         deltaT = construct_delta_variant(home,SimPeriod,RetroPeriod,Tdata,dateEN,pref,retro_lb,retro_ub,deltaT,delta,delta_average,var_initial,var_ss,var_infection_delta,var_growth,I);

        [var_share, var_prev,var_initial] = var_share_prev(Data,SimPeriod,var_ss,var_growth);
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

        [betaT,betaT_woAR1] = construct_beta(SimPeriod,retro_lb,retro_ub,...
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
            betaT = betaT.*var_boost;
        end


        betaT = beta_AR1(betaT_temp_ini, beta_rho, betaT, start_beta);

        alpha_off = mean(alpha((dateEN >= datetime(2020,2,7)) & (datetime(2020,2,28)>= dateEN ))); % output loss without the state of emergency
        InitialValues = [S(end),I(end),R(end),D(end),ICU(end)];


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [DMat(iTH),AlphaMat(iTH),AlphaPath(:,iTH),SimData(:,:,iTH),NPath(:,iTH),SimERN(:,iTH),THonPath(:,iTH),SimICU(:,iTH)] ...
                = Covid_projection_ICU4...
                (InitialValues,alpha_on,alpha_off,th_on1,th_on2,th_off1,th_off2,th_off3,betaT,gammaT,deltaT,deltaT(1),V,h,k,POP0,hconstant,DRi,state,ICU_inflow_avg,gamma_ICU);

        % Plot betaT
        if beta_figure_loop == 0
            if iTH == 1
                %figure(400) % Figure for BetaT with Variant Share
                f=figure('Name','BetaPath');
                set(gcf,'Position',[100,100,1200,500])
                if variant_India == 1
                    plot_beta_India(var_initial,var_initial2,var_share,var_shareInd,beta,beta_avg,betaT,betaT_Eng,betaT_woAR1,dateD,SimDate,Tdata,MonthWeekJP,MonthNumber,WeekNumber,fn)
                else
                    plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaT,betaT_woAR1,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
                end
                if figure_save == 1
                    saveas(f,[home 'Figures/' char(pref) '/beta_path' '.png']);
                end

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
            plot_beta(var_initial,var_share,beta,beta_avg*ones(length(SimDate),1),betaT,betaT_woAR1,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
            title(string(['βの推移（',sprintf(ft,TH(iTH)),'）']),'FontSize',20,'FontWeight','normal','FontName',fn)

            betaT_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT;
            betaT_woAR1_tilde = ((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k).*betaT_woAR1;
            beta_tilde_avg = beta_avg*((1+(h(2)/h(1)).*AlphaPath(:,iTH)).^k);

            %figure(440 + iTH)
            figure('Name',string(['BetaTildePath_', sprintf(ft,TH(iTH))]))
            set(gcf,'Position',[100,100,1200,500])
            plot_beta(var_initial,var_share,beta_tilde,beta_tilde_avg,betaT_tilde,betaT_woAR1_tilde,dateD,SimDate,Tdata,tt,MonthWeekJP,MonthNumber,WeekNumber,fn)
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
        BackDataICU = plot_ICU(TH,TH_index,ICU,SimICU,ICU_limit,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,20,fn,ft,l,th_off1);
        xlim([Tdata-7 Tdata+33])
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
