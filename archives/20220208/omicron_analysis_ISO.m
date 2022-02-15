clear variables
close all
iPC = 0; % 0 for Mac, 1 for Windows

if iPC == 1
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
    fn = 'Yu Gothic'; % Font style for xaxis, yaxis, title
else
%         home = '/Users/sohtakawawaki/Dropbox/fujii_nakata (1)/Website/Codes/';
    %     home = '/Users/okamotowataru/Dropbox/fujii_nakata/Website/Codes/';
    home = '/Users/ymaeda/Documents/ym_doc/Tokyo_Univrsity_MA/Reserach_Assistant/fujii_nakata/Codes/required_files/';
%     home = '/Users/shotaro/Dropbox/fujii_nakata/Website/Codes/';
    fn = 'YuGothic';
end

cd(home);
%====================== Program parameter values ======================%
pref = 'Tokyo';
prefGDP = 106;
figure_save = 1; % 0 = figures won't be saved, 1 = they will be saved
data_save = 0; % save back data

exogenous_omicron_share = 0;
alternative_scenario = 0; %=0 ... x loop is for relative BRN, =1 ... x loop is for relative VE
data_switch = 0; % Use I_data and gamma = mean(gamma_data(end-17+1:end))

fs = 12; % common font size for many figures
ldfs = 12; % legend font size for vaccine path
lgdfs = 10;
ldfs_main = 12;
axfs = 10;
ft = '%.1f';
yft = '%.0f';
language = {'EN', 'JP'};

%===================== Figure Names =================%
figname_main = string(['MainResults']);
figname_beta_tilde = 'Beta Tilde Path';
figname_beta = 'Beta Path';
figname_var = 'Variant_Share';
figname_delta = 'Death Rate transition';
figname_ICU_nation = 'ICU Rate transition National Standard';
figname_ICU_local = 'ICU Rate transition Local Standard';
figname_ERN = 'ERN transition';
figname_Hospital = 'Hospital Rate transition';
figname_dD = 'dD transition';
figname_BRN = 'BRN';
figname_omicron_share = 'Omicron share';


%================== Model Fixed Parameter Values ============================%
parameter
SimPeriod = 15;

%=============== Import data ============%
import_prefecture
YearMonth = [YearMonthEN, YearMonthJP];
YearMonthWeek = [YearMonthWeekEN, YearMonthWeekJP];
%Date and Figure parameter
% ind_first_omicron   = find(date ==  datetime(2022, 1, 13));
ind_first_omicron   = find(date ==  datetime(2021, 12, 23));
x_left_omicron = find(date == datetime(2021, 12, 9));
x_right_omicron = find(date == datetime(2022, 3, 3));

% ICU_nation(end) = 12; %Update Saturday values from Tokyo website : https://www.fukushihoken.metro.tokyo.lg.jp/iryo/kansen/corona_portal/info/kunishihyou.html
if data_switch == 1
    I = I_data;
    gamma_data = zeros(Tdata, 1);
    %     delta_data = zeros(Tdata,1);
    for i = 1:Tdata
        if I(i) > 0
            gamma_data(i) = dR_data(i) / I_data(i);
            %             delta_data(i) = dD_data(i)/I_data(i);
        end
    end
    gamma = mean(gamma_data(end - RetroPeriod + 1:end));
end


%============ Import vaccine data ============%
Vmat = readmatrix([home 'vaccination_Tokyo_newly.xlsx']);
vac1stweek = Vmat(1, 1); %column 1 = week, raw 1 = first week

V1_elderly = zeros(Tdata, 1);
V2_elderly = V1_elderly;
V1_others  = V1_elderly;
V2_others  = V1_elderly;
V1_elderly(vac1stweek:Tdata)    = Vmat(:, 4);
V2_elderly(vac1stweek:Tdata)    = Vmat(:, 5);
V1_others(vac1stweek:Tdata)     = Vmat(:, 7);
V2_others(vac1stweek:Tdata)     = Vmat(:, 8);

% Medical personels
vaccine_medical = readmatrix([home 'vaccine_daily_medical.xls']);
[V1_medical_past, V2_medical_past] = vaccine_daily_to_weekly_table(vaccine_medical, ps, dateEN, iPC);

M_first = cum_to_new(Data(:, 12));  M_second = cum_to_new(Data(:, 13));
M_ps_first = Data(:, 15);   M_ps_second = Data(:, 16);
M_first(find(date == datetime(2021, 8, 5)):end) = zeros(length(find(date == datetime(2021, 8, 5)):length(M_first)),1);
M_second(find(date == datetime(2021, 8, 5)):end) = zeros(length(find(date == datetime(2021, 8, 5)):length(M_second)),1);
indM1 = find(M_ps_first > 0, 1, 'first');   indM21 = find(M_ps_second > 0, 1, 'first');

V1_medical = (V1_medical_past / ps) * M_ps_first(indM1);
V1_medical(indM1 + 1:end) = M_first(indM1 + 1:end);
V2_medical = (V2_medical_past / ps) * M_ps_second(indM1);
V2_medical(indM21 + 1:end) = M_second(indM21 + 1:end);
cumsumPastV1 = cumsum(V1_elderly + V1_others + V1_medical);
cumsumPastV2 = cumsum(V2_elderly + V2_others + V2_medical);

%============ Simulated Vaccine Path ============%
% VT = zeros(SimPeriod,6);
% VT(1,2) = V1_elderly(end-2);
% VT(2,2) = V1_elderly(end-1);
% VT(3,2) = V1_elderly(end);
% VT(1,4) = V1_others(end-2);
% VT(2,4) = V1_others(end-1);
% VT(3,4) = V1_others(end);
% cumsumVT1           = cumsum(VT(:,1) + VT(:,3) + VT(:,5))+ cumsumPastV1(end);
% lagged_cumsumVT1    = [cumsumPastV1(end-2);cumsumPastV1(end-1);cumsumVT1(1:end-2)];
% cumsumVT2           = cumsum(VT(:,2) + VT(:,4) + VT(:,6))+ cumsumPastV2(end);
% lagged_cumsumVT2    = [cumsumPastV2(end-2);cumsumPastV2(end-1);cumsumVT2(1:end-2)];

    VT = zeros(SimPeriod,9);
    VT(1,2) = V1_elderly(end-2);
    VT(2,2) = V1_elderly(end-1);
    VT(3,2) = V1_elderly(end);
    VT(1,5) = V1_others(end-2);
    VT(2,5) = V1_others(end-1);
    VT(3,5) = V1_others(end);
    
    V3_elderly  = zeros(Tdata+SimPeriod,1);
    aa          = V2_elderly(find(V2_elderly>0,1,'first'):end);
    V3_elderly(indJan2022:length(aa)+indJan2022-1) = aa;
    V3_elderly(length(aa)+indJan2022:length(aa)+indJan2022+2) = VT(1:3,2);
    VT(:,3) = V3_elderly(Tdata+1:Tdata+SimPeriod);
    V3_elderly = V3_elderly(1:Tdata,1);
    
    V3_medical  = zeros(Tdata+SimPeriod,1);
    aa          = V2_medical(find(V2_medical>0,1,'first'):end);
    V3_medical(indDec2021:length(aa)+indDec2021-1) = aa;
    V3_medical(length(aa)+indDec2021:length(aa)+indDec2021+2) = VT(1:3,8);
    VT(:,9) = V3_medical(Tdata+1:Tdata+SimPeriod);
    V3_medical = V3_medical(1:Tdata,1);
    
    V3_others  = zeros(Tdata+SimPeriod,1);
    aa          = V2_others(find(V2_others>0,1,'first'):end);
    V3_others(indJan2022:length(aa)+indJan2022-1) = aa;
    V3_others(length(aa)+indJan2022:length(aa)+indJan2022+2) = VT(1:3,2);
    VT(:,6) = V3_others(Tdata+1:Tdata+SimPeriod);
    V3_others = V3_others(1:Tdata,1);
    
    cumsumVT1           = cumsum(VT(:,1) + VT(:,4) + VT(:,7))+ cumsumPastV1(end);
    lagged_cumsumVT1    = [cumsumPastV1(end-1);cumsumPastV1(end);cumsumVT1(1:end-2)];
    cumsumVT2           = cumsum(VT(:,2) + VT(:,5) + VT(:,8))+ cumsumPastV2(end);
    lagged_cumsumVT2    = [cumsumPastV2(end-1);cumsumPastV2(end);cumsumVT2(1:end-2)];
    
    cumsumPastV3        = cumsum(V3_elderly + V3_others + V3_medical);
    cumsumVT3           = cumsum(VT(:,3) + VT(:,6) + VT(:,9))+ cumsumPastV3(end);
    lagged_cumsumVT3    = [cumsumPastV3(end-1);cumsumPastV3(end);cumsumVT3(1:end-2)];

%============ Constructing the reference level of output ===========%
[potentialGDP, referenceGDP, alpha] = construct_GDP(GDP, TdataGDP);
[Malt, h_all, h_all_se, h, h_se] = estimate_h(M, alpha, TdataGDP, RetroH, hconstant);
% figname = string(['Mobility_GDP_' char(pref)]);
% f = figure('Name', figname);
% plot_mobility(Malt, alpha, Tdata, TdataGDP, YearMonthWeekJP, xtick1, fs, 16)
% if figure_save == 1
%     saveas(f, [home 'Figures/' char(pref) '/MobilityGDPLine_v.png']);
% end

%====== Common exogenous settings ====%


% Seasonality
seasonality = seasonal_adjustment(retro_lb, retro_ub, dateD, SimPeriod+DRi+1, seasonal_effect);

%===== Compute the history of S, I, R, D in the data period ====%
if data_switch == 0
    [S, I, R, D,cum_in_R] ...
        = SIRD(Tdata, POP0, N, E1, E2, E3,...
        V1_elderly, V1_medical, V1_others, V2_elderly, V2_medical, V2_others, ...
        V3_elderly,V3_medical,V3_others,...
        gamma, dD, TdataGDP, referenceGDP, alpha);
else
    pastV = zeros(Tdata, 1);
    pastV(3:end) = E1 * (V1_elderly(1:end - 2) + V1_medical(1:end - 2) + V1_others(1:end - 2)) ...
                 + (E2 - E1) * (V2_elderly(1:end - 2) + V2_medical(1:end - 2) + V2_others(1:end - 2)) ...
                 + (E3 - E2) * (V3_elderly(1:end - 2) + V3_medical(1:end - 2) + V3_others(1:end - 2)) ;
    S = S_data;
    S(2:end) = S(2:end) - cumsum(pastV);
    R = R_data;
    R(2:end) = R(2:end) + cumsum(pastV);
    D = D_data;
end

%========== Compute time series average ==========%
[delta, beta_tilde, ERN, beta, ...
    ICU_nation_inflow, ICU_pref_inflow, Hospital_inflow, ...
    gammaT, delta_average, delta_ICU_nation_average, delta_ICU_pref_average, delta_Hospital_average, ...
    ICU_nation_inflow_avg, ICU_pref_inflow_avg, Hospital_inflow_avg, ...
    simple_beta_avg, beta_se, ...
    delta_se, delta_ICU_nation_se, delta_ICU_pref_se, delta_Hospital_se] ...
    = Time_Series_Average(S, I, D, ICU_nation, ICU_pref, hospital, dD, N, ...
    Tdata, SimPeriod, RetroPeriod, POP0, ...
    hconstant, h_all, alpha, k, ...
    gamma, gamma_ICU_nation, gamma_ICU_pref, gamma_Hospital, ...
    ICU_nation_adjustment, ICU_pref_adjustment, Hospital_adjustment, ...
    RetroPeriodDelta, RetroPeriodICU_nation, RetroPeriodICU_pref, RetroPeriodHospital);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Simulaiton starts here %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=========== simulation setting =============%
nVariant = 2;
alpha_Aug = mean(alpha((dateEN >= datetime(2021, 8, 05)) & (datetime(2021, 8, 26) >= dateEN)));
alpha_Oct2021 = mean(alpha((dateEN >= datetime(2021, 10, 07)) & (datetime(2021, 10, 28) >= dateEN)));
% exogenous parameters
state               = 0;    %=1 if currently under SOE
th_on               = ICU_limit_pref_vec(end)*0.5; %8000*7; %ICU_limit_pref_vec(end)*0.5; %
th_off              = ICU_limit_pref_vec(end)*0.25; %1000 * 7;  %ICU_limit_pref_vec(end)*0.25; %
ori_beta_goal           = 2.1984; %3.75 2.0500; %3.5
if exogenous_omicron_share == 1
    relative_beta_SOE   = 1.2;
else
    relative_beta_SOE   = 0.25;
end
alpha_on            = alpha_Oct2021;
alpha_off           = 0.00;
alpha_jump          = 0.9 * alpha_Aug;
beta_jump           = 1.0; %0.6

InitialValues = [S(end), I(end), R(end), D(end), ICU_nation(end), ICU_pref(end), hospital(end)];

xmin = find(date == datetime(2021, 7, 1));
xmax = find(date == datetime(2022, 5, 26));
ymax_D = 50;

% Simulation parameters
nZ = 111;
omicron_relative_infectivity_vector = linspace(1,2,101);
omicron_E2_vector = linspace(0,1,101);    %[1, 0.6, 0.2]; %Relative
omicron_realtive_severity_vector = linspace(0,1.1,nZ);

% Directory for saved figures
if alternative_scenario == 0
%     xvec = omicron_relative_infectivity_vector;
    xvec = 1.25;
%     xvec = 1.0;
    target_ICU_vec = [75, 85, 100, 200, 500, 1000];
    yvec = omicron_E2_vector;
    xmin_ISO = 0; xmax_ISO = 1;
    zvec = omicron_realtive_severity_vector;
    nX = length(xvec);
    nY = length(yvec);
    nZ = length(zvec);
    
    
    figfolder = string(['Omicron_Analysis_Relative_Infectivity']);
    % Scenario Name and Line
    titlevec = {'relative BRN = ', '相対的基本再生産数　＝　'};
    figname_xvar = '_relativeBRN_';
elseif alternative_scenario == 1
    target_ICU_vec = [100, 152, 200, 500, 1000];
    xvec = 0.6;
    yvec = omicron_relative_infectivity_vector;
    xmin_ISO = 1; xmax_ISO = 2;
    zvec = omicron_realtive_severity_vector;
    nX = length(xvec);
    nY = length(yvec);
    nZ = length(zvec);
    figfolder = string(['Omicron_Analysis_Relative_VE']);
    titlevec = {'relative VE = ', '相対的ワクチン有効性　＝　'};
    figname_xvar = '_relativeVE_';
end
column_num_main = nY;

originalE1 = E1;
originalE2 = E2;
originalE3 = E3;

% initialize_matirix
DMat        = nan(nX, nY, nZ);
AlphaMat    = DMat;

SimData                 = nan(SimPeriod + 1, length(InitialValues), nX, nY, nZ);
SimData_endogenous      = zeros(SimPeriod + 1, length(InitialValues), nX, nY, nZ, nVariant);
beta_path_mat           = zeros(SimPeriod, nX, nY, nZ, nVariant);
beta_tilde_path_mat     = beta_path_mat;
BRN_path_mat            = beta_path_mat;
ERN_path_mat            = beta_path_mat;
AlphaPath               = nan(SimPeriod, nX, nY, nZ);
NPath                   = AlphaPath;
simD                    = AlphaPath;
SimERN                  = AlphaPath;
SimBRN                  = AlphaPath;
betaPath                = AlphaPath;
betaTildePath           = AlphaPath;
deltaPath               = AlphaPath;
delta_ICU_nationPath    = AlphaPath;
delta_ICU_prefPath      = AlphaPath;
delta_HospitalPath      = AlphaPath;
omicron_share_mat       = AlphaPath;
omicron_I_share_mat     = AlphaPath;
omicron_N_share_mat     = AlphaPath;
SimICU_nation           = nan(SimPeriod+1, nX, nY, nZ);
SimICU_pref             = SimICU_nation;
SimHospital             = SimICU_nation;
tic
for iX = 1:nX    

    for iY = 1:nY
        if alternative_scenario == 0
            omicron_relative_infectivity    = xvec(iX);
            relative_infectivity_path       = ones(SimPeriod,1); %vector
            omicronE3 = originalE3;
            omicronE2 = originalE2 * yvec(iY);
            omicronE1 = originalE1 * yvec(iY);
        elseif alternative_scenario == 1
            omicronE3 = originalE3;
            omicronE2 = originalE2 * xvec(iX);
            omicronE1 = originalE1 * xvec(iX);
            omicron_relative_infectivity    = yvec(iY);
            relative_infectivity_path       = ones(SimPeriod,1); %vector
        end
        
        for iZ = 1:nZ
            omicron_realtive_severity = omicron_realtive_severity_vector(iZ);
            relative_severity_path    = ones(SimPeriod,1); %vector
            past_omicron_share  = zeros(Tdata,1);
            if exogenous_omicron_share == 1
                %calculate_omicron_share
                logit_initial       = log(omicron_initial/(omicron_ss-omicron_initial)); % Logit of the variant share, most recently
                sim_omicron_share   = zeros(SimPeriod,1);
                sim_omicron_share(ind_first_omicron-Tdata)  =   omicron_initial;
                sim_omicron_share(ind_first_omicron-Tdata+1:end,1)  ...
                    = exp((1:length(sim_omicron_share(ind_first_omicron-Tdata+1:end,1)))'* omicron_growth + logit_initial).*omicron_ss ...
                    ./(1+exp((1:length(sim_omicron_share(ind_first_omicron-Tdata+1:end,1)))'*omicron_growth+logit_initial));
                sim_omicron_share(isnan(sim_omicron_share)) = 1;
                omicron_share       = [past_omicron_share; sim_omicron_share];
                %plot omicron share
                if iX == 1 && iY == 1 && iZ == 1
                    figure('Name', char(figname_omicron_share));
                    set(gcf, 'Position', [100, 100, 1200, 800])
                    jTitle = 'オミクロン株割合の推移';
                    plot(omicron_share)
                    xlim([x_left_omicron x_right_omicron])
                    ylim([0 1])
                    xticks(x_left_omicron:x_right_omicron)
                    xticklabels(YearMonthWeekJP(xticks))
                    title(char(jTitle),'FontSize',fs,'FontWeight','normal','FontName',fn)
                    ax = gca;
                    ax.YAxis.FontSize = axfs;   ax.XAxis.FontSize = axfs;   ax.YAxis.Exponent = 0;
                    xtickangle(45)
                    if figure_save == 1
                        saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_omicron_share) '.png']);
                    end
                end
                
                % calculate vaccine effectiveness
                E1 = originalE1*(1-sim_omicron_share) + omicronE1 * sim_omicron_share; %vector
                E2 = originalE2*(1-sim_omicron_share) + omicronE2 * sim_omicron_share; %vector
                E3 = originalE3*(1-sim_omicron_share) + omicronE3 * sim_omicron_share; %vector
                
                relative_infectivity_path       = 1 * (1-sim_omicron_share) ...
                    + omicron_relative_infectivity * sim_omicron_share; %vector
                
                relative_severity_path    = 1 * (1-sim_omicron_share) ...
                    + omicron_realtive_severity * sim_omicron_share; %vector
            end
            VE = E1.*(VT(:,1)+VT(:,4)+VT(:,7))...
               +(E2-E1).*(VT(:,2)+VT(:,5)+VT(:,8)) ...
               +(E3-E2).*(VT(:,3)+VT(:,6)+VT(:,9));
            VE_prev = originalE1.*(V1_elderly+V1_medical+V1_others)...
                    +(originalE2-originalE1).*(V2_elderly+V2_medical+V2_others)...
                    +(originalE3-originalE2).*(V3_elderly+V3_medical+V3_others);
            V = [VE_prev(end-1);VE_prev(end);VE(1:end-2)];
           
            VE_omicron      = omicronE1 * (VT(:,1)+VT(:,4)+VT(:,7))...
                + (omicronE2-omicronE1) * (VT(:,2)+VT(:,5)+VT(:,8)) ...
                + (omicronE3-omicronE2) * (VT(:,3)+VT(:,6)+VT(:,9));
            VE_prev_omicron = omicronE1 * (V1_elderly+V1_medical+V1_others) ...
                + (omicronE2-omicronE1) * (V2_elderly+V2_medical+V2_others) ...
                + (omicronE3-omicronE2) * (V3_elderly+V3_medical+V3_others);
            V_omicron       = [VE_prev_omicron(end-1); VE_prev_omicron(end); VE_omicron(1:end-2)];
            
            
            %Construct betapath
            beta_goal   = ori_beta_goal;
            beta_goal   = beta_goal * relative_infectivity_path;
            betaT       = beta_goal .* transpose(seasonality(1:SimPeriod));
            betaT_woAR             = betaT;
            betaT                  = beta_AR1(betaT_temp_ini, beta_rho, betaT, start_beta);
            betaBox     = [beta; betaT];
            
            %Construct alphapath
            alphaAfterSOE   = [alpha(Tdata):(alpha_off - alpha(Tdata)) / (DRi):alpha_off];
            alphaT          = [alphaAfterSOE'; alpha_off * ones(SimPeriod, 1)];
            alphaBox        = [alpha; alphaT];
            
            % death rate, severity rate, hospoital rate paths
            deltaT              = delta_average             * relative_severity_path;
            delta_ICU_nation    = delta_ICU_nation_average  * relative_severity_path;
            delta_ICU_pref      = delta_ICU_pref_average    * relative_severity_path;
            delta_Hospital      = delta_Hospital_average    * relative_severity_path;
            
            %AR1 adjustment
            deltaT_woAR             = deltaT;
            deltaT                  = beta_AR1(delta_temp_ini, delta_rho, deltaT, start_delta);
            
            delta_ICU_nation_woAR   = delta_ICU_nation;
            delta_ICU_nation        = beta_AR1(delta_ICU_nation_temp_ini, delta_ICU_nation_rho, delta_ICU_nation, start_delta_ICU_nation);
            
            delta_ICU_pref_woAR     = delta_ICU_pref;
            delta_ICU_pref          = beta_AR1(delta_ICU_pref_temp_ini, delta_ICU_pref_rho, delta_ICU_pref, start_delta_ICU_pref);
            
            delta_Hospital_woAR     = delta_Hospital;
            delta_Hospital          = beta_AR1(delta_Hospital_temp_ini, delta_Hospital_rho, delta_Hospital, start_delta_Hospital);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation SIRD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     simulationSIRD_(t,loop1,loop2,loop3)
            %             [DMat(1,iX,iY,iZ), AlphaMat(1,iX,iY,iZ), AlphaPath(:, iX,iY,iZ), ...
            %                 SimData(:, :, iX,iY,iZ), NPath(:, iX,iY,iZ), SimERN(:, iX,iY,iZ),...
            %                 SimICU_nation(:, iX,iY,iZ), SimICU_pref(:, iX,iY,iZ), SimHospital(:, iX,iY,iZ), betaPath(:,iX,iY,iZ)] ...
            %                 = Covid_projection_withSOE(InitialValues, alpha_on, alpha_off, th_on, th_off, ...
            %                 betaT, gammaT, deltaT, delta_ICU_nation, delta_ICU_pref, delta_Hospital, V, h, k, POP0, ...
            %                 hconstant, DRi, ICU_nation_inflow_avg, ICU_pref_inflow_avg, Hospital_inflow_avg, ...
            %                 gamma_ICU_nation, gamma_ICU_pref, gamma_Hospital, ...
            %                 relative_beta_SOE, beta_jump, beta_goal(1), seasonality, alphaBox(Tdata + 1:end), state, simple_beta_avg);
            if exogenous_omicron_share == 1
                [DMat(iX,iY,iZ), AlphaMat(iX,iY,iZ), AlphaPath(:, iX,iY,iZ), ...
                    SimData(:, :, iX,iY,iZ), NPath(:, iX,iY,iZ), SimERN(:, iX,iY,iZ),...
                    SimICU_nation(:, iX,iY,iZ), SimICU_pref(:, iX,iY,iZ), SimHospital(:, iX,iY,iZ), ...
                    betaPath(:,iX,iY,iZ),betaTildePath(:,iX,iY,iZ),SimBRN] = ...
                    Covid_projection_omicron_ISO(InitialValues,alpha_on,alpha_off,th_on,th_off,...
                    betaT,gammaT,deltaT,delta_ICU_nation,delta_ICU_pref,delta_Hospital,h,k,POP0,...
                    lagged_cumsumVT1, lagged_cumsumVT2, E1, E2, cum_in_R, ...
                    hconstant,DRi,ICU_nation_inflow_avg, ICU_pref_inflow_avg, Hospital_inflow_avg, ...
                    gamma_ICU_nation, gamma_ICU_pref, gamma_Hospital,...
                    relative_beta_SOE, beta_jump, beta_goal, seasonality, alphaBox(Tdata + 1:end), state, simple_beta_avg);
                simD(:, iX, iY, iZ) = squeeze(SimData(2:end,4,iX,iY,iZ)-SimData(1:end-1,4,iX,iY,iZ)+SimData(2:end,4,iX,iY,iZ)-SimData(1:end-1,4,iX,iY,iZ));
                
                deltaPath(:, iX,iY,iZ)  = deltaT;
                delta_ICU_nationPath(:, iX,iY,iZ) = delta_ICU_nation;
                delta_ICU_prefPath(:, iX,iY,iZ) = delta_ICU_pref;
                delta_HospitalPath(:, iX,iY,iZ) = delta_Hospital;
                
            elseif exogenous_omicron_share == 0
                SimData_endogenous(1,:,iX,iY,iZ,1) = InitialValues;
                SimData_endogenous(1,3,iX,iY,iZ,2) ...
                    = ( omicronE1 * lagged_cumsumVT1(1) ...
                    +   (omicronE2 - omicronE1) * lagged_cumsumVT2(1)       ...
                    +   (omicronE3 - omicronE2) * lagged_cumsumVT3(1)     )...
                    + omicron_immunity * cum_in_R;
                SimData_endogenous(1,1,iX,iY,iZ,2) ...
                    = (POP0- D(end)-iniI_omicron) - SimData_endogenous(1,3,iX,iY,iZ,2) ;
                SimData_endogenous(ind_first_omicron-Tdata+1,2,iX,iY,iZ,2) = iniI_omicron;
                POP0_omicron = (POP0- D(end));
                
                [   DMat(iX,iY,iZ), AlphaMat(iX,iY,iZ), AlphaPath(:, iX,iY,iZ), ...
                    SimData_endogenous(:, :, iX,iY,iZ,:), NPath(:, iX,iY,iZ), ERN_path_mat(:, iX,iY,iZ,:), ...
                    SimICU_nation(:, iX,iY,iZ), SimICU_pref(:, iX,iY,iZ), SimHospital(:, iX,iY,iZ), ...
                    beta_path_mat(:,iX,iY,iZ,:), beta_tilde_path_mat(:,iX,iY,iZ,:), ...
                    BRN_path_mat(:,iX,iY,iZ,:), omicron_share_mat(:, iX,iY,iZ)      ] ...
                    = Covid_projection_endogenous_omicron_ISO(SimData_endogenous(:, :, iX,iY,iZ,:),...
                    alpha_on,alpha_off,th_on,th_off,betaT,gammaT,...
                    deltaT,delta_ICU_nation,delta_ICU_pref,delta_Hospital,...
                    V,V_omicron,h,k,POP0,...
                    omicron_relative_infectivity, omicron_realtive_severity, omicron_immunity,  ...
                    hconstant,DRi,ICU_nation_inflow_avg, ICU_pref_inflow_avg, Hospital_inflow_avg, ...
                    gamma_ICU_nation, gamma_ICU_pref, gamma_Hospital, ...
                    relative_beta_SOE, beta_jump, beta_goal, seasonality, alphaBox(Tdata + 1:end), ...
                    state, simple_beta_avg);
                
                omicron_N_share_mat(:,iX,iY,iZ) = omicron_share_mat(:,iX,iY,iZ);
                omicron_I_share_mat(:,iX,iY,iZ) = SimData_endogenous(1:end-1,2,iX,iY,iZ,2)...
                                    ./(SimData_endogenous(1:end-1,2,iX,iY,iZ,1)+SimData_endogenous(1:end-1,2,iX,iY,iZ,2));
                omicron_share_mat   = omicron_I_share_mat;
                
                simD(:, iX, iY, iZ) = squeeze(SimData_endogenous(2:end,4,iX,iY,iZ,1)-SimData_endogenous(1:end-1,4,iX,iY,iZ,1)+SimData_endogenous(2:end,4,iX,iY,iZ,2)-SimData_endogenous(1:end-1,4,iX,iY,iZ,2));
                betaPath(:,iX,iY,iZ) = beta_path_mat(:,iX,iY,iZ,1) .* (1 - omicron_share_mat(:, iX,iY,iZ)) ...
                                     + beta_path_mat(:,iX,iY,iZ,2) .* omicron_share_mat(:, iX,iY,iZ);
                betaTildePath(:,iX,iY,iZ)   = beta_tilde_path_mat(:,iX,iY,iZ,1) .* (1 - omicron_share_mat(:, iX,iY,iZ)) ...
                                            + beta_tilde_path_mat(:,iX,iY,iZ,2) .* omicron_share_mat(:, iX,iY,iZ);
                SimBRN(:,iX,iY,iZ) = BRN_path_mat(:,iX,iY,iZ,1) .* (1 - omicron_share_mat(:, iX,iY,iZ)) ...
                                    + BRN_path_mat(:,iX,iY,iZ,2) .* omicron_share_mat(:, iX,iY,iZ);
                SimERN(:,iX,iY,iZ)  = ERN_path_mat(:,iX,iY,iZ,1) .* (1 - omicron_share_mat(:, iX,iY,iZ)) ...
                                    + ERN_path_mat(:,iX,iY,iZ,2) .* omicron_share_mat(:, iX,iY,iZ);
                
                deltaPath(:, iX,iY,iZ)  = deltaT(:,1) .* (1 - omicron_share_mat(:, iX,iY,iZ)) ...
                                    + deltaT(:,1) .* omicron_share_mat(:, iX,iY,iZ) * omicron_realtive_severity;
                delta_ICU_nationPath(:, iX,iY,iZ) ...
                    = delta_ICU_nation(:,1) .* (1 - omicron_share_mat(:, iX,iY,iZ)) ...
                    + delta_ICU_nation(:,1) .* omicron_share_mat(:, iX,iY,iZ) * omicron_realtive_severity;
                delta_ICU_prefPath(:, iX,iY,iZ) ... 
                    = delta_ICU_pref(:,1) .* (1 - omicron_share_mat(:, iX,iY,iZ)) ...
                    + delta_ICU_pref(:,1) .* omicron_share_mat(:, iX,iY,iZ) * omicron_realtive_severity;
                delta_HospitalPath(:, iX,iY,iZ) ... 
                    = delta_Hospital(:,1) .* (1 - omicron_share_mat(:, iX,iY,iZ)) ...
                    + delta_Hospital(:,1) .* omicron_share_mat(:, iX,iY,iZ) * omicron_realtive_severity;
           
            end
            
        end
        
    end
    
end
toc
%% Construct matrix for ISO curves
x_mat = yvec;
y_mat = zvec;
z_mat = squeeze(SimICU_pref(end,1,:,:));
columntable = [x_mat; z_mat'];

TICU = table(  [    ...
                [nan(1,1); y_mat'], columntable  ...
                ]  );
nT = length(target_ICU_vec);
ind_target = zeros(nY,nT);

for ii = 1:nT
    for yy = 1:nY
        target_ICU = target_ICU_vec(ii);
        ind_target(yy,ii) = find(abs(z_mat(yy,:) - target_ICU) == min(abs(z_mat(yy,:) - target_ICU )));
    end
end
if alternative_scenario == 0
    for ii = 1:nT
        ind = find(ind_target(:,ii) == nY,1,'first');
        if ind < nY
            ind_target(ind+1:end,ii) = nan;
        end
    end
elseif alternative_scenario == 1
    for ii = 1:nT
        ind = find(ind_target(:,ii) == nY,1,'last');
        if ind > 1
            ind_target(1:ind-1,ii) = nan;
        end
    end
end
if alternative_scenario == 0
    title_name = {'Iso-ICU curve (relative BRN = ', '等重症者数曲線 (相対的BRN = '};
    x_name = {'Relative VE', '相対的VE'};
    save_name = {'Iso_curve_fixed_BRN'};
    lgdloc = 'northwest';
elseif alternative_scenario == 1
    title_name = {'Iso-ICU curve (relative VE = ', '等重症者数曲線 (相対的VE = '};
    x_name = {'Relative BRN', '相対的BRN'};
    save_name = 'Iso_curve_fixed_VE';
    lgdloc = 'northeast';
end
y_name = {'Relative Severity', '相対的重症化率'};
z_name = {'# of Severe Cases (End of March, Tokyo Standard) = ', '重症者数 (3月末, 東京基準) = '};
figname = {'Iso-ICU curve', '等重症者数曲線'}; 

for l = 1:2
    f=figure('Name',[figname{l}, '_', language{l}]);
    set(gcf, 'Position', [50*l, 50*l, 1200, 800])
    for ii = 1:nT
%         plot(yvec, ind_target(:,ii),'LineWidth',2.0, 'Color','r', 'DisplayName', [z_name{l} , num2str(target_ICU_vec(ii))])
        plot(yvec, ind_target(:,ii),'LineWidth',2.0, 'DisplayName', [z_name{l} , num2str(target_ICU_vec(ii))])
        hold on
    end
    lgd = legend;
    lgd.Location = lgdloc;
    lgd.FontSize = 14;
    ax = gca;
    ax.YAxis.FontSize = 14;
    ax.XAxis.FontSize = 14;
    ax.YAxis.Exponent = 0;
    xlim([xmin_ISO xmax_ISO])
    xlabel(x_name{l})
    yticks(11:10:101)
    yticklabels(zvec(yticks))
    ylim([1 find(abs(zvec - 1) < 0.001)])
    ylabel(y_name{l})
    title([title_name{l}, num2str(xvec), ')'],'FontSize',20,'FontWeight','normal','FontName',fn)
    lng = language{l};
    if figure_save == 1
        saveas(f, [home 'Figures/' char(pref) '/' char(save_name) '_' num2str(xvec) '_' char(lng) '.png']);
    end
end


% SimI = SimData_endogenous(:,2,1,:,:,1) + SimData_endogenous(:,2,1,:,:,2);
% % z_mat = round(squeeze(SimI(end,1,1,:,:)),0);
% z_mat = round(squeeze(SimData_endogenous(end,2,1,:,:,2)),0);
% columntable = [x_mat; z_mat'];
% 
% TI_omicron = table(  [    ...
%                 [nan(1,1); y_mat'], columntable  ...
%                 ]  );
% 
% z_mat = round(squeeze(SimData_endogenous(end,2,1,:,:,1)),0);
% columntable = [x_mat; z_mat'];
% 
% TI_delta = table(  [    ...
%                 [nan(1,1); y_mat'], columntable  ...
%                 ]  );
            