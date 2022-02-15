clear variables
close all
iPC = 0; % 0 for Mac, 1 for Windows

if iPC == 1
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
    fn = 'Yu Gothic'; % Font style for xaxis, yaxis, title
else
%         home = '/Users/sohtakawawaki/Dropbox/fujii_nakata (1)/Website/Codes/';
    %     home = '/Users/okamotowataru/Dropbox/fujii_nakata/Website/Codes/';
    home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Covid19OutputJapan.github.io/archives/20211223/Codes_omicron/';
%     home = '/Users/shotaro/Dropbox/fujii_nakata/Website/Codes/';
    fn = 'YuGothic';
end

cd(home);
%====================== Program parameter values ======================%
pref = 'Tokyo';
prefGDP = 106;
figure_save = 0; % 0 = figures won't be saved, 1 = they will be saved
data_save = 0; % save back data

exogenous_omicron_share = 0;
no_omicron = 0;     %
alternative_scenario = 1; %=0 ... x loop is for relative BRN, =1 ... x loop is for relative VE
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
if no_omicron == 1
    iniI_omicron = 0;
end

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

% xmin = find(date == datetime(2021, 11, 04));
% xmax = find(date == datetime(2022, 12, 01));
xmax = find(date == datetime(2022, 3, 31));
xmin = find(date == datetime(2021, 7, 1));
ymax_D = 50;

% Simulation parameters
omicron_relative_infectivity_vector = [1, 1.25, 1.5];
omicron_E2_vector = [1, 0.6, 0.2]; %Relative
omicron_realtive_severity_vector = [1, 0.5, 0.25];
if no_omicron == 1
    omicron_relative_infectivity_vector = [1,1,1,1];
    omicron_E2_vector = [1,1,1,1];
    omicron_realtive_severity_vector = [1,1];
end
% xvec = 1;
% yvec = 1;
% zvec = 1;

% Directory for saved figures
if alternative_scenario == 0
    xvec = omicron_relative_infectivity_vector;
    yvec = omicron_E2_vector;
    zvec = omicron_realtive_severity_vector;
    nX = length(xvec);
    nY = length(yvec);
    nZ = length(zvec);
    
    
    figfolder = string(['Omicron_Analysis_Relative_Infectivity']);
    % Scenario Name and Line
    titlevec = {'relative BRN = ', '相対的基本再生産数　＝　'};
    figname_xvar = '_relativeBRN_';
    Scenario = ["相対的基本再生産数　1倍"; "相対的基本再生産数　1.25倍"; "相対的基本再生産数　1.5倍"];
    ScenarioEN = [ "Relative BRN　1倍"; "Relative BRN　1.25倍"; "Relative BRN　1.5倍"];
    Scenario_vec = [ScenarioEN,Scenario];
    linecolor = {"b", "b", "b", "b"; ...
                 "k", "k", "k", "k"; ...
                 "r", "r", "r", "r"};
    LineStyles = {"-", "--", "-.", ":"; ...
                  "-", "--", "-.", ":"; ...
                  "-", "--", "-.", ":"};
    lineWidth = [1.5, 1.5, 1.5, 1.5; 1.5, 1.5, 1.5, 1.5; 1.5, 1.5, 1.5, 1.5];
    markertype = {'o','s','d','p';'o','s','d','p';'o','s','d','p'};
    lineNameJP = strings(nY,nZ);
    lineNameEN = strings(nY,nZ);

    for iY = 1:nY
        for iZ = 1:nZ
            lineNameJP(iY,iZ) = ['相対的VE = ', num2str(yvec(iY)), ', 相対的重症化率 = ', num2str(zvec(iZ))];
            lineNameEN(iY,iZ) = ['Rel. VE = ', num2str(yvec(iY)), ', Rel. Sev. = ', num2str(zvec(iZ))];
        end
    end

elseif alternative_scenario == 1
    xvec = omicron_E2_vector;
    yvec = omicron_relative_infectivity_vector;
    zvec = omicron_realtive_severity_vector;
    nX = length(xvec);
    nY = length(yvec);
    nZ = length(zvec);
    
    figfolder = string(['Omicron_Analysis_Relative_VE']);
    titlevec = {'relative VE = ', '相対的ワクチン有効性　＝　'};
    figname_xvar = '_relativeVE_';
    Scenario = ["相対的ワクチン有効性 1.0倍"; "相対的ワクチン有効性 0.6倍"; "相対的ワクチン有効性 0.2倍"];
    ScenarioEN = ["Relative Vaccine Effectiveness 1.0"; ...
                  "Relative Vaccine Effectiveness 0.6"; "Relative Vaccine Effectiveness 0.2"];
    Scenario_vec = [ScenarioEN,Scenario];
    linecolor = {"b", "b", "b"; "k", "k", "k"; "r", "r", "m"; "m", "m", "m"};
    LineStyles = {"-", "--", ":"; "-", "--", ":"; "-", "--", ":"; "-", "--", ":"};
    lineWidth = [1.5, 1.5, 1.5; 1.5, 1.5, 1.5; 1.5, 1.5, 1.5; 1.5, 1.5, 1.5];
    markertype = {'o','s','d','p';'o','s','d','p';'o','s','d','p'};
    lineNameJP = strings(nY,nZ);
    lineNameEN = strings(nY,nZ);
    for iY = 1:nY
        for iZ = 1:nZ
            lineNameJP(iY,iZ) = ['相対的BRN = ', num2str(yvec(iY)), ', 相対的重症化率 = ', num2str(zvec(iZ))];
            lineNameEN(iY,iZ) = ['Rel. BRN = ', num2str(yvec(iY)), ', Rel. Sev. = ', num2str(zvec(iZ))];
        end
    end
end

if no_omicron == 1
    figfolder = string(['No_Omicron']);
%    linecolor = {'b','k','r','m'};
    beta_goal_vec = [1.7580, 2.1984, 2.9300, 3.5175]; % [1.7580, 2.0500,2.1984, 2.3430, 2.9300, 3.5175] ... 3, 3.5, 3.75, 4, 5, 6
    yvec = beta_goal_vec;
    nX = 1;
    nY = length(yvec);
    nZ = 1;
%     linecolor = {'k'};
    linecolor = {"b","k", "r", "m"}';
    LineStyles = {"-", "-", "-", "-"}';
    lineWidth = [1.5, 1.5, 1.5, 1.5;1.5, 1.5, 1.5, 1.5]';
    markertype = {'o','s','d','p';'o','s','d','p';'o','s','d','p'};
    yname = [3, 3.75, 5, 6];
    lineNameJP = strings(nY,nZ);
    lineNameEN = strings(nY,nZ);
    for iY = 1:nY
        for iZ = 1:nZ
            lineNameJP(iY,iZ) = ['BRN = ', num2str(yname(iY))];
            lineNameEN(iY,iZ) = ['BRN = ', num2str(yname(iY))];
        end
    end
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
Sim_dD                    = AlphaPath;
SimERN                  = AlphaPath;
SimBRN                 = AlphaPath;
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

for iX = 1:nX    

    for iY = 1:nY

        if alternative_scenario == 0
            omicron_relative_infectivity    = omicron_relative_infectivity_vector(iX);
            relative_infectivity_path       = ones(SimPeriod,1); %vector
            omicronE3 = originalE3;
            omicronE2 = originalE2 * omicron_E2_vector(iY);
            omicronE1 = originalE1 * omicron_E2_vector(iY);
        elseif alternative_scenario == 1
            omicronE3 = originalE3;
            omicronE2 = originalE2 * omicron_E2_vector(iX);
            omicronE1 = originalE1 * omicron_E2_vector(iX);
            omicron_relative_infectivity    = omicron_relative_infectivity_vector(iY);
            relative_infectivity_path       = ones(SimPeriod,1); %vector
        end
        if no_omicron == 1
            ori_beta_goal  = beta_goal_vec(iY);
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
            
            betaT               = beta_AR1(betaT_temp_ini, beta_rho, betaT, start_beta);  %2021/12/23 kawawaki
            
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
                    Covid_projection_omicron(InitialValues,alpha_on,alpha_off,th_on,th_off,...
                    betaT,gammaT,deltaT,delta_ICU_nation,delta_ICU_pref,delta_Hospital,h,k,POP0,...
                    lagged_cumsumVT1, lagged_cumsumVT2, E1, E2, cum_in_R, ...
                    hconstant,DRi,ICU_nation_inflow_avg, ICU_pref_inflow_avg, Hospital_inflow_avg, ...
                    gamma_ICU_nation, gamma_ICU_pref, gamma_Hospital,...
                    relative_beta_SOE, beta_jump, beta_goal, seasonality, alphaBox(Tdata + 1:end), state, simple_beta_avg);
                Sim_dD(:, iX, iY, iZ) = squeeze(SimData(2:end,4,iX,iY,iZ)-SimData(1:end-1,4,iX,iY,iZ)+SimData(2:end,4,iX,iY,iZ)-SimData(1:end-1,4,iX,iY,iZ));
                
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
                    = Covid_projection_endogenous_omicron(SimData_endogenous(:, :, iX,iY,iZ,:),...
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
                
                Sim_dD(:, iX, iY, iZ) = squeeze(SimData_endogenous(2:end,4,iX,iY,iZ,1)-SimData_endogenous(1:end-1,4,iX,iY,iZ,1)+SimData_endogenous(2:end,4,iX,iY,iZ,2)-SimData_endogenous(1:end-1,4,iX,iY,iZ,2));
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
%%

% %==================== Plot and Backdata ====================%
yft = '%.2f';
lgdfs = 12;
axfs = 12;
column_num = 3;
lgdLocation = 'northwest';
l = 2;
if l == 1
    lineName = lineNameEN;
elseif l == 2
    lineName = lineNameJP;
end
% Plot Omicron Share
% if exogenous_omicron_share == 0
%     for iX = 1:nX
%         fig=figure('Name', [char(figname_omicron_share) ' ' num2str(iX) ]   );
%         set(gcf, 'Position', [100, 100, 1200, 800])
%         title_vec = ["Transition of the Share of Omicron Variant, ", "オミクロン株割合の推移, "];
%         titlevec = {'relative infectivity = ', '相対感染力　＝　'};
%         plot_4Dfunction(past_omicron_share, omicron_share_mat(:,:,:,:), iX, ...
%             WeekNumber, YearMonth, xmin, xmax, ...
%             fn, fs, lgdfs, axfs,yft,...
%             lgdLocation, column_num, l, title_vec, ...
%             lineWidth,linecolor, LineStyles, lineName)
%         title([title_vec{l},titlevec{l}, num2str(xvec(iX),'%.2f')],'FontSize',20)
%         ylim([0 1])
%         ax = gca;
%         ax.YAxis.FontSize = axfs;   ax.XAxis.FontSize = axfs;   ax.YAxis.Exponent = 0;
%         xtickangle(45)
%         if figure_save == 1
%             saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' ...
%                 [char(figname_omicron_share), '_rel_inf_' num2str(xvec(iX)) ] '.png']);
%         end
%     end
% end
% Save variable names and values for backdata files
minAlpha = alpha_off; % 経済損失0の基準
AlphaM = AlphaMat;
% AlphaM = AlphaMat(~isnan(AlphaMat));
AlphaM = (AlphaM - minAlpha) * prefGDP * 10000 * SimPeriod / 52;
DM = DMat;
% DM = DMat(~isnan(DMat));
% BackDataDA(1:nX*nY*nZ) = [round(AlphaM') / 10000, round(DM'), Scenario'];

%--- Record how many times on and off are triggered ---%
waves = zeros(nX,nY,nZ);
for ixx = 1:nX
    for iyy = 1:nY
        for izz = 1:nZ
            svec = zeros(SimPeriod - 1, 1);
            
            for t = 1:SimPeriod - 1
                svec(t) = AlphaPath(t + 1, ixx,iyy,izz) - AlphaPath(t, ixx,iyy,izz);
            end
            
            waves(ixx,iyy,izz) = sum(svec > 0);
        end
    end
end

for l = 1:2 %1:2 when english version needed
    
    for iX = 1:nX
        if l == 1
            lineName = lineNameEN;
        elseif l == 2
            lineName = lineNameJP;
        end
        
        % Generate graphs for the website
        lng = language{l};
        figname = [char(figname_main)  char(figname_xvar) num2str(xvec(iX),'%.2f') '_' char(lng)];
        if no_omicron == 1
            figname = [char(figname_main) '_' char(lng)];
        end
        f = figure('Name', figname);
        t = tiledlayout(3,3, 'TileSpacing', 'compact');
        if no_omicron == 0
            title(t,[titlevec{l}, num2str(xvec(iX),'%.2f')],'FontSize',20)
        end
        f.WindowState = 'maximized';
        lgdLocation = 'NorthWest';
        yft = '%.0f';
        lgdfs = 6;
        axfs = 8;
        column_num = 1;
        
        %--- Number of people who are in ICU (local standard)---%
        %         subplot(3, 3, 1)
        nexttile
        
        title_vec = ["ICU (Tokyo Standard)", "重症患者数（都基準）"];
        plot_4Dfunction(ICU_pref(2:end), SimICU_pref(2:end,:,:,:), iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            fn, fs, lgdfs, axfs,yft,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
%         Handle = legend;
%         set(Handle, 'Visible', 'off');
        hold on
        plot([ICU_limit_pref_vec; ones(SimPeriod, 1) * ICU_limit_pref_vec(end)], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        text(xmin, ICU_limit_pref_vec(xmin) * 0.85, '100%', 'FontSize', fs)
        hold on
        plot([ICU_limit_pref_vec * 0.5; ones(SimPeriod, 1) * ICU_limit_pref_vec(end) * 0.5], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        hold on
%         plot([nan(Tdata,1); ones(SimPeriod, 1) * th_on], '-r', 'HandleVisibility', 'off', 'LineWidth', 1.0)
%         hold on
%         plot([nan(Tdata,1); ones(SimPeriod, 1) * th_off], '-b', 'HandleVisibility', 'off', 'LineWidth', 1.0)
        text(xmin, ICU_limit_pref_vec(xmin) * 0.4, '50%', 'FontSize', fs)
        ylim([0 max(ICU_limit_pref_vec)*1.5])
        
        %--- Number of people who are in ICU (Naitonal standard) ---%
        %         subplot(3, 3, 2)
        nexttile
        title_vec = ["ICU (National Standard)", "重症患者数（国基準）"];
        plot_4Dfunction(ICU_nation(2:end), SimICU_nation(2:end,:,:,:), iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            fn, fs, lgdfs, axfs,yft,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        hold on
        plot([BED; ones(SimPeriod, 1) * BED(end)], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        text(xmin, BED(xmin) * 0.85, '100%', 'FontSize', fs)
        hold on
        plot([BED * 0.5; ones(SimPeriod, 1) * BED(end) * 0.5], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        text(xmin, BED(xmin) * 0.4, '50%', 'FontSize', fs)
        ylim([0 max(BED)*1.5])
        
        %--- Number of newly hospitalized ---%
        %         subplot(3, 3, 3)
        nexttile
        title_vec = ["Hospitalized Patients", "入院患者数"];
        plot_4Dfunction(hospital(2:end), SimHospital(2:end,:,:,:), iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            fn, fs, lgdfs, axfs,yft,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        hold on
        plot([Hospital_limit_vec; ones(SimPeriod, 1) * Hospital_limit_vec(end)], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        text(xmin, Hospital_limit_vec(xmin) * 0.85, '100%', 'FontSize', fs)
        hold on
        plot([Hospital_limit_vec * 0.5; ones(SimPeriod, 1) * Hospital_limit_vec(end) * 0.5], '--k', 'HandleVisibility', 'off', 'LineWidth', 1.5)
        text(xmin, Hospital_limit_vec(xmin) * 0.4, '50%', 'FontSize', fs)
        ylim([0 max(Hospital_limit_vec)*1.5])
        
        %         subplot(3, 3, 4)
        nexttile
        title_vec = ["New Deaths (Daily Average)", "新規死亡者数（1日平均）"];
        plot_4Dfunction(dD / 7, Sim_dD/ 7, iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            fn, fs, lgdfs, axfs,yft,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        ylim([0 ymax_D])
        %
        %--- Number of new cases ---%
        %         subplot(3, 3, 5)
        nexttile
        title_vec = ["New Cases (Daily Average)", "新規感染者数（1日平均）"];
        plot_4Dfunction(N/7, NPath/7, iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            fn, fs, lgdfs, axfs,yft,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        ylim([0,10000])
        %--- GDP Path ---%
        yft = '%.2f';
        %         subplot(3, 3, 6)
        nexttile
        title_vec = ["GDP", "GDP"];
        plot_4Dfunction(100*(1-alpha), 100*(1-AlphaPath), iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            fn, fs, lgdfs, axfs,yft,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        ylim([90 100])
        Handle = legend;
        set(Handle, 'Visible', 'off');
        
        %         subplot(3, 3, 7)
        nexttile
        title_vec = ["Transitions of ERN", "実効再生産数"];
        plot_4Dfunction(ERN, SimERN, iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            fn, fs, lgdfs, axfs,yft,...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        
        
        %         subplot(3, 3, 8) %Cumulative vaccinated
        nexttile
        
        VT_23 = VT;
        VT_23(:,[1,4,7]) = [];
        plot_vaccinepath_separate_percentage(2, VT_23, V2_medical, V3_medical, V2_elderly, V3_elderly, V2_others, V3_others, ps, YearMonthWeek(:, l), WeekNumber, Tdata, fs, 8, fn, xmin, xmax, l, POP0);
        yticks(0:10:100)
        grid on
        xlim([xmin, xmax])
        %     xticks(find(WeekNumber == 1 & abs(1 - mod(MonthNumber, 4)) < 0.01))
        xticks(find(WeekNumber == 1))
        %     xticklabels(YearMonthWeekJP(xticks))
        xticklabels(YearMonthJP(xticks))
        lgd = legend;
        lgd.FontSize = 12;
        lgd.Location = 'Southeast';
        ax = gca;
        ax.YAxis.FontSize = axfs;
        ax.XAxis.FontSize = axfs;
        %     xline(xline_ind, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        %--- Trade-off Curve ---%
        %         subplot(3, 3, 9)
        nexttile
        title_vec = ["Transition of the Share of Omicron Variant, ", "オミクロン株割合の推移, "];
        plot_4Dfunction(past_omicron_share, omicron_N_share_mat(:,:,:,:), iX, ...
            WeekNumber, YearMonth, xmin, xmax, ...
            fn, fs, lgdfs, axfs,'%.1f',...
            lgdLocation, column_num, l, title_vec, ...
            lineWidth,linecolor, LineStyles, lineName)
        Handle = legend;
        set(Handle, 'Visible', 'off');
        ylim([0 1])
        grid on
        ax = gca;
        box on
        xtickangle(45)

        
%         for iyy = 1:nY
%             for izz = 1:nZ
%                 scatter(AlphaM(iX,iyy,izz),DM(iX,iyy,izz),250,...
%                     markertype{iyy,izz}, linecolor{iyy,izz},'filled');
%                 hold on
%             end
%         end
%         if l == 1
%             xlabel('Output Loss (hundred million yen)','FontSize',fs)
%             ylabel('Cumulative Deaths','FontSize',fs)
%             title('Relationship between Covid-19 and Output(within 10 years)','FontSize',fs,'FontWeight','normal')
%         elseif l == 2
%             xlabel('経済損失 (兆円)','FontSize',fs,'FontName',fn)
%             ylabel('累計死亡者数','FontSize',fs,'FontName',fn)
%             title('コロナ感染と経済の関係(今後10年)','FontSize',fs,'FontWeight','normal','FontName',fn)
%         end
%         xlim([0,inf])
%         xtickangle(45)
%         grid on
%         ax = gca;
%         ax.YAxis.FontSize = 12;
%         ax.XAxis.FontSize = 12;
%         ax.YAxis.Exponent = 0;
%         ax.XAxis.Exponent = 0;
%         ytickformat('%,6.0f')
%         box on
%         grid on
%         
        lgd = legend(nexttile(1));
        lgd.Layout.Tile  = 'south';
        lgd.NumColumns = column_num_main;
        lgd.FontSize = ldfs_main;
        %     xticks(0:2.5:30)
        %     ylim([max(0, mean(DM(iX,:,:),[2,3]) - 7500) mean(DM(iX,:,:),[2,3]) + 7500])
        
        if figure_save == 1
            %         saveas(f, [home 'Figures/' char(pref) '/' char(figname_main) char(lng) '.png']);
            saveas(f, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname) '.png']);
        end
        
        
    end
end %End of language loop = figure loop

%% Other Plot
lgdfs = 12;
axfs = 12;
lineWidth = 2.0*ones(nY,nZ);
omicron_plot_parameter


%% Backdata
BRNpast = beta_tilde ./ (gamma + delta);
BackDataN           = zeros(8 + SimPeriod, nX, nY, nZ);
BackDataAlpha       = BackDataN;
BackDataERN         = BackDataN;
BackDataBRN         = BackDataN;
BackDatadD          = BackDataN;
BackDataICU_pref    = BackDataN;
BackDataICU_nation  = BackDataN;
BackDataHospital    = BackDataN;
for iX = 1:nX
    for iY = 1:nY
        for iZ = 1:nZ
            BackDataN(:, iX, iY, iZ)            = [N(end - 7:end);          NPath(:, iX, iY, iZ)];
            BackDataAlpha(:, iX, iY, iZ)        = [alpha(end - 7:end);      AlphaPath(:, iX, iY, iZ)];
            BackDataERN(:, iX, iY, iZ)          = [ERN(end - 7:end);        SimERN(:, iX, iY, iZ)];
            BackDataBRN(:, iX, iY, iZ)          = [BRNpast(end - 7:end);    SimBRN(:, iX, iY, iZ)];
            BackDatadD(:, iX, iY, iZ)           = [dD(end - 7:end);         Sim_dD(:, iX, iY, iZ)];
            BackDataICU_pref(:, iX, iY, iZ)     = [ICU_pref(end - 7:end);   SimICU_pref(2:end, iX, iY, iZ)];
            BackDataICU_nation(:, iX, iY, iZ)   = [ICU_nation(end - 7:end); SimICU_nation(2:end, iX, iY, iZ)];
            BackDataHospital(:, iX, iY, iZ)     = [hospital(end - 7:end);   SimHospital(2:end, iX, iY, iZ)];
        end
    end
end
%%
% if data_save == 1
    titleN = strings(1, 1 + nZ * 8);
    titleN(1) = "週";
    for ti = 1:nZ
        titleN(1, 1 + ti) = append("新規感染者数（", lineNameJP(iY,ti), "）");
        titleN(1, 1 + nZ + ti) = append("経済活動（", lineNameJP(iY,ti), "）");
        titleN(1, 1 + nZ * 2 + ti) = append("実効再生産数（", lineNameJP(iY,ti), "）");
        titleN(1, 1 + nZ * 3 + ti) = append("基本再生産数（", lineNameJP(iY,ti), "）");
        titleN(1, 1 + nZ * 4 + ti) = append("入院患者数（", lineNameJP(iY,ti), "）");
        titleN(1, 1 + nZ * 5 + ti) = append("重症者数_国基準（", lineNameJP(iY,ti), "）");
        titleN(1, 1 + nZ * 6 + ti) = append("重症者数_都基準（", lineNameJP(iY,ti), "）");
        titleN(1, 1 + nZ * 7 + ti) = append("新規死亡者数（", lineNameJP(iY,ti), "）");
    end
         TN = table([
                titleN; 
                YearMonthWeekJP(Tdata - 7:end - 1), ...
                squeeze(round(BackDataN(:, iX, iY, :) / 7)), ...
                squeeze(round(100 * (1 - BackDataAlpha(:, iX, iY, :)), 1)), ...
                squeeze(round(BackDataERN(:, iX, iY, :), 2)), ...
                squeeze(round(BackDataBRN(:, iX, iY, :), 2)), ...
                squeeze(round(BackDataHospital(:, iX, iY, :))), ...
                squeeze(round(BackDataICU_nation(:, iX, iY, :))), ...
                squeeze(round(BackDataICU_pref(:, iX, iY, :))), ...
                squeeze(round(BackDatadD(:, iX, iY, :) / 7))
                ]);
%%
% end

% %%
% if data_save == 1
%     titleN = strings(1, 1 + length(TH_index) * 8);
%     titleN(1) = "週";
%
%     TN = table([titleN; YearMonthWeekJP(Tdata - 7:end - 1), ...
%                 round(BackDataN(:, 1:length(TH_index)) / 7), ...
%                 round(100 * (1 - BackDataAlpha(:, 1:length(TH_index))), 1), ...
%                 round(BackDataERN(:, 1:length(TH_index)), 2), ...
%                 round(BackDataICU_nation(:, 1:length(TH_index))), ...
%                 round(BackDataICU_pref(:, 1:length(TH_index))), ...
%                 round(BackDatadD(:, 1:length(TH_index)) / 7), ...
%                 round(BackDataHospital(:, 1:length(TH_index))),...
%                 round(BackDataBRN(:, 1:length(TH_index)), 2)]);
%
%     titleAD = ["経済損失（兆円）", "死亡者数", "ケース"];
%     TAD = table([titleAD; BackDataDA(1:length(TH), :)]);
%
%     %     TVAC = table([titleVAC; ...
%     %         MonthWeekJP(x_left-3:x_right),...
%     %         round(BackData_Area_nv_all(x_left-3:x_right,:))*10000]);
%     TVAC_new = table([titleVAC; ...
%                 YearMonthWeekJP(x_left - 3:x_right), ...
%                     round(BackData_Area_nv_all(x_left - 3:x_right, :) * 10000)]);
%
%     TVAC_cum = table([titleVAC; ...
%                     YearMonthWeekJP(x_left - 3:x_right), ...
%                     round(BackData_Area_cv_all(x_left - 3:x_right, :) * 100000000)]);
%
%     writetable(TN, [home 'Figures/' char(pref) '/' char(figfolder) '/BackData_' char(figname_main) char(pref) '.xls'], 'Sheet', '新規感染者数（1日平均）', 'WriteVariableNames', false);
%     writetable(TVAC_new, [home 'Figures/' char(pref) '/' char(figfolder) '/BackData_' char(figname_main) char(pref) '.xls'], 'Sheet', 'ワクチン新規接種パス', 'WriteVariableNames', false);
%     writetable(TVAC_cum, [home 'Figures/' char(pref) '/' char(figfolder) '/BackData_' char(figname_main) char(pref) '.xls'], 'Sheet', 'ワクチン累計接種パス', 'WriteVariableNames', false);
%     writetable(TAD, [home 'Figures/' char(pref) '/' char(figfolder) '/BackData_' char(figname_main) char(pref) '.xls'], 'Sheet', '経済損失と死亡者数', 'WriteVariableNames', false);
% end
%
