%====================== Program parameter values ======================%
fs = 16;            % common font size for many figures
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
linecolor = {'red','blue','black','green'};
language = {'EN','JP'};
th_wave = 1;


%================== Model Fixed Parameter Values ============================%
% Declare Parameters
SimPeriod = 52;        % simulation period in weeks
hconstant = 1;         % 0 = without intercept, 1 = with intercept for h regression
medical_start_date = datetime(2021,3,18);
elderly_start_date = datetime(2021,5,13);
RetroPeriod = 17;      % retroactive periods used to estimate gamma and delta
tt = 12; % Showing previous t periods for the plot
% Parameters for beta
beta_rho = 0.85;
retro_ub = 17; % Control the moving average of beta (beta_avg = sum_{t = lb}^{ub} (1/(ub-lb + 1) sum_{x=1}^t (1/t) beta_t)
retro_lb = 17;
% Parameters for mobility estimation
retroH_switch = 1; %If retroH_switch == 1, retroH = TdataGDP - 4, else = retroH
RetroH = 15;
% Population
POP_jp = 125710000;
medical_jp = 4700000;
elderly_jp = 36000000;
ordinary_jp = (POP_jp-elderly_jp-medical_jp);
accept_share = 0.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = 7/12;          % recovery rate from Covid % Should change this to 7/12 (4/25 Kohei Machi)
k = 2;                 % exponent of (1-h*alpha)
% vaccine pace
PF = 1; % 0 for AZ, 1 for PF
if PF == 0  % AZ
    E1 = 0.615;
    E2 = 0.64;
    D1 = 0.8;
    D2 = 0.85;
else   % PF
    E1 = 0.625;
    E2 = 0.895;
    D1 = 0.8;
    D2 = 0.94;
end
DRi = 6; % 経済回復速度
% Parameters for Vaccine Path
paces_ori = 3600000;
gradual_paces = 10;
sw_vacpath = 0;
% Parameters for variants
var_infection = 0.5;  % Baseline Relative Infection Rate of Variants
var_ss = 1.0;         % steady-state share of variant
var_growth = 0.47;    % growth of share of variant

%================== Parameter Values (Prefecture Specific) ============================%
% 緊急事態宣言の発令基準
th_on_vector = [750,500,400,350,1000,350,350]; % present.
th_on_vector2 = [2000,1000,800,700,2000,350,350]; % 高齢者がうち終わったあとの基準
% 緊急事態宣言の解除基準
th_off_vector = [500,100,110,80,825,50,60]; %1回目の緊急事態宣言の解除基準
th_off_vector2 = [400,100,110,80,500,50,60]; %2回目の緊急事態宣言の解除基準
th_off_vector3 = [400,100,110,80,500,50,60]; %3回目の緊急事態宣言の解除基準
% 緊急事態宣言の強さの基準: Simulation開始時のERNをいくつにするか
ERN_on_vector =  [0.55,0.5,0.5,0.5,0.65]; % 緊急事態宣言下のERN
% Size of AR(1) shock for beta process
betaT_temp_ini_vec = [0.2,0,0,0,0.15,0,0,0];
% Initial Share of Variants (Need to change these values every week)
var_initial_vector = [0.2531, 0.2909, 0.1335, 0.0909, 0.7721, 0, 0, 0,0]; % 4/26
share_index = 2; % = 2 for Monday , = 1 after Wednesday
% extrapolate var_initial_vector %
var_initial_vector = ini2now_infection_rate(var_initial_vector,var_growth,SimPeriod,share_index);
%=====================================================%

% ------------------------- Experiment Vectors ----------------------- %
% different threshold for declaring the state of emergency
THON_vector = {[900:300:1500,2000:1000:4000],400:100:800,200:100:600,300:50:500,300:100:1000,200:50:400,200:50:400,300:50:500};
THON_index_vector = {[1500,3000],[600,800],[300,600],[350,500],[400,600],[400,600],[300,600],[350,500]};

% Duration of Recovery
DR_vector = {4:2:16,0:2:16,0:2:16,0:2:16,6};
DR_index = [6,12]; %色付けして強調するものを決める

% Different ERN during the state of emergency
ERNON_vector = {0.3:0.02:0.7,[],[],[],0.3:0.02:0.7};

% Variant Infection
var_infection_vec = [0.3, 0.5];

% Vaccine Paces
paces_vector = [3600000,7000000];

% Different threshold for lifting the state of emergency
TL_vector = {100:50:500,70:10:100,80:10:130,50:10:100,100:50:500,20:10:100,10:10:100,50:10:100}; % 解除基準分析をコントロールしている Cell array
TL_index_vector = {[250,500,100],[100,80],[130,110],[100,80],[250,500,100],[50,20],[60,10]}; %色付けして強調するものを決める
