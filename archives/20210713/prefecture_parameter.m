PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka','Aichi','Fukuoka','Hyogo'};
GDPVector = [106,36,23,21,40,41,20,20]; % 兆円, one trillion yen (chou-yen)
% 緊急事態宣言の発令基準
th_on_vector = [1250,500,400,350,1000,350,350]; % present.
th_on_vector2 = [1500,1000,800,700,1500,350,350]; % 高齢者がうち終わったあとの基準

% 緊急事態宣言の解除基準
th_off_vector = [800,100,110,80,100,50,60]; %1回目の緊急事態宣言の解除基準
th_off_vector2 = [1000,100,80,60,100,50,60]; %2回目の緊急事態宣言の解除基準
th_off_vector3 = [1000,100,80,60,100,50,60]; %3回目の緊急事態宣言の解除基準

% 緊急事態宣言の強さの基準: Simulation開始時のERNをいくつにするか
alpha_on_vector = [0.2,0.1,0.1,0.1,0.2];
% Size of AR(1) shock for beta process
betaT_temp_ini_vec = [0,0,0,0,0,0,0,0];%[-0.15,0,0,0,0,0,0,0];%[0.1,0,0,0,0,0,0,0];%-0.1
start_beta_vec = [3,1,1,1,3];
beta_rho_vec = [0.5,0,0,0,0];%[0.95,0,0,0,0];
% Initial Share of Variants (Need to change these values every week)
% var_initial_vector = [0.2531, 0.2909, 0.1335, 0.0909, 0.7721, 0, 0, 0,0]; % 4/26
% var_initial_vector = [0.4035, 0.3280, 0.3261, 0.3008, 0.8224, 0, 0, 0,0]; % 5/3 (4/11-4/18のデータ)
% var_initial_vector = [0.556244, 0.450617, 0.514451, 0.422360, 0.824487, 0, 0, 0,0]; % 5/10 (4/19-4/25のデータ)
% var_initial_vector = [0.631038, 0.596774, 0.578947, 0.56,     0.856401, 0, 0, 0, 0];% 5/17 (4/26-5/2のデータ)
% var_ss_vector = [1, 1, 1, 1, 0.85, 0, 0, 0,0]; % 5/10
var_ss_vector = [1, 1, 1, 1, 0.90, 0, 0, 0,0]; % 5/17
share_index = 2; % = 2 for Monday, = 1 after Wednesday

% Parameters for Vaccine Path
paces_ori_vec = [7000000,0,0,0,7000000];
total_paces=[9800000,0,0,0,9800000];
gradual_paces_vec = [1,0,0,0,1];
sw_vacpath_vec = [0,0,0,0,0];
VT3share_vec = [0.6,0,0,0,0.6];
lag_vec = [3,0,0,0,3];
medical_duration_vec = [1,0,0,0,1];
ind_date_vec = [datetime(2021,7,22),datetime(2021,7,15),datetime(2021,7,15),datetime(2021,7,15),datetime(2021,7,15)];
paces2_ori_vec = total_paces-paces_ori_vec;
gradual_paces2_vec = [1,0,0,0,2];
paces3_ori_vec = [6300000,0,0,0,6300000]; %[4200000,0,0,0,4200000];
lagged_VT2share = 0.9;
