SimPeriod = 52;        % simulation period in weeks
gamma = 7/12;          % recovery rate from Covid % Should change this to 7/12 (4/25 Kohei Machi)
k = 2;                 % exponent of (1-h*alpha)
hconstant = 1;         % 0 = without intercept, 1 = with intercept for h regression
medical_start_date = datetime(2021,3,18);
elderly_start_date = datetime(2021,5,13);
RetroPeriod = 17;      % retroactive periods used to estimate gamma 
RetroPeriodDelta = 17;      % retroactive periods used to estimate delta
tt = 12; % Showing previous t periods for the plot
% 経済回復速度
DRi = 17; %10;
% Parameters for beta
retro_ub = 17; % Control the moving average of beta (beta_avg = sum_{t = lb}^{ub} (1/(ub-lb + 1) sum_{x=1}^t (1/t) beta_t)
retro_lb = 17;
% Parameters for mobility estimation
retroH_switch = 1; %If retroH_switch == 1, retroH = TdataGDP - 4, else = retroH
RetroH = 15;
% Parameters for variants
var_infection = 0.3;
var_infection_delta = 0.4;
var_growth = 0.47;
% Population
POP_jp = 125710000;
medical_jp = 4700000;
elderly_jp = 36000000;
ordinary_jp = (POP_jp-elderly_jp-medical_jp);
accept_share = 0.8;
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

% parameters for ICU
gamma_ICU = 7/28; % Recovery rate from ICU
ICU_limit_vec = [373,0,0,0,659]; % Define ICU limit for each prefecture
ICU_adjustment = 0.85; %0.8

% Indian Variant Parameters
var_initial2 = 0.14;
var_growth2 = 0.5;     % weekly growth parameter for logit model
var_infection2 = 0.3;
var_start = 1;         % time when the Indian variant starts spreading
