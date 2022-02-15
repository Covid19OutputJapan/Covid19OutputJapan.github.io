SimPeriod =  520;%156;% 208;   %156      % simulation period in weeks
gamma = 7/12;          % recovery rate from Covid 
k = 2;                 % exponent of (1-h*alpha)
hconstant = 1;         % 0 = without intercept, 1 = with intercept for h regression
medical_start_date = datetime(2021,3,18);
elderly_start_date = datetime(2021,5,13);
RetroPeriod = 17;      % retroactive periods used to estimate gamma (gamma???? maybe beta? 7/28 Mori)
RetroPeriodDelta = 17; %8; %2;      % retroactive periods used to estimate delta (Death rate)
RetroPeriodICU_pref = 17; %2; % retroactive periods used to estimate delta_ICU_pref (Tokyo standard Severe case)
RetroPeriodICU_nation = 17; %2;      % retroactive periods used to estimate delta_ICU_nation (National standard Severe case)
RetroPeriodHospital = 17; %13; %2;      % retroactive periods used to estimate delta_Hospital (Hospital)

lambda_delta = 0.3 * 0.95; %0.95;
lambda_delta_ICU_nation = 0.25 * 0.75; %0.5; %0.825
lambda_delta_ICU_pref = 0.375 * 0.75; %0.5; %0.825
lambda_delta_Hospital = 0.3 * 0.75; %0.5; %0.825
vac_end_date = 11;

lambda_delta =[transpose(linspace(1,lambda_delta,vac_end_date));ones(SimPeriod-vac_end_date,1)*lambda_delta];
lambda_delta_ICU_nation = [transpose(linspace(1,lambda_delta_ICU_nation,vac_end_date));ones(SimPeriod-vac_end_date,1)*lambda_delta_ICU_nation]; % 0.825;
lambda_delta_ICU_pref = [transpose(linspace(1,lambda_delta_ICU_pref,vac_end_date));ones(SimPeriod-vac_end_date,1)*lambda_delta_ICU_pref]; %0.825;
lambda_delta_Hospital = [transpose(linspace(1,lambda_delta_Hospital,vac_end_date));ones(SimPeriod-vac_end_date,1)*lambda_delta_Hospital]; %0.825;

tt = 12; % Showing previous t periods for the plot
% 経済回復速度
DRi= 26;%17;%26 %10;
% Parameters for beta
retro_ub = 17; % Control the moving average of beta (beta_avg = sum_{t = lb}^{ub} (1/(ub-lb + 1) sum_{x=1}^t (1/t) beta_t)
retro_lb = 17;
% Parameters for mobility estimation
retroH_switch = 1; %If retroH_switch == 1, retroH = TdataGDP - 4, else = retroH
RetroH = 15;
% Parameters for variants
var_infection = 0.3; %Relative increase of infectiousness (alpha varaint)
var_infection_delta = 0.4; %Relative increase of death rate (alpha varaint)
var_growth = 0.47;
% Population
POP_jp = 125710000;
medical_jp = 4700000;
medical_tokyo = 5500000;
elderly_jp = 36000000;
elderly_tokyo = 3100000;  %%% 65 and after: 3,122,068; Tokyo (2020Jan)
ordinary_jp = (POP_jp - elderly_jp - medical_jp);
working_tokyo = 9100000; %%%15 ~ 64 years old: 9,109,812; Tokyo (2020Jan)
children_tokyo = 1600000; %%%0 ~ 15 years old: 1,603,044; Tokyo (2020Jan)

%基本 or 悲観
% accept_share = 0.9; 
% accept_share_ordinary = 0.6449; %so that an accept share of age 13-64 = 75%   

%希望
% accept_share = 0.95; 
% accept_share_ordinary = 0.7739; %so that an accept share of age 13-64 = 90%

accept_share = 0.92;
%accept_share_ordinary = 0.6879; %so that an accept share of age 13-64 = 80%  
accept_share_ordinary = 0.7; %so that an accept share of age 13-64 = 80%  


accept_share_working_tokyo = 0.8;
% vaccine pace 
PF = 1; % 0 for AZ, 1 for PF
if PF == 0  % AZ
    % Old data
    E1 = 0.365; %0.615;
    E2 = 0.625; %0.64;
    D1 = 0.675; %0.8;85
    D2 = 0.905; %0.85;
else   % PF
    %(based on SPI_M_O, October 13th https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1027851/S1383_SPI-M-O_Summary_autumn_winter_scenarios.pdf)
    E1 = (0.33 + 0.62)/2 ;          % = 0.4750 Mean value of ICL and LSHTM assumpitons for the reducito in risk of infection
    E2 = 0.7; %感染予防効果減退を考慮 ... (0.85 + 0.80 + 0.85)/3  = 0.8333 Mean value of ICL, LSHTM, and Warwick assumpitons for the reducito in risk of infection
    E3 = 0.85; % ... (0.85 + 0.80 + 0.85)/3  = 0.8333 Mean value of ICL, LSHTM, and Warwick assumpitons for the reducito in risk of infection
    D1 = (0.85 + 0.92)/2 ;          % = 0.8850  Mean value of ICL and LSHTM assumpitons for reduciton in risk of death
    D2 = (0.95 + 0.96 + 0.98)/3 ;   % =0.9633 Mean value of ICL and LSHTM assumpitons for reduciton in risk of death
    D3 = D2 ;   % =0.9633 Mean value of ICL and LSHTM assumpitons for reduciton in risk of death
   
    
    %(based on SPI_M_O, July 7th)
%     E1 = 0.45;  
%     E2 = 0.75; %感染予防効果減退を考慮 %E2 = 0.815; from SPI_M_O July 7th
%     D1 = 0.865; 
%     D2 = 0.96; 
    %(based on SPI_M_O, June 9th)
    %     E1 = 0.38; E2 = 0.795; D1 = 0.68; D2 = 0.925;

end
D1 = (D1 - E1)/(1-E1);      % Find reduction of death conditoinal on infection after first does
D2 = (D2 - E2)/(1-E2);      % Find reduction of death conditoinal on infection after second does
D3 = (D3 - E3)/(1-E3);      % Find reduction of death conditoinal on infection after second does
gapThirdVac = 34+1; %8 months (from May 6th to Jan 6th)


% parameters for ICU
gamma_ICU_nation = 7/12; %7/21; %7/26 % Recovery rate from ICU
gamma_ICU_pref = 7/12; %7/21; % 7/28 Recovery rate from ICU
ICU_nation_adjustment = 1; %0.8
ICU_pref_adjustment = 1; %0,8

% parameters for Hospitalizaiton
gamma_Hospital = 7/10; %7/10;
Hospital_adjustment = 1; % 0.5;
Hospital_limit = 6406; %5882;

% Indian Variant Parameters
% var_initial2 = 0.50;
var_growth2 = 0.75;     % weekly growth parameter for logit model
var_infection2 = 0.5; %0.3;%0.2; % relative infection rate for delta variant compared to alpha variant
var_infection_delta2 = 0.1; % relative increase of death rate for delta variant compared to alpha variant
var_start = 1;         % time when the Indian variant starts spreading

seasonal_effect = 0.1;   % degree of seasonality
% seasonal_effect = 0.4;   % degree of seasonality
% seasonal_effect = 0;   % degree of seasonality

beta_goal = 1.3;

% Beta Ar1 shock parameter
betaT_temp_ini = 0.55; 
start_beta = 1;
beta_rho = 0.75; 

% Delta AR1 shock parameter
delta_temp_ini = 0;
start_delta = 1;
delta_rho = 0.5;

delta_ICU_nation_temp_ini = 0;
start_delta_ICU_nation = 1;
delta_ICU_nation_rho = 0.95;

delta_ICU_pref_temp_ini =  0;
start_delta_ICU_pref = 1;
delta_ICU_pref_rho = 0.95;

delta_Hospital_temp_ini = 0;
start_delta_Hospital = 1;
delta_Hospital_rho = 0.95;

% Omicron paramters
omicron_initial = 0.01;
omicron_ss      = 1;
omicron_growth  = 1.75;
iniI_omicron    = 10;
omicron_immunity = 0.2;%1;

