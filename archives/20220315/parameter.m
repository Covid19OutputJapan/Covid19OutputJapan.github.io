SimPeriod                   = 520;  %156;% 208;   %156      % simulation period in weeks
DR                          = 26;   %Duration of Economic Revovery after the end of the State of Emergency (weeks)
k                           = 2;    % exponent of (1-h*alpha)
hconstant                   = 1;    % 0 = without intercept, 1 = with intercept for h regression
gamma                       = 7/12; % recovery rate from Covid past_omicron_share
gamma_ICU_nation            = 7/12; %7/21; %7/26 % Recovery rate from ICU
gamma_ICU_pref              = 7/12; %7/21; % 7/28 Recovery rate from ICU
gamma_newICU_pref           = 7/12; 
gamma_Hospital              = 7/10; %7/10;
seasonal_effect             = 0.1;  % degree of seasonality (number of new cases is greater in winter & lower in summer)

ICU_nation_adjustment       = 1; %0.8
ICU_pref_adjustment         = 1; %0,8
newICU_pref_adjustment      = 1;
Hospital_adjustment         = 1; % 0.5;

% Population
POP_jp              = 125710000;
medical_jp          = 4700000;
elderly_jp          = 36000000;
ordinary_jp         = (POP_jp - elderly_jp - medical_jp);
medical_tokyo       = 5500000;
elderly_tokyo       = 3100000; % 65 and after: 3,122,068; Tokyo (2020Jan)
working_tokyo       = 9100000; % 15 ~ 64 years old: 9,109,812; Tokyo (2020Jan)
children_tokyo      = 1600000; % 0 ~ 15 years old: 1,603,044; Tokyo (2020Jan)

RetroPeriod             = 17; % retroactive periods used to estimate gamma (gamma???? maybe beta?)
RetroPeriodDelta        = 17; %8; %2;      % retroactive periods used to estimate delta (Death rate)
RetroPeriodICU_pref     = 17; %2; % retroactive periods used to estimate delta_ICU_pref (Tokyo standard Severe case)
RetroPeriodICU_nation   = 17; %2;      % retroactive periods used to estimate delta_ICU_nation (National standard Severe case)
RetroPeriodHospital     = 17; %13; %2;      % retroactive periods used to estimate delta_Hospital (Hospital)
RetroH                  = 15; % Duration periods of estimation of h (in weeks)
retroH_switch           = 1; %If retroH_switch == 1, retroH = TdataGDP - 4, else = retroH
retro_lb                = 17;
retro_ub                = 17;

accept_share            = 0.95; % How many the 2nd vaccine-dose takers recieve 3rd dose among eldelry and medical personnels
accept_share_ordinary   = 0.85; % How many the 2nd vaccine-dose takers recieve 3rd dose among eldelry and medical personnels

% vaccine effectgs
PF = 1; % 0 for AZ, 1 for PF
if PF == 0  % AZ
    E1 = 0.365; E2 = 0.625; D1 = 0.675; D2 = 0.905;     % Old data
else   % PF
    %(based on SPI_M_O, October 13th https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1027851/S1383_SPI-M-O_Summary_autumn_winter_scenarios.pdf)
    E1 = (0.33 + 0.62)/2 ;          % = 0.4750 Mean value of ICL and LSHTM assumpitons for the reducito in risk of infection
    E2 = 0.5;  %感染予防効果減退を考慮 ... (0.85 + 0.80 + 0.85)/3  = 0.8333 Mean value of ICL, LSHTM, and Warwick assumpitons for the reducito in risk of infection
    E3 = 0.85; % ... (0.85 + 0.80 + 0.85)/3  = 0.8333 Mean value of ICL, LSHTM, and Warwick assumpitons for the reducito in risk of infection
    D1 = (0.85 + 0.92)/2 ;          % = 0.8850  Mean value of ICL and LSHTM assumpitons for reduciton in risk of death
    D2 = (0.95 + 0.96 + 0.98)/3 ;   % =0.9633 Mean value of ICL and LSHTM assumpitons for reduciton in risk of death
    D3 = D2 ;   % =0.9633 Mean value of ICL and LSHTM assumpitons for reduciton in risk of death
end

% AR1 shock parameters 
betaT_temp_ini                  = 0.0; %0.55 
start_beta                      = 1;
beta_rho                        = 0.75; 

delta_temp_ini                  = -0.75;
start_delta                     = 1;
delta_rho                       = 0.95;

delta_ICU_nation_temp_ini       = 0.8;
start_delta_ICU_nation          = 1;
delta_ICU_nation_rho            = 0.75;

delta_ICU_pref_temp_ini         =  0;
start_delta_ICU_pref            = 1;
delta_ICU_pref_rho              = 0.95;

delta_newICU_pref_temp_ini      = 0;
start_delta_newICU_pref         = 1;
delta_newICU_pref_rho           = 0.95;

delta_Hospital_temp_ini         = 3.0;
start_delta_Hospital            = 1;
delta_Hospital_rho              = 0.99;


% Omicron paramters
omicron_initial = 0.999;
omicron_ss      = 1;
omicron_growth  = 1.75;
iniI_omicron    = 10;
omicron_immunity = 0.2;%1;

