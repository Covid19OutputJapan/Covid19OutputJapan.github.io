% %--- Import Prefecture-specifc data ---%
pref = PrefVector{pindex};        % prefecture to be analyzed
prefGDP = GDPVector(pindex);

disp('Current analysis is on')
disp(pref)

% Covid data are recorded at weekly frequencies (Mon-Sun)
% The first week start on January 20 (Mon), 2020
if iPC==1
    covid = importdata([home '\Covid_weekly_newV.csv']);  % Import weekly Covid data by prefecture
else
    covid = importdata([home 'Covid_weekly_newV.csv']);  % Import weekly Covid data by prefecture
end
Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
% Columns: 1 = date, 2 = new positive, 3 = death, 4 = mobility,
% 5 = GDP, 6 = population, 7 = GDO
dateD = Data(:,1) + 21916;
N = Data(:,2);
dD = Data(:,3);
M = Data(:,4);
% GDP = Data(:,5);
POP = Data(:,6);
GDP = Data(:,7);
GDP = GDP/GDP(1)*100;   % Normalization: 100 = 2020JAN
Tdata= size(dateD,1);    % Data period in weeks

%--- Import ICU data ---%
ICU_nation = zeros(Tdata + 1, 1);
ICU_pref = zeros(Tdata + 1, 1);
BED = zeros(Tdata, 1);
ICU_nation(2:Tdata + 1, 1) = Data(:, 22);
BED(1:Tdata, 1) = Data(:, 23);
ICU_pref(2:Tdata + 1, 1) = Data(:, 30); %Data(:, 21);
%--- Import hospitalized patient data ---%
hospital = zeros(Tdata + 1, 1);
hospital(2:end) = Data(:, 25);
hospital(isnan(hospital)) = 0;
%--- Import SIRD data ---%
I_data = zeros(Tdata+1,1);
I_data(2:end) = Data(:,27);
R_data = zeros(Tdata+1,1);
R_data(2:end)= Data(:,28);
D_data = zeros(Tdata+1,1);
D_data(2:end) = Data(:,29);
S_data = [POP(1);POP] - (I_data + R_data + D_data);
dI_data = I_data(2:end) - I_data(1:end-1);
dR_data = R_data(2:end) - R_data(1:end-1);
dD_data = D_data(2:end) - D_data(1:end-1);
dS_data = S_data(2:end) - S_data(1:end-1);
N_pref = Data(:,26);
Severe_pref = zeros(Tdata+1,1);
Severe_pref(2:end) = Data(:,30);

POP0 = POP(1);          % initial population
ps = POP0/POP_jp;    % population share
xtick1 = 1:13:Tdata;
dateEN = datetime(dateD,'ConvertFrom','excel');
SimDate = dateD(end)+7:7:dateD(end)+7*SimPeriod;
SimDateEN = datetime(SimDate,'ConvertFrom','excel');
%--- Create Month-Week labels ---%
dateP = dateD(end)+7:7:dateD(end)+7*(SimPeriod+1);
date = [dateD;dateP'];
date = datetime(date,'ConvertFrom','excel');
MonthNumber = month(date);
WeekNumber = zeros(length(MonthNumber),1);
WeekNumber(1:2) = [3;4];
for i = 3:length(WeekNumber)
    if MonthNumber(i)~=MonthNumber(i-1)
        WeekNumber(i) = 1;
    else
        WeekNumber(i) = WeekNumber(i-1) + 1;
    end
end
MonthWeekJP = strings([length(MonthNumber),1]);
MonthWeekEN = strings([length(MonthNumber),1]);
for i = 1:length(MonthNumber)
    MonthWeekJP(i) = [num2str(date(i).Month) '月第' num2str(WeekNumber(i)) '週'];
    MonthWeekEN(i) = [datestr(date(i),'mmm') '-' num2str(WeekNumber(i)) 'w'];
end
M = 1+0.01*M;
TdataGDP = Tdata-sum(isnan(GDP));
if retroH_switch == 1
    RetroH = TdataGDP -4;
end
if isempty(find(SimDateEN == medical_start_date,1)) == 0
    medical_start = find(SimDateEN == medical_start_date);
else
    medical_start = 1;
end
elderly_start = find(SimDateEN == elderly_start_date);
VacStart = find(SimDateEN == datetime(2021,4,1));
End2020 = find(dateEN == datetime(2021,1,7));

betaT_temp_ini = betaT_temp_ini_vec(pindex);
start_beta = start_beta_vec(pindex);
beta_rho = beta_rho_vec(pindex);

delta_temp_ini = delta_temp_ini_vec(pindex);
delta_ICU_nation_temp_ini = delta_ICU_nation_temp_ini_vec(pindex);
delta_ICU_pref_temp_ini = delta_ICU_pref_temp_ini_vec(pindex);
delta_Hospital_temp_ini = delta_Hospital_temp_ini_vec(pindex);

start_delta = start_delta_vec(pindex);
start_delta_ICU_nation = start_delta_ICU_nation_vec(pindex);
start_delta_ICU_pref = start_delta_ICU_pref_vec(pindex);
start_delta_Hospital = start_delta_Hospital_vec(pindex);

delta_rho = delta_rho_vec(pindex);
delta_ICU_nation_rho = delta_ICU_nation_rho_vec(pindex);
delta_ICU_pref_rho = delta_ICU_pref_rho_vec(pindex);
delta_Hospital_rho = delta_Hospital_rho_vec(pindex);

var_ss = var_ss_vector(pindex);

th_on1 = th_on_vector(pindex)*7;         % threshold to place the state of emergency
th_on2 = th_on_vector2(pindex)*7;         % threshold to place the state of emergency
th_off1 = th_off_vector(pindex)*7;
th_off2 = th_off_vector2(pindex)*7;
th_off3 = th_off_vector3(pindex)*7;
alpha_on = alpha_on_vector(pindex); %(((ERN_on*(POP0/S(end))*((gammaT(1)+delta_average)/beta_avg)).^(1/k))-1)*(h(1)/h(2));
ICU_limit = ICU_limit_vec(pindex);

% Parameters for Vaccine Path
paces_ori = paces_ori_vec(pindex);
gradual_paces = gradual_paces_vec(pindex);
sw_vacpath = sw_vacpath_vec(pindex);
VT3share = VT3share_vec(pindex);
lag = lag_vec(pindex);
medical_duration = medical_duration_vec(pindex);
ind_date = find(date == ind_date_vec(pindex));

paces2_ori = paces2_ori_vec(pindex);
paces2 = paces2_ori*ps;
gradual_paces2 = gradual_paces2_vec(pindex);
ind_date2 = find(date == datetime(2021,6,24)); %when 職域接種 starts
date_slowdown = find(date == datetime(2021,8,26)) - Tdata; % when the vaccine paces slow down
date_slowdown2 = find(date == datetime(2021,9,30)) - Tdata;
date_slowdown3 = find(date == datetime(2021,10,28)) - Tdata;
paces3_ori = paces3_ori_vec(pindex);
paces3 = paces3_ori*ps;

paces4_ori = paces4_ori_vec(pindex);
paces4 = paces4_ori*ps;
paces5_ori = paces5_ori_vec(pindex);
paces5 = paces5_ori*ps;

ICU_limit_pref_vec = zeros(Tdata,1);
ICU_limit_pref_vec(find(dateEN == datetime(2020,12,03))) = 150;
ICU_limit_pref_vec(find(dateEN == datetime(2020,12,10)):find(dateEN == datetime(2020,12,10))+1) = 200;
ICU_limit_pref_vec(find(dateEN == datetime(2021,1,07))) = 220;
ICU_limit_pref_vec(find(dateEN == datetime(2021,1,7)):find(dateEN == datetime(2021,1,28))-1) = 250;
ICU_limit_pref_vec(find(dateEN == datetime(2021,1,28))) = 265;
ICU_limit_pref_vec(find(dateEN == datetime(2021,2,04)):find(dateEN == datetime(2021,2,18))-1) = 315;
ICU_limit_pref_vec(find(dateEN == datetime(2021,2,18)):find(dateEN == datetime(2021,3,18))-1) = 330;
ICU_limit_pref_vec(find(dateEN == datetime(2021,3,18)):find(dateEN == datetime(2021,4,29))-1) = 332;
ICU_limit_pref_vec(find(dateEN == datetime(2021,4,29)):find(dateEN == datetime(2021,7,08))-1) = 373;
ICU_limit_pref_vec(find(dateEN == datetime(2021,7,08)):Tdata) = 392;


Hospital_limit_vec = zeros(Tdata,1);
Hospital_limit_vec(1:find(dateEN == datetime(2021,4,29))-1) = 5048;
Hospital_limit_vec(find(dateEN == datetime(2021,4,29)):find(dateEN == datetime(2021,7,15))-1) = 5594;
Hospital_limit_vec(find(dateEN == datetime(2021,7,15))) = 5882;
Hospital_limit_vec(find(dateEN == datetime(2021,7,22)):Tdata) = 5967;

% BED_limit_vec = zeros(Tdata,1);



