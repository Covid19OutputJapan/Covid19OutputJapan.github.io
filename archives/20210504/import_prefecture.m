%--- Import Prefecture-specifc data ---%
pref = PrefVector{pindex};        % prefecture to be analyzed
prefGDP = GDPVector(pindex);

% Covid data are recorded at weekly frequencies (Mon-Sun)
% The first week start on January 20 (Mon), 2020
if iPC==1
    covid = importdata([home '\Covid_weekly.csv']);  % Import weekly Covid data by prefecture
else
    covid = importdata([home 'Covid_weekly.csv']);  % Import weekly Covid data by prefecture
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
Tdata= size(Data,1);    % Data period in weeks
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

th_on1 = th_on_vector(pindex)*7;         % threshold to place the state of emergency
th_on2 = th_on_vector2(pindex)*7;         % threshold to place the state of emergency
th_off1 = th_off_vector(pindex)*7;
th_off2 = th_off_vector2(pindex)*7;
th_off3 = th_off_vector3(pindex)*7;
ERN_on = ERN_on_vector(pindex);
betaT_temp_ini = betaT_temp_ini_vec(pindex);
var_initial = var_initial_vector(pindex);

