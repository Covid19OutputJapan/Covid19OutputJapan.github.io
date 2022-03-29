% %--- Import Prefecture-specifc data ---%
% Covid data are recorded at weekly frequencies (Mon-Sun)
% The first week start on January 20 (Mon), 2020
table = readtable([DATA_PATH 'Covid_weekly_NewV.csv']);
tablepref = table(strcmp(table.prefecture, PREF), :);
clear table

dateD       = tablepref.date + 21916;
N           = tablepref.positive; %New cases %or tablepref.newpositive 
dD          = tablepref.death; %New deaths
M           = tablepref.mobility; %Mobility
POP         = tablepref.population; %Population
GDP         = tablepref.gop/(tablepref.gop(1)) * 100; % GDO, Normalization: 100 = 2020JAN
Tdata       = size(dateD, 1); % Data period in weeks

ICU_nation                  = zeros(Tdata + 1, 1);
ICU_pref                    = zeros(Tdata + 1, 1);
hospital                    = zeros(Tdata + 1, 1);
ICU_nation(2:Tdata+1, 1)    = tablepref.ICU_NHK;
ICU_nation(112)=619;
ICU_pref(2:Tdata + 1, 1)    = tablepref.severe; %Data(:, 21);
hospital(2:end)             = tablepref.hospital;
hospital(isnan(hospital))   = 0;
BED                         = tablepref.all_bed;
BED(111)=BED(111-1);

M = 1 + 0.01 * M;
TdataGDP = Tdata - sum(isnan(GDP));
if retroH_switch == 1
    RetroH = TdataGDP -4;
end

%--- Import SIRD data ---%
I_data          = zeros(Tdata + 1, 1);
R_data          = zeros(Tdata + 1, 1);
D_data          = zeros(Tdata + 1, 1);
I_data(2:end)   = tablepref.I;
R_data(2:end)   = tablepref.R;
D_data(2:end)   = tablepref.D;
S_data          = [POP(1); POP] - (I_data + R_data + D_data);
dI_data         = I_data(2:end) - I_data(1:end - 1);
dR_data         = R_data(2:end) - R_data(1:end - 1);
dD_data         = D_data(2:end) - D_data(1:end - 1);
dS_data         = S_data(2:end) - S_data(1:end - 1);

POP0            = POP(1); % initial population
ps              = POP0 / POP_jp; % population share
xtick1          = 1:13:Tdata;
dateEN          = datetime(dateD, 'ConvertFrom', 'excel');
SimDate         = dateD(end) + 7:7:dateD(end) + 7 * SimPeriod;
SimDateEN       = datetime(SimDate, 'ConvertFrom', 'excel');
%--- Create Month-Week labels ---%
dateP           = dateD(end) + 7:7:dateD(end) + 7 * (SimPeriod + 1);
date            = [dateD; dateP'];
date            = datetime(date, 'ConvertFrom', 'excel');
MonthNumber     = month(date);
WeekNumber      = zeros(length(MonthNumber), 1);
WeekNumber(1:2) = [3; 4];

for i = 3:length(WeekNumber)
    if MonthNumber(i) ~= MonthNumber(i - 1)
        WeekNumber(i) = 1;
    else
        WeekNumber(i) = WeekNumber(i - 1) + 1;
    end
end
YearMonthWeekJP = strings([length(MonthNumber), 1]);
YearMonthWeekEN = strings([length(MonthNumber), 1]);

for i = 1:length(MonthNumber)
    YearMonthWeekJP(i) = [num2str(date(i).Year) '年' num2str(date(i).Month) '月第' num2str(WeekNumber(i)) '週'];
    YearMonthWeekEN(i) = [datestr(date(i), 'yyyy') '-' datestr(date(i), 'mmm') '-' num2str(WeekNumber(i)) 'w'];
end

YearMonthJP = strings([length(MonthNumber), 1]);
YearMonthEN = strings([length(MonthNumber), 1]);

for i = 1:length(MonthNumber)
    YearMonthJP(i) = [num2str(date(i).Year) '年' num2str(date(i).Month) '月'];
    YearMonthEN(i) = [datestr(date(i), 'yyyy') '-' datestr(date(i), 'mmm')];
end

YearMonth       = [YearMonthEN,     YearMonthJP];
YearMonthWeek   = [YearMonthWeekEN, YearMonthWeekJP];

% Vector for maximum number of beds for severe cases
ICU_limit_pref_vec = zeros(Tdata,1);
ICU_limit_pref_vec(find(dateEN == datetime(2020,12,03))) = 150;
ICU_limit_pref_vec(find(dateEN == datetime(2020,12,10)):find(dateEN == datetime(2021,1,7))- 1) = 200;
ICU_limit_pref_vec(find(dateEN == datetime(2021,1,7))) = 220;
ICU_limit_pref_vec(find(dateEN == datetime(2021,1,7)):find(dateEN == datetime(2021,1,28))-1) = 250;
ICU_limit_pref_vec(find(dateEN == datetime(2021,1,28))) = 265;
ICU_limit_pref_vec(find(dateEN == datetime(2021,2,04)):find(dateEN == datetime(2021,2,18))-1) = 315;
ICU_limit_pref_vec(find(dateEN == datetime(2021,2,18)):find(dateEN == datetime(2021,3,18))-1) = 330;
ICU_limit_pref_vec(find(dateEN == datetime(2021,3,18)):find(dateEN == datetime(2021,4,29))-1) = 332;
ICU_limit_pref_vec(find(dateEN == datetime(2021,4,29)):find(dateEN == datetime(2021,7,08))-1) = 373;
ICU_limit_pref_vec(find(dateEN == datetime(2021,7,08)):find(dateEN == datetime(2021,9,9))-1) = 392;
ICU_limit_pref_vec(find(dateEN == datetime(2021,9,9)):find(dateEN == datetime(2021,9,16))-1) = 492;
ICU_limit_pref_vec(find(dateEN == datetime(2021,9,16)):find(dateEN == datetime(2021,10,28))-1) = 503;
ICU_limit_pref_vec(find(dateEN == datetime(2021,10,28)):find(dateEN == datetime(2021,11,04))-1) = 372;
ICU_limit_pref_vec(find(dateEN == datetime(2021,11,04)):find(dateEN == datetime(2021,11,11))-1) = 366;
ICU_limit_pref_vec(find(dateEN == datetime(2021,11,11)):find(dateEN == datetime(2021,12,02))-1) = 356;
ICU_limit_pref_vec(find(dateEN == datetime(2021,12,02)):find(dateEN == datetime(2021,12,09))-1) = 352;
ICU_limit_pref_vec(find(dateEN == datetime(2021,12,02)):find(dateEN == datetime(2021,12,16))-1) = 349;
ICU_limit_pref_vec(find(dateEN == datetime(2021,12,16)):find(dateEN == datetime(2022,1,6))-1) = 346;
ICU_limit_pref_vec(find(dateEN == datetime(2022,1,6)):find(dateEN == datetime(2022,1,20))-1)   = 353;
ICU_limit_pref_vec(find(dateEN == datetime(2022,1,20)):find(dateEN == datetime(2022,1,27))-1) = 358;
ICU_limit_pref_vec(find(dateEN == datetime(2022,1,27)):find(dateEN == datetime(2022,2,3))-1) = 370;
ICU_limit_pref_vec(find(dateEN == datetime(2022,2,3)):find(dateEN == datetime(2022,2,10))-1) = 377;
ICU_limit_pref_vec(find(dateEN == datetime(2022,2,10)):find(dateEN == datetime(2022,3,3))-1) = 397;
ICU_limit_pref_vec(find(dateEN == datetime(2022,3,3)):Tdata) = 471;

Hospital_limit_vec = zeros(Tdata,1);
Hospital_limit_vec(1:find(dateEN == datetime(2021,4,29))-1) = 5048;
Hospital_limit_vec(find(dateEN == datetime(2021,4,29)):find(dateEN == datetime(2021,7,15))-1) = 5594;
Hospital_limit_vec(find(dateEN == datetime(2021,7,15))) = 5882;
Hospital_limit_vec(find(dateEN == datetime(2021,7,22)):find(dateEN == datetime(2021,9,9))-1) = 5967;
Hospital_limit_vec(find(dateEN == datetime(2021,9,9)):find(dateEN == datetime(2021,9,16))-1) = 6319;
Hospital_limit_vec(find(dateEN == datetime(2021,9,16)):find(dateEN == datetime(2021,10,7))-1) = 6583;
Hospital_limit_vec(find(dateEN == datetime(2021,10,7)):find(dateEN == datetime(2021,10,28))-1) = 6651;
Hospital_limit_vec(find(dateEN == datetime(2021,10,28)):find(dateEN == datetime(2021,11,04))-1) = 4836;
Hospital_limit_vec(find(dateEN == datetime(2021,11,04)):find(dateEN == datetime(2021,11,11))-1) = 4868;
Hospital_limit_vec(find(dateEN == datetime(2021,11,11)):find(dateEN == datetime(2021,11,18))-1) = 4834;
Hospital_limit_vec(find(dateEN == datetime(2021,11,18)):find(dateEN == datetime(2021,11,25))-1) = 4823;
Hospital_limit_vec(find(dateEN == datetime(2021,11,25)):find(dateEN == datetime(2021,12,02))-1) = 4820;
Hospital_limit_vec(find(dateEN == datetime(2021,12,02)):find(dateEN == datetime(2021,12,09))-1) = 4723;
Hospital_limit_vec(find(dateEN == datetime(2021,12,09)):find(dateEN == datetime(2021,12,16))-1)   = 4703;
Hospital_limit_vec(find(dateEN == datetime(2021,12,16)):find(dateEN == datetime(2021,12,23))-1)   = 4657;
Hospital_limit_vec(find(dateEN == datetime(2021,12,23)):find(dateEN == datetime(2021,12,30))-1)   = 4669;
Hospital_limit_vec(find(dateEN == datetime(2021,12,23)):find(dateEN == datetime(2022,1,6))-1)   = 4669;
Hospital_limit_vec(find(dateEN == datetime(2022,1,6)):find(dateEN == datetime(2022,1,13))-1)   = 4839;
Hospital_limit_vec(find(dateEN == datetime(2022,1,13)):find(dateEN == datetime(2022,1,20))-1)   = 4863;
Hospital_limit_vec(find(dateEN == datetime(2022,1,20)):find(dateEN == datetime(2022,1,27))-1)   = 5015;
Hospital_limit_vec(find(dateEN == datetime(2022,1,27)):find(dateEN == datetime(2022,2,3))-1)   = 6189;
Hospital_limit_vec(find(dateEN == datetime(2022,2,3)):find(dateEN == datetime(2022,2,10))-1)   = 6415;
Hospital_limit_vec(find(dateEN == datetime(2022,2,10)):find(dateEN == datetime(2022,2,24))-1) = 6529;
Hospital_limit_vec(find(dateEN == datetime(2022,2,24)):find(dateEN == datetime(2022,3,3))-1) = 6599;
Hospital_limit_vec(find(dateEN == datetime(2022,3,3)):find(dateEN == datetime(2022,3,17))-1) = 6815;
Hospital_limit_vec(find(dateEN == datetime(2022,3,17)):Tdata) = 6946;

newICU_limit_pref_vec = zeros(Tdata,1);
newICU_limit_pref_vec(1: find(dateEN == datetime(2022,3,3))-1) = 750;
newICU_limit_pref_vec(find(dateEN == datetime(2022,3,3)) : Tdata) = 804;

% Data Source: https://www.bousai.metro.tokyo.lg.jp/taisaku/saigai/1013388/index.html 
% 02　感染状況・医療提供体制の分析




