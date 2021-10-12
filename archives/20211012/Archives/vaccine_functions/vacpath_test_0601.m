
% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata
%
% To investigate the effects of Indian Variant, make "variant_India" = 1.
%


clear variables
close all
iPC=0; % 0 for Mac, 1 for Windows
if iPC==1
     home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
    home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/';
    home = '/Users/koheimachi/Dropbox/fujii_nakata/Website/Codes/';
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 1;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 1;      % save back data
vaccine_figure_loop = 1; % =0 appear only once; =1 appear every loop;
beta_figure_loop = 0; % =0 appear only once; =1 appear every loop;
vaccine_disp_switch = 1; % =0 not display the summary of # of the vaccinated ; =1 display
variant_India = 1; % = 0 ... without Indian variant, = 1 ... with Indian variant
ICU_nation = 1; % = 1 use national definition (NHK data), = 0 use data from Tokyo Keizai
% in the "Figure" folder
fs = 20;            % common font size for many figures
ldfs = 8;           % legend font size for vaccine path
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
linecolor = {'red','blue','black','green'};
language = {'EN','JP'};
%======================================================================%

%================== Model Fixed Parameter Values ============================%
parameter

gradual_paces = 1;
lag = 6;


%================== Parameter Values (Prefecture Specific) ============================%
prefecture_parameter



%    covid = importdata([home '\Covid_weekly_newV.csv']);  % Import weekly Covid data by prefecture
   
pref = 'Tokyo';

covid = importdata([home 'Covid_weekly_newV.csv']);  % Import weekly Covid data by prefecture
Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);

Tdata= size(Data,1);    % Data period in weeks

dateD = Data(:,1) + 21916;
dateEN = datetime(dateD,'ConvertFrom','excel');

pref = 'Tokyo';

ps = 0.11;    % population share
POP0 = 14000000;

delta_average = 0.2;

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


    [V1_medical_ori, V1_medical, V2_medical_ori, V2_medical,...
    V1_elderly_ori, V1_elderly, V2_elderly_ori, V2_elderly, ...
    vs_MT,vs_ET,vs_M1,vs_E1,vs_M2,vs_E2] ...
    = ImportVaccineData(home,iPC,Data,pref,dateEN,ps,vaccine_disp_switch);
%---------------------

paces = ps*paces_ori; %3600000;
vacpath_ori = zeros(SimPeriod,1);
vacpath_ori(1+sw_vacpath:gradual_paces) = (paces/(gradual_paces-sw_vacpath)):(paces/(gradual_paces-sw_vacpath)):paces;
vacpath_ori(gradual_paces+1:end) = paces*ones(SimPeriod-gradual_paces,1);
elderly_total = ps*elderly_jp;
medical_total = ps*medical_jp;
ordinary_total = ps*ordinary_jp;
%medical = medical_total*accept_share;
medical = medical_total;
elderly = elderly_total*accept_share;
ordinary = ordinary_total*accept_share;
elderly = elderly - (sum(V1_elderly));


T = length(vacpath_ori);
VT = zeros(T,6);

VT(1:lag,6) = VT(1:lag,6)+(sum(V1_medical)-sum(V2_medical) )/lag;
VT(1:medical_duration,5) =VT(1:medical_duration,5) + (medical - sum(V1_medical))/medical_duration;
VT(1+lag:medical_duration+lag,6) = VT(1+lag:medical_duration+lag,6)+  (medical-sum(V1_medical))/medical_duration;
VT(1:lag,2) = VT(1:lag,2)+(sum(V1_elderly)-sum(V2_elderly) )/lag;

sumV = sum(VT,2);

%umV(1) = 1.0e+05 * 3;

for i = 1:medical_duration + lag
    if sumV(i) > vacpath_ori(i) 
        error('Reported vaccine > vacpath_ori')
    end
end

vacpath = vacpath_ori- sumV;


for t = 1:T
    CV = cumsum(VT);
    if t <= lag
        if CV(t,1) + vacpath(t) <= elderly
            VT(t,1) = vacpath(t);
        elseif CV(t,1) + vacpath(t) > elderly && CV(t,1) < elderly
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1);
        elseif CV(t,1) >= elderly
            if CV(t,3) + vacpath(t) <= ordinary
                VT(t,3) = vacpath(t);
            elseif CV(t,3) + vacpath(t) > ordinary && CV(t,3) < ordinary
                VT(t,3) = ordinary - CV(t,3);
            end
        end
    elseif CV(t,1) + vacpath(t) <= elderly
        if VT(t-lag,1) == 0
            VT(t,1) = vacpath(t);
        elseif VT(t-lag,1) > 0
            VT(t,2) = VT(t-lag,1);
            VT(t,1) = vacpath(t) - VT(t,2);
        end
    elseif CV(t,1) + vacpath(t) > elderly && CV(t,1) + vacpath(t) - VT(t-lag,1) <= elderly
        if VT(t-lag,1) == 0
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1);
        elseif VT(t-lag,1) > 0
            VT(t,2) = VT(t-lag,1);
            VT(t,1) = vacpath(t) - VT(t,2);
        end
    elseif CV(t,1) + vacpath(t) - VT(t-lag,1) >= elderly && CV(t,1) < elderly
        if VT(t-lag,1) == 0
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1);
        elseif VT(t-lag,1) > 0
            VT(t,2) = VT(t-lag,1);
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1) - VT(t,2);
        end
    elseif CV(t,1) >= elderly
        if VT(t-lag,1) > 0
            if VT(t-lag,3) == 0
                VT(t,2) = VT(t-lag,1);
                VT(t,3) = vacpath(t) - VT(t,2);
            else
                VT(t,2) = VT(t-lag,1);
                VT(t,4) = VT(t-lag,3);
                VT(t,3) = vacpath(t) - VT(t,2) - VT(t,4);
            end
        elseif CV(t,3) < ordinary
            if VT(t-lag,3) > 0
                VT(t,4) = VT(t-lag,3);
                VT(t,3) = vacpath(t) - VT(t,4);
            else
                VT(t,3) = vacpath(t);
            end
        elseif CV(t,3) >= ordinary
            if VT(t-lag,3) > 0
                VT(t,4) = VT(t-lag,3);
            end
        end
    end
end

%VT(1:lag,6) = VT(1:lag,6)+(sum(V1_medical)-sum(V2_medical) )/lag;
%VT(1:medical_duration,5) =VT(1:medical_duration,5) + (medical - sum(V1_medical))/medical_duration;
%VT(1+lag:medical_duration+lag,6) = VT(1+lag:medical_duration+lag,6)+  (medical-sum(V1_medical))/medical_duration;
%VT(1:lag,2) = VT(1:lag,2)+(sum(V1_elderly)-sum(V2_elderly) )/lag;
V1 = E1*(VT(:,1)+VT(:,3)+VT(:,5))+(E2-E1)*(VT(:,2)+VT(:,4)+VT(:,6));
V2 = D1*VT(:,1)+(D2-D1)*VT(:,2);
V2 = [0;0;V2(1:end-2)];
V = [0;0;V1(1:end-2)];
% V_ord = D1*(VT(:,3)+VT(:,5))+D2*(VT(:,4)+VT(:,6));
V_ord = D1*(VT(:,3)+VT(:,5))+(D2-D1)*(VT(:,4)+VT(:,6));
V_ord = [0;0;V_ord(1:end-2)];

delta_ss = delta_average*(0.1063/1.53);
% delta_ss = delta_average*(0.09/1.28);

deltaT = (1-(cumsum(V2)/elderly_total))*(delta_average-delta_ss)+(1-(cumsum(V_ord)/(POP0-elderly_total)))*delta_ss;



plot_vaccinepath(200,VT,V1_medical,V2_medical,V1_elderly,V2_elderly,SimPeriod,ps,MonthWeekJP,WeekNumber,Tdata,fs,ldfs,fn);