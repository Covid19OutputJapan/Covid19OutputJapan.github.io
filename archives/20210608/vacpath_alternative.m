%%%%%%%%%%% Vaccine Distribution Specification %%%%%%%%%%%%
clear variables; close all;
% home = '/Users/ymaeda/Documents/ym_doc/Tokyo_Univrsity_MA/Reserach Assistant/fujii_nakata/Codes/';
home = '/Users/ymaeda/Dropbox/fujii_nakata/Website/Codes/';
cd(home);
iPC = 0;

parameter
%========================Parameter to control =======================%
% path_choice = "government";
path_choice = "baseline";

if path_choice == "baseline"
    paces_ori = 4200000;
    gradual_paces = 1;
    lag = 4;
    medical_duration = 2;
elseif path_choice == "government"
    paces_ori = 7000000;
    gradual_paces = 3;
    lag = 3;
    medical_duration = 2;
else
    error("Not correct specificaiton of vaccine path name")
end
%========================Parameter to control =======================%



%====================== Program parameter values ======================%
vaccine_disp_switch = 0;
figure_save = 0;
fs = 20;            % common font size for many figures
ldfs = 8;           % legend font size for vaccine path
ICU_nation = 1;
fn = 'YuGothic';
%================== Model Fixed Parameter Values ============================%

%================== Parameter Values (Prefecture Specific) ============================%
pref = 'Tokyo';

%--- Import data ---%
% Covid data are recorded at weekly frequencies (Mon-Sun)
% The first week start on January 20 (Mon), 2020
import_prefecture

%--- Construct weekly vaccine data ---%　
[V1_medical_ori, V1_medical, V2_medical_ori, V2_medical,...
V1_elderly_ori, V1_elderly, V2_elderly_ori, V2_elderly, ...
vs_MT,vs_ET,vs_M1,vs_E1,vs_M2,vs_E2] ...
= ImportVaccineData(home,iPC,Data,pref,dateEN,ps,vaccine_disp_switch);


%--- Constructing the reference level of output ---%
[potentialGDP, referenceGDP, alpha] = construct_GDP(GDP,TdataGDP);

%--- Regress mobility on alpha to estimate the elasticity h ---%　関数化？
%[Malt,h_all,h_all_se,h,h_se] = estimate_h(M,alpha,TdataGDP,RetroH,hconstant);
[Malt,h_all,h_all_se,h_ori,h_se] = estimate_h(M,alpha,TdataGDP,RetroH,hconstant);


%--- Import ICU data (Option 2) ---%
ICU = zeros(Tdata+1,1);
BED = zeros(Tdata,1);
if ICU_nation == 1
    ICU(2:Tdata+1,1) = Data(:,22);
    BED(1:Tdata,1) = Data(:,23);
else
    ICU(2:Tdata+1,1) = Data(:,21);
end


%--- Compute the history of S, I, R, D in the data period ---%
[S,I,R,D]...
    = SIRD(Tdata,POP0,N,E1,E2,...
            V1_elderly,V1_medical,V2_elderly,V2_medical,...
            gamma,dD,TdataGDP,referenceGDP,alpha);



[delta,beta_tilde,ERN,beta,ICU_inflow,...
        gammaT,delta_average,ICU_inflow_avg]...
            = Time_Series_Average(S,I,D,ICU,dD,N,Tdata,SimPeriod,...
                RetroPeriod,POP0,gamma,hconstant,h_all,alpha,k,...
                gamma_ICU,ICU_adjustment);








%---------------------
paces_now = V1_elderly(end)+V2_elderly(end)+V1_medical(end)+V2_medical(end);
paces = ps*paces_ori; %3600000;
vacpath_ori = zeros(SimPeriod,1);
% vacpath_ori(1+sw_vacpath:gradual_paces) = (paces/(gradual_paces-sw_vacpath)):(paces/(gradual_paces-sw_vacpath)):paces;
vacpath_ori(1+sw_vacpath:gradual_paces) = paces_now+((paces-paces_now)/(gradual_paces-sw_vacpath)):((paces-paces_now)/(gradual_paces-sw_vacpath)):paces;
vacpath_ori(gradual_paces+1:end) = paces*ones(SimPeriod-gradual_paces,1);
elderly_total = ps*elderly_jp;
medical_total = ps*medical_jp;
ordinary_total = ps*ordinary_jp;
medical = medical_total;
elderly = elderly_total*accept_share;
ordinary = ordinary_total*accept_share;
elderly = elderly - (sum(V1_elderly));

T = length(vacpath_ori);
VT = zeros(T,6);

VT(1:lag,6) ...
    = VT(1:lag,6)+(sum(V1_medical)-sum(V2_medical) )/lag; %MEDICAL second shots in the first lag periods
VT(1:medical_duration,5) ...
    =max(0,VT(1:medical_duration,5) ...
    + (medical - sum(V1_medical))/medical_duration); %MEDICAL first shots in the first lag periods
VT(lag+1:medical_duration+lag,6) ...
    = max(0,(medical-sum(V1_medical))/medical_duration); %distribute the remaining of MEDICAL-second shots
VT(1:lag,2) = V1_elderly(end-lag+1:end); % second ELDERLY shots in the first few weeks of simulation
sumV = sum(VT,2); % Current total shots at each week after calculating # of vaccination dependent on before-simulation values


% VT(1:lag,6) = VT(1:lag,6)+(sum(V1_medical)-sum(V2_medical) )/lag;
% VT(1:medical_duration,5) =VT(1:medical_duration,5) + (medical - sum(V1_medical))/medical_duration;
% VT(1+lag:medical_duration+lag,6) = VT(1+lag:medical_duration+lag,6)+  (medical-sum(V1_medical))/medical_duration;
% VT(1:lag,2) = VT(1:lag,2)+(sum(V1_elderly)-sum(V2_elderly) )/lag;
% 
% sumV = sum(VT,2);

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

% For debugging purpose
disp('Total Vaccinated Elderly for the First Shot')
disp(sum(VT(:,1))+sum(V1_elderly))

disp('Total Vaccinated Elderly for the Second Shot')
disp(sum(VT(:,2))+sum(V2_elderly))

disp('# of elderly * accept_share')
disp(elderly_total*accept_share)

disp('Total Vaccinated Medical for the First Shot')
disp(sum(VT(:,5))+sum(V1_medical))

disp('Total Vaccinated Medical for the Second Shot')
disp(sum(VT(:,6))+sum(V2_medical))

disp('# of medical')
disp(medical)



if any(VT(:) < 0)
    error('Negative Vaccination Value Exists')
end
plot_deltapath(201,delta,deltaT,deltaT(1),MonthWeekJP,WeekNumber,Tdata,fs,fn,1);
ylim([0 1.1])
plot_vaccinepath_ym(200,VT,V1_medical,V2_medical,V1_elderly,V2_elderly,date,ps,MonthWeekJP,WeekNumber,Tdata,fs,ldfs,fn);

if figure_save == 1
    saveas(figure(200), [home 'Figures/'  'Vaccine_path_Test_Alternative_' num2str(paces_ori) '.png']);
end

