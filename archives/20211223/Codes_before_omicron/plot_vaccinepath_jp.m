function plot_vaccinepath_jp(paces_ori,sw_vacpath,gradual_paces,V1_medical_ori,V2_medical_ori,medical_jp,V1_elderly_ori,V2_elderly_ori,elderly_jp,ordinary_jp,accept_share,E1,E2,D1,D2,lag,medical_duration,vaccine_figure_loop,SimPeriod,MonthWeekJP,MonthWeekEN,WeekNumber,Tdata,fs,ldfs,fn,language,iTH,TH,ft)

% if ~(iTH ~= 1 && vaccine_figure_loop == 0)
    ps = 1;
    paces = ps*paces_ori;
    vacpath = zeros(SimPeriod,1);
    vacpath(1+sw_vacpath:gradual_paces) = (paces/(gradual_paces-sw_vacpath)):(paces/(gradual_paces-sw_vacpath)):paces;
    vacpath(gradual_paces+1:end) = paces*ones(SimPeriod-gradual_paces,1);
    elderly_total = ps*elderly_jp;
    medical_total = ps*medical_jp;
    ordinary_total = ps*ordinary_jp;
    medical = medical_total*accept_share;
    elderly = elderly_total*accept_share;
    ordinary = ordinary_total*accept_share;
    elderly = elderly - (sum(V1_elderly_ori));
    delta_average = 0;
    POP_jp = 125710000;
    [~,~,VT] = vaccine_distribution_medical_ori(vacpath,medical,V1_medical_ori,V2_medical_ori,elderly,V1_elderly_ori,V2_elderly_ori, ordinary,elderly_total,delta_average,E1,E2,D1,D2,ps,POP_jp,lag,medical_duration);
% end

%% figure for vaccine path
if vaccine_figure_loop == 0
    if iTH == 1
        l = 1; % 1 for EN, 2 for JP
        lng = language{l};
        figname = string(['VaccineDistribution_' char(lng)]);
        f = figure('Name',figname);
        if l == 1
            plot_vaccinepath(figname,VT,V1_medical_ori,V2_medical_ori,V1_elderly_ori,V2_elderly_ori,SimPeriod,ps,MonthWeekEN,WeekNumber,Tdata,fs,ldfs,fn);
        elseif l ==2
            plot_vaccinepath(figname,VT,V1_medical_ori,V2_medical_ori,V1_elderly_ori,V2_elderly_ori,SimPeriod,ps,MonthWeekJP,WeekNumber,Tdata,fs,ldfs,fn);
        end
    end
elseif vaccine_figure_loop == 1
    l = 1; % 1 for EN, 2 for JP
    lng = language{l};
    figname = string(['VaccineDistribution_' char(lng)]);
    f = figure('Name',figname);
    if l == 1
        plot_vaccinepath(30+iTH,VT,V1_medical_ori,V2_medical_ori,V1_elderly_ori,V2_elderly_ori,SimPeriod,ps,MonthWeekEN,WeekNumber,Tdata,fs,ldfs,fn);
        subplot(1,2,1)
        title(string(['Number of New Vaccinations (',sprintf(ft,TH(iTH)),')']), 'FontSize',fs,'FontName',fn);
        subplot(1,2,2)
        title(string(['Number of Cumulative Vaccinations (',sprintf(ft,TH(iTH)),')']), 'FontSize',fs,'FontName',fn);
    elseif l == 2
        plot_vaccinepath(40+iTH,VT,V1_medical_ori,V2_medical_ori,V1_elderly_ori,V2_elderly_ori,SimPeriod,ps,MonthWeekJP,WeekNumber,Tdata,fs,ldfs,fn);
        subplot(1,2,1)
        title(string(['新規ワクチン接種本数（週ごと）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
        subplot(1,2,2)
        title(string(['累計ワクチン接種本数（週ごと）（',sprintf(ft,TH(iTH)),'）']), 'FontSize',fs,'FontName',fn);
    end
end
