function [Vplot] = plot_vaccinepath(fignum,VT,V1_medical,V2_medical,V1_elderly,V2_elderly,SimPeriod,ps,MonthWeek,WeekNumber,Tdata,fs,ldfs,fn)

Vplot = figure(fignum);
set(gcf,'Position',[100,100,1200,500])
Data_vac = VT(:,:);
medStuff1_nv = [V1_medical; Data_vac(:,5)]; % In future, you should change V1_w to V1past_med.
medStuff2_nv = [V2_medical; Data_vac(:,6)]; % In future, you should change V2_w to V2past_med.
elderly1_nv = [V1_elderly; Data_vac(:,1)];
elderly2_nv = [V2_elderly; Data_vac(:,2)];
others1_nv = [zeros(Tdata,1); Data_vac(:,3)];
others2_nv = [zeros(Tdata,1); Data_vac(:,4)];

% Find the cumulative number of vaccines for each agent
medStuff1_cv = cumsum(medStuff1_nv);
medStuff2_cv = cumsum(medStuff2_nv);
elderly1_cv = cumsum(elderly1_nv);
elderly2_cv = cumsum(elderly2_nv);
others1_cv = cumsum(others1_nv);
others2_cv = cumsum(others2_nv);

Area_nv = [medStuff1_nv,medStuff2_nv,elderly1_nv,elderly2_nv,others1_nv,others2_nv];
Area_cv = [medStuff1_cv,medStuff2_cv,elderly1_cv,elderly2_cv,others1_cv,others2_cv];
Area_nv = Area_nv / (ps*10000);
Area_cv = Area_cv / (ps*100000000);


for ifig = 1:1:2
    if ifig == 1
        subplot(1,2,1)
        Areagraph = bar(Area_nv,'stacked','LineStyle','none');
    else
        subplot(1,2,2)
        Areagraph = bar(Area_cv,'stacked','LineStyle','none');
    end
    
%     ldfs = fs;
    if ifig == 1
        ylim([0 850])
        ylabel('新規ワクチン接種本数（万本）','FontName',fn)
        legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3),Areagraph(2),Areagraph(1)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','医療従事者２本目','医療従事者１本目','FontSize',ldfs,'Location','NorthEast','FontName',fn);
        title('新規ワクチン接種本数（週ごと）', 'FontSize',fs,'FontName',fn);
    else
        ylim([0 2.5])
        ylabel('累計ワクチン接種本数（億本）','FontName',fn)
        legend([Areagraph(6),Areagraph(5),Areagraph(4),Areagraph(3),Areagraph(2),Areagraph(1)], 'その他２本目','その他１本目','高齢者２本目','高齢者１本目','医療従事者２本目','医療従事者１本目','FontSize',ldfs,'Location','NorthWest','FontName',fn);
        title('累計ワクチン接種本数（週ごと）', 'FontSize',fs,'FontName',fn);
    end
    xtickangle(45)
    xticks(find(WeekNumber==1))
    xlim([Tdata-7 Tdata+56])
    xticklabels(MonthWeek(WeekNumber==1))
    
    ax = gca;
    ax.YAxis.Color = 'k';
    ax.YAxis.FontSize = fs;
    ax.XAxis.FontSize = fs;
    Areagraph(3).FaceColor = [0.2 0.6 0.8];
    Areagraph(1).FaceColor = [0.6 0.6 0.6];
    Areagraph(5).FaceColor = [0.4 0.4 0.8];
    Areagraph(4).FaceColor = [0.2 0.6 0.8];
    Areagraph(2).FaceColor = [0.6 0.6 0.6];
    Areagraph(6).FaceColor = [0.4 0.4 0.8];
    Areagraph(1).FaceAlpha = 1;
    Areagraph(3).FaceAlpha = 0.7;
    Areagraph(5).FaceAlpha = 1;
    Areagraph(2).FaceAlpha = 0.8;
    Areagraph(4).FaceAlpha = 0.5;
    Areagraph(6).FaceAlpha = 0.8;
end
