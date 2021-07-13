function plot_function2(TH,TH_index,pastPath1,SimPath1,pastPath2,SimPath2,...
    MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,xline_ind,linecolor,lineName,...
    lineName2,eTitle,jTitle,fs,lgfs,fn,ft,yft,column,l,show_other,lgdLocation)
%plot Past and Simulated Paths, two variables
th_wave = 1;
for i = 1:length(TH)
    if abs(TH(i) - TH_index(th_wave)) < 0.0001
        plot([pastPath1;SimPath1(:,i)],linecolor{th_wave},'LineWidth',2,'DisplayName',lineName{i})
        hold on
        plot([pastPath2;SimPath2(:,i)],linecolor{th_wave},'LineStyle','--','LineWidth',2,'DisplayName',lineName2{i})
        th_wave = th_wave + 1;
        if th_wave > length(TH_index)
            th_wave = length(TH_index);
        end
    elseif show_other == 1
        plot([pastPath1;SimPath1(:,i)],'--','LineWidth',0.3,'DisplayName',lineName{i})
        hold on
        plot([pastPath2;SimPath2(:,i)],'--','LineWidth',0.3,'DisplayName',lineName{i})
    end
end
plot(pastPath1,'k','LineWidth',2,'HandleVisibility','off')
plot(pastPath2,'k','LineWidth',2,'LineStyle','--','HandleVisibility','off')
xline(Tdata,'LineWidth',2.0,'HandleVisibility','off');
xline(xline_ind,'LineWidth',1.0,'HandleVisibility','off');
ax = gca;
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
ax.YAxis.Exponent = 0;
ytickformat(yft)
xticks(find(WeekNumber==1))
if l == 1
    title(char(eTitle),'FontSize',fs,'FontWeight','normal')
    xticklabels(MonthWeekEN(WeekNumber==1))
    %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
elseif l == 2
    title(char(jTitle),'FontSize',fs,'FontWeight','normal','FontName',fn)
    xticklabels(MonthWeekJP(WeekNumber==1))
end
lgd = legend;
lgd.NumColumns = column;
lgd.FontSize = lgfs;
lgd.Location = lgdLocation;
xtickangle(45)
% xlim([Tdata-7 Tdata+31])
