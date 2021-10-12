function plot_SimN(TH,TH_index,N,NPath,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,fs,fn,ft,l)
%plot Simulated New Cases
th_wave = 1;
for i = 1:length(TH)
    if abs(TH(i) - TH_index(th_wave)) < 0.0001
        plot([N;NPath(:,i)]/7,linecolor{th_wave},'LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
        th_wave = th_wave + 1;
        if th_wave > length(TH_index)
            th_wave = length(TH_index);
        end
    else
        plot([N;NPath(:,i)]/7,'--','LineWidth',0.3,'DisplayName',sprintf(ft,TH(i)))
    end
    hold on
end
plot(N/7,'k','LineWidth',2,'HandleVisibility','off')
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
xticks(find(WeekNumber==1))
if l == 1
    title('Projected path of new cases','FontSize',fs,'FontWeight','normal')
    xticklabels(MonthWeekEN(WeekNumber==1))
    %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
elseif l == 2
    title('新規感染者数の推移','FontSize',fs,'FontWeight','normal','FontName',fn)
    xticklabels(MonthWeekJP(WeekNumber==1))
end
lgd = legend;
lgd.NumColumns = 2;
xtickangle(45)
xlim([Tdata-7 Tdata+36])

