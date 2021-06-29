function [BackDataN, BackDataAlpha, BackDataERN]...
    = plot_SimN(TH,TH_index,N,NPath,alpha,AlphaPath,ERN,SimERN,THonPath,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,fs,fn,ft,l,fs_legend)
%plot Simulated New Cases
th_wave = 1;
for i = 1:length(TH)
    if abs(TH(i) - TH_index(th_wave)) < 0.0001
        plot([N;NPath(:,i)]/7,linecolor{th_wave},'LineWidth',2)
        hold on 
%         plot([N;NPath(:,i)]/7,linecolor{th_wave},'LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
        BackDataN(:,th_wave) = [N(end-7:end);NPath(:,i)];
        BackDataAlpha(:,th_wave) = [alpha(end-7:end);AlphaPath(:,i)];
        BackDataERN(:,th_wave) = [ERN(end-7:end);SimERN(:,i)];
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
hold on 
% plot([nan(length(N),1);THonPath(:,i)]/7,'-k','LineWidth',0.1)
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
xline(74,'LineWidth',1.5,'HandleVisibility','off');
ax = gca;
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
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
    legend('4ヶ月','6ヶ月','8ヶ月','FontSize',fs_legend,'FontName',fn,'Location','northwest');
end
% lgd = legend;
% lgd.NumColumns = 2;
% lgd.Location = 'northwest';
xtickangle(45)
xlim([Tdata-7 Tdata+25])

