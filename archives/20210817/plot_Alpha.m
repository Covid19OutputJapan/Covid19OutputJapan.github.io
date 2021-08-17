function plot_Alpha(alpha,AlphaPath,TH,TH_index,MonthWeek,WeekNumber,MonthNumber,Tdata,linecolor,ft,fs,fn,l)
%plot Alpha
th_wave = 1;
for i = 1:length(TH)
    if abs(TH(i) - TH_index(th_wave)) < 0.0001
        if i == 2
            plot([1-alpha;1-AlphaPath(:,i)]*100,linecolor{th_wave},'LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
        else
            plot([1-alpha;1-AlphaPath(:,i)]*100,'b-','LineWidth',0.5,'DisplayName',sprintf(ft,TH(i)))
        end
        th_wave = th_wave + 1;
        if th_wave > length(TH_index)
            th_wave = length(TH_index);
        end
    else
%         plot([1-alpha;1-AlphaPath(:,i)]*100,'--','LineWidth',0.3,'DisplayName',sprintf(ft,TH(i)))
    end
    hold on
end
plot([1-alpha;nan(length(AlphaPath(:,1)),1)]*100,'-k','LineWidth',2,'HandleVisibility','off')
grid on
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
ax = gca;
ax.YAxis.FontSize = 12;
ax.XAxis.FontSize = 12;
xticks(find(WeekNumber==1 & abs(1 - mod(MonthNumber,3)) < 0.01))
% xticks(find(WeekNumber==1))
xtickangle(45)
ytickformat('%,3.0f')
xticklabels(MonthWeek(xticks))
xlim([0 inf])
ylim([80 100])