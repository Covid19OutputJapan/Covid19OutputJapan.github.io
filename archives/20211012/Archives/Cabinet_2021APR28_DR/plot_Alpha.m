function plot_Alpha(alpha,AlphaPath,TH,TH_index,MonthWeek,WeekNumber,Tdata,linecolor,ft,fs,fn,l)
%plot Alpha
th_wave = 1;
for i = 1:length(TH)
    if abs(TH(i) - TH_index(th_wave)) < 0.0001
        plot([1-alpha;1-AlphaPath(:,i)]*100,linecolor{th_wave},'LineWidth',2,'DisplayName',sprintf(ft,TH(i)))
        th_wave = th_wave + 1;
        if th_wave > length(TH_index)
            th_wave = length(TH_index);
        end
    else
        plot([1-alpha;1-AlphaPath(:,i)]*100,'--','LineWidth',0.3,'DisplayName',sprintf(ft,TH(i)))
    end
    hold on
end
grid on
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
ax = gca;
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
xticks(find(WeekNumber==1))
xtickangle(45)
ytickformat('%,3.0f')
xticklabels(MonthWeek(WeekNumber==1))
xlim([0 inf])
ylim([70 100])
if l == 1
    xlabel('Weeks)','FontSize',fs)
    ylabel('1-alpha','FontSize',fs)
    title('1-alpha Path','FontSize',fs,'FontWeight','normal')
elseif l == 2
    xlabel('Weeks','FontSize',fs,'FontName',fn)
    ylabel('alpha','FontSize',fs,'FontName',fn)
    title('1-alpha path','FontSize',fs,'FontWeight','normal','FontName',fn)
end
