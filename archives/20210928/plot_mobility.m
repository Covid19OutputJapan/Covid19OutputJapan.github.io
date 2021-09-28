function plot_mobility(Malt,alpha,Tdata,TdataGDP,MonthWeekJP,xtick1,fs,fs_legend)
yyaxis left
plot(Malt,'k','LineWidth',1.5)
ylabel('Mobility')
yyaxis right
hold on
plot((1-alpha)*100,'r-.','LineWidth',1.5)
hold on
plot((1-alpha(1:TdataGDP))*100,'b-.','LineWidth',1.5)
ylabel('GDP')
xlim([1 Tdata])
xticks(xtick1)
xticklabels(MonthWeekJP(xtick1))
legend('Mobility (L axis)','Imputed GDP (R axis)', 'Data GDP (R axis)','FontSize',fs_legend,'Location','southeast');
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
ax.YAxis(1).FontSize = fs;
ax.YAxis(2).FontSize = fs;
ax.XAxis.FontSize = fs;
xtickangle(45)

