function [delta_plot] = plot_deltapath(fignum,delta,deltaT,delta_average,MonthWeek,WeekNumber,Tdata,fs,fn)

figure(fignum)
plot([delta; deltaT]./delta_average,'-k','LineWidth',1.5);
xtickangle(45)
xticks(find(WeekNumber==1))
xlim([Tdata+1 Tdata+56])
xticklabels(MonthWeek(WeekNumber==1))
title('致死率（現在のレベルで標準化）','FontName',fn)
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
legend;