function [delta_plot] = plot_deltapath(fignum,delta,deltaT,deltaT_ini,MonthWeek,WeekNumber,Tdata,fs,fn,iTH)

figure(fignum)
if iTH == 1
    plot([delta; deltaT]./deltaT_ini,'-k','LineWidth',1.5,'DisplayName',"基本見通し");
else
    plot([delta; deltaT]./deltaT_ini,'-b','LineWidth',1.5,'DisplayName',"希望見通し");
end
hold on
xtickangle(45)
xticks(find(WeekNumber==1))
xlim([Tdata+1 Tdata+50])
xticklabels(MonthWeek(WeekNumber==1))
title('致死率（現在のレベルで標準化）','FontName',fn)
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
if iTH == 2
    lgd = legend;
    lgd.FontSize = fs;
    lgd.FontName = fn;
end