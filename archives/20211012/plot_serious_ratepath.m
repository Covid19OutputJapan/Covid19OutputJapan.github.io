function [deltaICU_plot] = plot_serious_ratepath(fignum,delta,deltaICU,deltaT_ini,MonthWeek,WeekNumber,Tdata,fs,fn,iTH)

figure(fignum)
if iTH == 1
    plot([delta; deltaICU]./deltaT_ini,'-k','LineWidth',1.5,'DisplayName',"基本見通し");
elseif iTH == 2
    plot([delta; deltaICU]./deltaT_ini,'-r','LineWidth',1.5,'DisplayName',"悲観見通し");
else
    plot([delta; deltaICU]./deltaT_ini,'-b','LineWidth',1.5,'DisplayName',"希望見通し");
end
hold on
xtickangle(45)
xticks(find(WeekNumber==1))
xlim([Tdata+1 Tdata+50])
xticklabels(MonthWeek(WeekNumber==1))
title('重症化率（現在のレベルで標準化）','FontName',fn)
ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
if iTH == 2
    lgd = legend;
    lgd.FontSize = fs;
    lgd.FontName = fn;
end