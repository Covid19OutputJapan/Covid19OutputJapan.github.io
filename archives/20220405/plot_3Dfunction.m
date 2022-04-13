function plot_3Dfunction(past_path, sim_path, iX, ...
                       WeekNumber, YearMonth, xmin, xmax, ...
                       fn, fs, lgdfs, axfs,yft,...
                       lgdLocation, column_num, l, title_vec, ...
                       lineWidth,linecolor, LineStyles, lineName)

Tdata = length(past_path(:,1));
SimPeriod = length(sim_path(:,1,1));
nY = length(sim_path(1,1,:));

for iY = 1:nY
        plot([past_path;sim_path(:,iX,iY)],'LineWidth',lineWidth(iY),...
            'Color',linecolor{iY}, 'LineStyle', LineStyles{iY},'DisplayName',lineName{iY})
        hold on
end
plot([past_path; nan(SimPeriod,1)],'k', 'LineStyle', '-','LineWidth',2.0,'HandleVisibility','off')
hold on
xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
lgd = legend;
lgd.FontSize = lgdfs;
lgd.Location = lgdLocation;
lgd.NumColumns = column_num;
title(title_vec(l),'FontSize',fs,'FontWeight','normal','FontName',fn)
ax = gca;
ax.YAxis.FontSize = axfs;
ax.XAxis.FontSize = axfs;
ax.YAxis.Exponent = 0;
ytickformat(yft)
xticks(find(WeekNumber == 1))
xtickangle(45)
xlim([xmin, xmax])
xticklabels(YearMonth(xticks,l))
Simmax = max(sim_path(1:(xmax - Tdata),iX, :,:), [], 'all');
Datamax = max(past_path(xmin:end));
ymax = max(Simmax, Datamax) * 1.1;
if ymax>0
    ylim([0 ymax]) 
end