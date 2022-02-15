function plot_4Dfunction(past_path, sim_path, iX, ...
                       WeekNumber, YearMonth, xmin, xmax, ...
                       fn, fs, lgdfs, axfs,yft,...
                       lgdLocation, column_num, l, title_vec, ...
                       lineWidth,linecolor, LineStyles, lineName)

Tdata = length(past_path(:,1));
SimPeriod = length(sim_path(:,1,1,1));
nY = length(sim_path(1,1,:,1));
nZ = length(sim_path(1,1,1,:));

for iY = 1:nY
    for iZ = 1:nZ
        plot([past_path;sim_path(:,iX,iY,iZ)],'LineWidth',lineWidth(iY,iZ),...
            'Color',linecolor{iY,iZ}, 'LineStyle', LineStyles{iY,iZ},'DisplayName',lineName{iY,iZ})
        hold on
    end
end
plot([past_path; nan(SimPeriod,1)],'k', 'LineStyle', '-','LineWidth',2.0,'HandleVisibility','off')
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
xticks(find(WeekNumber == 1))
xticklabels(YearMonth(xticks,l))
Simmax = max(sim_path(1:(xmax - Tdata),iX, :,:), [], 'all');
Datamax = max(past_path(xmin:end));
ylim([0 max(Simmax, Datamax) * 1.1]) %ylim([0 ICU_limit * 1.1])