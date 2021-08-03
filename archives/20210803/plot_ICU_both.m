function [BackDataICU_nation, BackDataICU_pref] = ...
    plot_ICU_both(TH,TH_index,ICU_nation,SimICU_nation,ICU_limit_nation,...
    ICU_pref,SimICU_pref,ICU_limit_pref,MonthWeekEN,MonthWeekJP,WeekNumber,...
    Tdata,linecolor,fs,fn,ft,l)
th_wave = 1;
linename_nation = {'Naitonal Standard','国基準'};
linename_pref = {'Local Standard','地方基準'};
for i = 1:length(TH)
    if abs(TH(i) - TH_index(th_wave)) < 0.0001
        yyaxis left
        plot([ICU_nation(2:end);SimICU_nation(2:end,i)],linecolor{th_wave},...
        'LineStyle','-','LineWidth',1.5,'DisplayName',[linename_nation{l} sprintf(ft,TH(i))]);
        hold on

        BackDataICU_nation(:,th_wave) = [ICU_nation(end-7:end);SimICU_nation(2:end,i)];
        th_wave = th_wave + 1;
        if th_wave > length(TH_index)
            th_wave = length(TH_index);
        end
    else
        plot([ICU_nation(2:end);SimICU_nation(2:end,i)], ':','LineWidth',0.3)
        hold on
        colororder('default')
    end
    hold on
end
plot(ICU_nation(2:end),'-k','LineWidth',1.5,'HandleVisibility','off');
yline(ICU_limit_nation/2,'--k',"50%",'LineWidth',1.0,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
yline(ICU_limit_nation,'--k',"100%",'LineWidth',1.0,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
ylim([0 ICU_limit_nation*1.1])

th_wave = 1;
for i = 1:length(TH)
    if abs(TH(i) - TH_index(th_wave)) < 0.0001
        yyaxis right
        plot([ICU_pref(2:end);SimICU_pref(2:end,i)],linecolor{th_wave},'LineStyle','--',...
        'LineStyle','-.','LineWidth',1.5,'DisplayName',[linename_pref{l} sprintf(ft,TH(i))]);

        BackDataICU_pref(:,th_wave) = [ICU_pref(end-7:end);SimICU_pref(2:end,i)];
        th_wave = th_wave + 1;
        if th_wave > length(TH_index)
            th_wave = length(TH_index);
        end
    else
        plot([ICU_pref(2:end);SimICU_pref(2:end,i)], ':','LineWidth',0.3)
        colororder('default')
%         yline(TH(i),'-',TH(i))
%         colororder('default')
    end
    hold on
end
plot(ICU_pref(2:end),'-k','LineWidth',1.5,'HandleVisibility','off');
yline(ICU_limit_pref/2,':k',"50%",'LineWidth',1.0,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
yline(ICU_limit_pref,':k',"100%",'LineWidth',1.0,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');

xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
% xline(74,'LineWidth',1.5,'HandleVisibility','off');
ylim([0 ICU_limit_pref*1.1])


ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
ax.YAxis(1).FontSize = fs;
ax.YAxis(2).FontSize = fs;
ax.XAxis.FontSize = fs;
ax.YAxis(1).Exponent = 0;
ax.YAxis(2).Exponent = 0;

%ytickformat('%,6.0f')
xticks(find(WeekNumber==1))
if l == 1
    title('Projected Path of ICU','FontSize',fs,'FontWeight','normal')
    xticklabels(MonthWeekEN(WeekNumber==1))
    %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
elseif l == 2
    title('重症患者数の推移','FontSize',fs,'FontWeight','normal','FontName',fn)
    xticklabels(MonthWeekJP(WeekNumber==1))
end
lgd = legend;
lgd.NumColumns = 2;
xtickangle(45)
xlim([Tdata-7 Tdata+25])
