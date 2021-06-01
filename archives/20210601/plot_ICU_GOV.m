function BackDataICU = ...
    plot_ICU_GOV(TH,TH_index,ICU,SimICU,ICU_limit,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,fs,fn,ft,l)
th_wave = 1;
for i = 1:length(TH)
    if abs(TH(i) - TH_index(th_wave)) < 0.0001        
        p = plot([ICU(2:end);SimICU(2:end,i)],linecolor{th_wave},'LineWidth',1.5,'DisplayName',sprintf(ft,TH(i)));
        p.LineStyle = '-';
%         hold on 
%         yline(TH(i),linecolor{th_wave},TH(i))

        BackDataICU(:,th_wave) = [ICU(end-7:end);SimICU(:,i)];
        th_wave = th_wave + 1;
        if th_wave > length(TH_index)
            th_wave = length(TH_index);
        end
    else
        plot([ICU(2:end);SimICU(2:end,i)], '--','LineWidth',0.3)
        hold on
        colororder('default')
%         yline(TH(i),'-',TH(i))
%         colororder('default')
    end
    hold on
    plot(ICU(2:end),'-k','LineWidth',1.5);
end

xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
xline(74,'LineWidth',1.5,'HandleVisibility','off');

% yline(ICU_limit,'-g',num2str(ICU_limit),'LineWidth',2.5,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
yline(ICU_limit/2,'--k',"50%",'LineWidth',1.0,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
yline(ICU_limit,'--k',"100%",'LineWidth',1.0,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
yline(ICU_limit/5,'--k',"20%",'LineWidth',1.0,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');

% yline(th_off1,'-r',num2str(th_off1),'LineWidth',1.5,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
% yline(th_off2,'-b',num2str(th_off2),'LineWidth',1.5,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
% yline(th_off3,'-k',num2str(th_off3),'LineWidth',1.5,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');

ax = gca;
ax.YAxis.Color = 'k';
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
ax.YAxis.Exponent = 0;

%ytickformat('%,6.0f')
xticks(find(WeekNumber==1))
if l == 1
    title('Projected path of ICU','FontSize',fs,'FontWeight','normal')
    xticklabels(MonthWeekEN(WeekNumber==1))
    %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
elseif l == 2
    title('重症患者数の推移','FontSize',fs,'FontWeight','normal','FontName',fn)
    xticklabels(MonthWeekJP(WeekNumber==1))
end
% lgd = legend;
lgd.NumColumns = 2;
xtickangle(45)
xlim([Tdata-7 Tdata+31])
ylim([0 ICU_limit*1.2])
