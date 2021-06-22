function plot_ICU_N(TH,TH_index,N,NPath,ICU,SimICU,ICU_limit,MonthWeekEN,MonthWeekJP,WeekNumber,Tdata,linecolor,fs,fn,ft,l,th_off1,th_off2,th_off3)

th_wave = 1;
for i = 1:length(TH)
    if abs(TH(i) - TH_index(th_wave)) < 0.0001
        yyaxis left
        p1 = plot([ICU;SimICU(2:end,i)],linecolor{th_wave},'LineWidth',2,'DisplayName',sprintf(ft,TH(i)));
        p1.LineStyle = ':';
        hold on
        
        yyaxis right
        p2 = plot([N;NPath(:,i)]/7,linecolor{th_wave},'LineWidth',2,'DisplayName',sprintf(ft,TH(i)));
        p2.LineStyle = '-';
        hold on
        %         hold on
        %         yline(TH(i),linecolor{th_wave},TH(i))
        
        BackDataICU(:,th_wave) = [ICU(end-7:end);SimICU(:,i)];
        th_wave = th_wave + 1;
        if th_wave > length(TH_index)
            th_wave = length(TH_index);
        end
    else
        yyaxis left
        plot([ICU;SimICU(2:end,i)], ':','LineWidth',0.3)
        colororder('default')
        hold on
        
        yyaxis right
        plot([N;NPath(2:end,i)]/7, ':','LineWidth',0.3)
        colororder('default')
        %         yline(TH(i),'-',TH(i))
        %         colororder('default')
    end
    hold on
end

xline(Tdata,'LineWidth',1.5,'HandleVisibility','off');
xline(71,'LineWidth',1.5,'HandleVisibility','off');

yyaxis left
yline(ICU_limit,'-g',num2str(ICU_limit),'LineWidth',2.5,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
yline(ICU_limit/2,'-k',num2str(floor(ICU_limit/2)),'LineWidth',2.5,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
% yline(th_off1,'-r',num2str(th_off1),'LineWidth',1.5,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
% yline(th_off2,'-b',num2str(th_off2),'LineWidth',1.5,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');
% yline(th_off3,'-k',num2str(th_off3),'LineWidth',1.5,'HandleVisibility','off','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom');


ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'b';
ax.YAxis(1).FontSize = fs;
ax.YAxis(2).FontSize = fs;
ax.XAxis.FontSize = fs;


%ytickformat('%,6.0f')
xticks(find(WeekNumber==1))
if l == 1
    title('ICU and New cases','FontSize',fs,'FontWeight','normal')
    xticklabels(MonthWeekEN(WeekNumber==1))
    %xticklabels(datestr(MonthWeekEN(xtick1), 'mmm-yy'))
elseif l == 2
    title('重症患者数と新規感染者数の推移','FontSize',fs,'FontWeight','normal','FontName',fn)
    xticklabels(MonthWeekJP(WeekNumber==1))
end
% lgd = legend;
lgd.NumColumns = 2;
xtickangle(45)
xlim([Tdata-7 Tdata+36])

