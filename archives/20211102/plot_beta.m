function plot_beta(var_initial,var_share,beta,beta_avg,betaT,betaT_woAR1,dateD,SimDate,Tdata,tt,MonthWeek,MonthNumber,WeekNumber,fn)
% plot beta

plot_var_share(1,1) = var_initial;
plot_var_share(2:length(SimDate)+1,1) = var_share(:,1);
tt = 12; %Showing previous t periods for the plot

plot_betaT(1:tt,1) = beta(end-tt+1:end);
plot_betaT(tt+1:length(SimDate)+tt,1) = betaT(:);

plot_betaT_woAR1(1:tt,1) = beta(end-tt+1:end);
plot_betaT_woAR1(tt+1:length(SimDate)+tt,1) = betaT_woAR1(:);

plot_betaT3(1:tt,1) = beta(end-tt+1:end);
plot_betaT3(tt+1:length(SimDate)+tt,1) = beta_avg;

xaxis_vec = 0:1:length(SimDate);
xaxis_vec2 = 1:1:length(SimDate)+tt;
xaxis_vec3 = 1:1:Tdata+tt;

MonthWeek_Sim = MonthWeek(length(dateD)+1:end,1);
MonthNumber_Sim = MonthNumber(length(dateD)+1:end,1);
WeekNumber_Sim = WeekNumber(length(dateD)+1:end,1);

MonthWeek_Sim2 = MonthWeek(length(dateD)-tt+1:end,1);
MonthNumber_Sim2 = MonthNumber(length(dateD)-tt+1:end,1);
WeekNumber_Sim2 = WeekNumber(length(dateD)-tt+1:end,1);

subplot(1,2,1)
plot(xaxis_vec,plot_var_share(:,1)*100,'-r','LineWidth',2)
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
ax.YAxis.Exponent = 0;
ytickformat('%,3.0f')
xticks(find(WeekNumber_Sim==1))
title('新規感染者数の推移','FontSize',20,'FontWeight','normal','FontName',fn)
title('変異株シェアの推移','FontSize',20,'FontWeight','normal','FontName',fn)
xticklabels(MonthWeek_Sim(WeekNumber_Sim==1))
xlabel('週','FontSize',20,'FontName',fn)
ylabel('変異株割合(%)','FontSize',20,'FontName',fn)
xtickangle(45)
xline(1,'--','LineWidth',1.5,'HandleVisibility','off');
xlim([0,length(SimDate)]);

subplot(1,2,2)
plot(xaxis_vec2(tt+1:end), plot_betaT(tt+1:end,1),'-r','LineWidth',2)
hold on
plot(xaxis_vec2(tt+1:end), plot_betaT_woAR1(tt+1:end,1),'-b','LineWidth',2.0)
plot(xaxis_vec2,plot_betaT3,'-k','LineWidth',2)
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
ax.YAxis.Exponent = 0;
ytickformat('%,0.2f')
xticks(find(WeekNumber_Sim2==1))
xtickangle(45)
legend('感染率(変異株成長ケース) W/ AR1','感染率(変異株成長ケース) W/O AR1','過去の平均','FontSize',10,'FontName',fn)
xticklabels(MonthWeek_Sim2(WeekNumber_Sim2==1))
title('βの推移','FontSize',20,'FontWeight','normal','FontName',fn)
xlabel('週','FontSize',20)
ylabel('感染率','FontSize',20,'FontName',fn)
xline(tt+1,'--','LineWidth',1.5,'HandleVisibility','off');
xlim([0,length(SimDate)+tt]);
lgd = legend;
lgd.Location = 'southeast';
