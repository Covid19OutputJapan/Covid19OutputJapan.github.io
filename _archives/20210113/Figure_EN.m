clear all;
close all;

load Main_Japan_all.mat;

%%% Time series figures (variables and parameters) %%%

Past = [N,GDP,ERN,zeros(Tdata,1)];
VarList2 = ["N","GDP","ERN","V"];
VarName21 = ["Newly Infected","Output","Effective Reproduction","Newly Vaccinated"];
VarName22 = ["Persons",       "",      "Number",                "Persons"];
for i = 1:length(AlphaIndex)
    alphaT = flip(0:0.01*AlphaIndex(i)*2/(SimPeriod-1):(0.01*AlphaIndex(i)*2))';
    %     alphaT = 0.01*AlphaIndex(i)*ones(SimPeriod,1);
    if hconstant == 0
        ERNT = (AlphaIndexVariables(:,1,i)/POP0).*(((1 - h*alphaT).^k).*betaT)./(gammaT+deltaT);
    elseif hconstant == 1
        ERNT = (AlphaIndexVariables(:,1,i)/POP0).*(((1 - 0.01*h(1) - h(2)*alphaT).^k).*betaT)./(gammaT+deltaT);
    end
    ProjectedData = [AlphaIndexVariables(:,5,i),AlphaIndexVariables(:,6,i),ERNT,V];
    CombinedData = [Past;ProjectedData];
    figure(5)
    for j = 1:length(VarList2)
        subplot(2,2,j)
        plot(CombinedData(:,j),'k','LineWidth',2)
        xlim([1 Tdata+SimPeriod])
        title({VarName21(j);VarName22(j)},'FontWeight','normal')
        xline(Tdata);
        xticks([1 27 53 79 98])
        xticklabels({'Jan-20','Jul-20','Jan-21','Jul-21','Dec-21'})
        ax = gca;
        ax.XAxis.FontSize = 8; %10
        xtickangle(45)
        if j == 3
        ylim([0 3]);
        yline(1);
        end
        hold on
    end
    
%     alpha2 = [alpha;alphaT];
%     beta2 = [beta;betaT];
%     beta_tildeT = ((1 - h*alphaT).^k).*betaT;
%     beta_tilde2 = [beta_tilde;beta_tildeT];
%     ERN2 = [ERN;ERNT];
%     delta2 = [delta;deltaT];
%     V2 = [zeros(Tdata,1);V];
%     figure(101)
%     %set(gcf,'Position',[100,100,800,500])
%     for j = 1:length(ParamList2)
%         subplot(2,2,j)
%         plot(eval(ParamList2(j)),'k','LineWidth',2)
%         xlim([1 Tdata+SimPeriod])
%         title(ParamName(j),'Interpreter','latex','FontSize',11,'FontWeight','normal')
%         xline(Tdata);
%         xticks([1 27 53 79 98])
%         xticklabels( {'Jan-20','Jul-20','Jan-21','Jul-21','Dec-21'} )
%         xtickangle(45)
%         hold on
%     end
end

%--- Trade-off figure (lag) ---%
%figure;
figure(8)
plot(100*AverageAlpha,CumD,'k','LineWidth',2,'DisplayName','This week')
%plot(100*AverageAlpha,CumD,'k','LineWidth',2,'DisplayName','baseline')
legend('This week')
hold on
for i = 1:length(wl)
    wlag = wl(i);
    if wlag==1
            plot(100*LagResults(2,:,i),LagResults(1,:,i),'r','LineWidth',2,'DisplayName',sprintf('%.0f week ago',wlag))
        else
            plot(100*LagResults(2,:,i),LagResults(1,:,i),'b','LineWidth',2,'DisplayName',sprintf('%.0f weeks ago',wlag))
    end
%     plot(100*LagResults(2,:,i),LagResults(1,:,i),'LineWidth',2,'DisplayName',sprintf('%.0f weeks ago',wlag))
    hold on
end
hold off
xlabel('Output Loss (%)','FontSize',fs)
ylabel('Cumlative Deaths','FontSize',fs)
xlim(xlim_tradeoff);
ax = gca;
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
lgd = legend;
lgd.NumColumns = 1;
lgd.FontSize = 15;
grid on

xtick1_F=[StartDateF+4:8:Tdata Tdata];
xtick4_F=[StartDateF+4:8:Tdata Tdata];
fs_F=10;
fs_legend_F=10;
fs_title=10;
%fs_title=12;
%figure
figure(10)
subplot(2,2,1)
plot(StartDateF:(Tdata),dNForecast(StartDateF:Tdata,1),'r',StartDateF:(Tdata),dNActual(StartDateF:Tdata,1),'k','LineWidth',1.5)  
title({'Conditional Forecast vs. Actual';'(one-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF (Tdata)])
xticks(xtick1_F)
xticklabels(Month(xtick1_F))
ax = gca;
ax.YAxis.Exponent = 0;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ytickformat('%,6.0f')
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,2)
plot(StartDateF+3:(Tdata),dNForecast(StartDateF:Tdata-3,2),'r',StartDateF+3:(Tdata),dNActual(StartDateF:Tdata-3,2),'k','LineWidth',1.5)  
title({'Conditional Forecast vs. Actual';'(four-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+3 (Tdata)])
xticks(xtick4_F)
xticklabels(Month(xtick4_F))
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,3)
plot(StartDateF:(Tdata),dNForecast(StartDateF:Tdata,1)-dNActual(StartDateF:Tdata,1),'k','LineWidth',1.5)  
title({'Conditional Forecast Errors',' (one-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF (Tdata)])
xticks(xtick1_F)
xticklabels(Month(xtick1_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
subplot(2,2,4)
plot(StartDateF+3:(Tdata),dNForecast(StartDateF:Tdata-3,2)-dNActual(StartDateF:Tdata-3,2),'k','LineWidth',1.5)  
title({'Conditional Forecast Errors',' (four-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+3 (Tdata)])
xticks(xtick4_F)
xticklabels(Month(xtick4_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

%figure;
figure(11)
subplot(2,2,1)
plot(StartDateF:(Tdata),dDForecast(StartDateF:Tdata,1),'r',StartDateF:(Tdata),dDActual(StartDateF:Tdata,1),'k','LineWidth',1.5)  
title({'Conditional Forecast vs. Actual';'(one-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
xlim([StartDateF (Tdata)])
xticks(xtick1_F)
xticklabels(Month(xtick1_F))
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,2)
plot(StartDateF+3:(Tdata),dDForecast(StartDateF:Tdata-3,2),'r',StartDateF+3:(Tdata),dDActual(StartDateF:Tdata-3,2),'k','LineWidth',1.5)  
title({'Conditional Forecast vs. Actual';'(four-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+3 (Tdata)])
xticks(xtick4_F)
xticklabels(Month(xtick4_F))
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
legend('Forecast','Actual','FontSize',fs_legend_F,'Location','Northwest')
subplot(2,2,3)
plot(StartDateF:(Tdata),dDForecast(StartDateF:Tdata,1)-dDActual(StartDateF:Tdata,1),'k','LineWidth',1.5)  
title({'Conditional Forecast Errors',' (one-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF (Tdata)])
xticks(xtick1_F)
xticklabels(Month(xtick1_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')
subplot(2,2,4)
plot(StartDateF+3:(Tdata),dDForecast(StartDateF:Tdata-3,2)-dDActual(StartDateF:Tdata-3,2),'k','LineWidth',1.5)  
title({'Conditional Forecast Errors',' (four-week horizon)'},'FontSize',fs_title,'FontWeight','normal')
xlim([StartDateF+3 (Tdata)])
xticks(xtick4_F)
xticklabels(Month(xtick4_F))
yline(0)
ax = gca;
ax.YAxis.FontSize = fs_F;
ax.XAxis.FontSize = fs_F;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

%--- Trade-off figure with UB (baseline) ---%
%figure;
figure(6)
AreaUB=area(X1UB,Y1UB,0);
set(AreaUB(1),'FaceColor',[1 1 1])
set(AreaUB(2),'FaceColor',[0.9 0.9 0.9])
set(AreaUB(3),'FaceColor',[0.8 0.8 0.8])
set(AreaUB(4),'FaceColor',[0.7 0.7 0.7])
set(AreaUB(5),'FaceColor',[0.65 0.65 0.65])
set(AreaUB(6),'FaceColor',[0.65 0.65 0.65])
set(AreaUB(7),'FaceColor',[0.7 0.7 0.7])
set(AreaUB(8),'FaceColor',[0.8 0.8 0.8])
set(AreaUB(9),'FaceColor',[0.9 0.9 0.9])
hold on
for iUB=1:9
	set(AreaUB(iUB),'LineStyle','none')
end
plot(X1UB,CumD,'k','LineWidth',1.5)
xlabel('Output Loss (%)','FontSize',16)
ylabel('Cumlative Deaths','FontSize',16)
xlim([2,5]);
yline(D(end),'r--','LineWidth',1.5);
grid on
ax = gca;
ax.YAxis.FontSize = 16;
ax.XAxis.FontSize = 16;
ax.YAxis.Exponent = 0;
ytickformat('%,6.0f')

%--- Save all figures ---%
%if figure_save == 1
%     saveas(figure(1),[home 'Figures/NewCases.png']);
%     saveas(figure(2),[home 'Figures/MobilityGDPLine.png']);
%     saveas(figure(3),[home 'Figures/Parameters.png']);
%     saveas(figure(4),[home 'Figures/Variables.png']);
    saveas(figure(5),[home 'Figures/VariablesProjection.png']);
    saveas(figure(6),[home 'Figures/BaselineTradeoffUB.png']);
%     saveas(figure(7),[home 'Figures/Sensitivity.png']);
    saveas(figure(8),[home 'Figures/LaggedTradeoff.png']);
%     saveas(figure(9),[home 'Figures/CF_multiplicative.png']);
    saveas(figure(10),[home 'Figures/ForecastErrorsN.png']);
    saveas(figure(11),[home 'Figures/ForecastErrorsD.png']);
    
%     saveas(figure(100),[home 'Figures/MobilityGDPScatter.png']);
%     saveas(figure(101),[home 'Figures/ParametersProjection.png']);
%     saveas(figure(102),[home 'Figures/BaselineTradeoff.png']);
%end