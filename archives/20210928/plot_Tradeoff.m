function plot_Tradeoff(AlphaM,DM,waves,TH,TH_index,l,linecolor,fs,fn)

plot(AlphaM(waves==0),DM(waves==0),'-go','LineWidth',2,'MarkerSize',10);
hold on
plot(AlphaM(waves==1),DM(waves==1),'-mo','LineWidth',2,'MarkerSize',10);
hold on
plot(AlphaM(waves==2),DM(waves==2),'-bo','LineWidth',2,'MarkerSize',10);
hold on
plot(AlphaM(waves==3),DM(waves==3),'-ko','LineWidth',2,'MarkerSize',10);
% text(AlphaM,DM,string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
hold on
for i = 1:length(TH_index)
    scatter(AlphaM(abs(TH - TH_index(i)) < 0.0001),DM(abs(TH - TH_index(i)) < 0.0001),250,linecolor{i},'filled'); %for loop にする
    hold on
end
if l == 1
    xlabel('Output Loss (hundred million yen)','FontSize',20)
    ylabel('Cumulative Deaths','FontSize',20)
    title('Relationship between Covid-19 and output','FontSize',20,'FontWeight','normal')
elseif l == 2
    xlabel('経済損失 (億円)','FontSize',20,'FontName',fn)
    ylabel('累計死亡者数','FontSize',20,'FontName',fn)
    title('コロナ感染と経済の関係','FontSize',20,'FontWeight','normal','FontName',fn)
end
xlim([0,inf])
xtickangle(45)
grid on
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
ax.YAxis.Exponent = 0;
ax.XAxis.Exponent = 0;
ytickformat('%,6.0f')