function plot_Tradeoff2(AlphaM,DM,waves,TH,TH_index,l,linecolor,fs,fn)

% text(AlphaM,DM,string(TH),'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',14);
hold on
for i = 1:length(TH_index)
%     scatter(AlphaM(abs(TH - TH_index(i)) < 0.0001),DM(abs(TH - TH_index(i)) < 0.0001),250,linecolor{i},'filled'); %for loop にする
    scatter(AlphaM(i),DM(i),250,linecolor{i},'filled'); %for loop にする
    hold on
end
if l == 1
    xlabel('Output Loss (hundred million yen)','FontSize',fs)
    ylabel('Cumulative Deaths','FontSize',fs)
    title('Relationship between Covid-19 and output','FontSize',fs,'FontWeight','normal')
elseif l == 2
    xlabel('経済損失 (億円)','FontSize',fs,'FontName',fn)
    ylabel('累計死亡者数','FontSize',fs,'FontName',fn)
    title('コロナ感染と経済の関係','FontSize',fs,'FontWeight','normal','FontName',fn)
end
xlim([0,inf])
xtickangle(45)
grid on
ax = gca;
ax.YAxis.FontSize = fs;
ax.XAxis.FontSize = fs;
ax.YAxis.Exponent = 0;
ax.XAxis.Exponent = 0;
ytickformat('%,6.0f')