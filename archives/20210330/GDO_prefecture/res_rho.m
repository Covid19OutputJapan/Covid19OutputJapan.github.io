VarVector_sort = horzcat(VarVector1,VarVector2,VarVector3);

for varindex = 7:7
    close all
    VarVector_sort{varindex};
    GDO_ = GDO;
    
    predicted = [ones(length(Xs_1(:,varindex)),1), Xs_1(:,varindex)]*ss_1(:,varindex);
    residual = GDO_ - predicted;
    
    GDO_(13) = GDO_preds_1(varindex);
    figure(1)
    yyaxis left
    plot(residual)
    yyaxis right
    hold on
    plot(GDO_)
    title(VarVector_sort{varindex});
    
    legend(strcat(VarVector_sort{varindex},'-res'),'GDObythivar');
    
    saveas(gcf,strcat(VarVector_sort{varindex},'_res'),'fig');
end

residual_1 = residual(1:11,1);
residual_2 = residual(2:12,1);
corrcoef(residual_1,residual_2)