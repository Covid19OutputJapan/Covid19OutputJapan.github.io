VarVector_sort = horzcat(VarVector1,VarVector2,VarVector3);

for varindex = 1:13
    close all
    VarVector_sort{varindex};
    GDO_ = GDO;
    
    GDO_(13) = GDO_preds_1(varindex)
    figure(1)
    yyaxis left
    plot(Xs_1(:,varindex)*100)
    yyaxis right
    hold on
    plot(GDO_)
    title(VarVector_sort{varindex});
    
    legend(VarVector_sort{varindex},'GDO');
    
    saveas(gcf,strcat(VarVector_sort{varindex}),'fig');
end