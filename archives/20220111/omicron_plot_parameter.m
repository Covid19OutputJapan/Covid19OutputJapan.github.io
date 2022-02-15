% Transitions of Beta
figure('Name', char(figname_beta));
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of \beta", "\betaの推移"];
yft = '%.3f';
lineNameJP = strings(nX,nZ);
lineNameEN = strings(nX,nZ);
LineStyles = {"-",   "-",   "-"; ...
    "--",  "--", "--"; ...
    ":",   ":",  ":"};

for iX = 1:nX
    for iZ = 1:nZ
        lineNameJP(iX,iZ) = ['相対重症化率 = ', num2str(xvec(iX)), ', シナリオ ... ', char(Scenario(iZ))];
        lineNameEN(iX,iZ) = ['Rel. Inf = ', num2str(xvec(iX)), ', Rel. BRN = ', char(ScenarioEN(iZ))];
    end
end

for iX = 1:nX  
    l = 2;
    if l == 1
        lineName = lineNameEN;
    elseif l == 2
        lineName = lineNameJP;
    end

    plot_4Dfunction(beta, betaPath(:,:,:,:), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,linecolor, LineStyles(iX,:), lineName(iX,:))
end

if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_beta) '.png']);
end

% Transitions of Beta Tilde
figure('Name', char(figname_beta_tilde));
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of \beta tilde", "\beta tildeの推移"];
yft = '%.3f';
for iX = 1:nX
 
    if l == 1
        lineName = lineNameEN;
    elseif l == 2
        lineName = lineNameJP;
    end
    plot_4Dfunction(beta_tilde, betaTildePath(:,:,:,:), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,linecolor, LineStyles(iX,:), lineName(iX,:))
end

if figure_save == 1
    %     saveas(gcf, [home 'Figures/' char(pref) '/' char(figname_beta_tilde) '.png']);
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_beta_tilde) '.png']);
end

% Transitions of ERN
figure('Name', char(figname_ERN));
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of ERN", "実効再生産数の推移"];
for iX = 1:nX

    if l == 1
        lineName = lineNameEN;
    elseif l == 2
        lineName = lineNameJP;
    end

    plot_4Dfunction(ERN, SimERN(:,:,:,:), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,linecolor, LineStyles(iX,:), lineName(iX,:))
end

if figure_save == 1
    %     saveas(gcf, [home 'Figures/' char(pref) '/' char(figname_ERN) '.png']);
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_ERN) '.png']);
end

% Transitions of BRN
BRNpast = beta_tilde ./ (gamma + delta);

figure('Name', char(figname_BRN));
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of BRN", "基本再生産数の推移"];
for iX = 1:nX

    if l == 1
        lineName = lineNameEN;
    elseif l == 2
        lineName = lineNameJP;
    end

    plot_4Dfunction(BRNpast, SimBRN(:,:,:,:), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,linecolor, LineStyles(iX,:), lineName(iX,:))
end

if figure_save == 1
    %     saveas(gcf, [home 'Figures/' char(pref) '/' char(figname_ERN) '.png']);
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_BRN) '.png']);
end
%%
% Transitions of Death Rate
figure('Name', char(figname_delta));
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of Death Rate", "死亡率の推移"];
yft = '%.3f';
for iX = 1:nX
    if l == 1
        lineName = lineNameEN;
    elseif l == 2
        lineName = lineNameJP;
    end

    plot_4Dfunction(delta, deltaPath(:,:,:,2), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,{"k"}, LineStyles(iX,2), lineName(iX,2))
end
if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_delta) '.png']);
end

% Transitions of ICU Rate
figure('Name', char(figname_ICU_nation));
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of ICU Rate (National Standard)", "重症化率の推移(国基準)"];
for iX = 1:nX
    
    plot_4Dfunction(ICU_nation_rate, SimICU_nation_rate_Path(:,:,:,2), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,{"k"}, LineStyles(iX,2), lineName(iX,2))
end
if figure_save == 1
    %     saveas(gcf, [home 'Figures/' char(pref) '/' char(figname_ICU_nation) '.png']);
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_ICU_nation) '.png']);
end

% Transitions of ICU Rate

figure('Name', char(figname_ICU_local));
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of ICU Rate (Local Standard)", "重症化率の推移(都基準)"];
for iX = 1:nX
plot_4Dfunction(ICU_pref_rate, SimICU_pref_rate_Path(:,:,:,2), iX, ...
    WeekNumber, YearMonth, xmin, xmax, ...
    fn, fs, lgdfs, axfs,yft,...
    lgdLocation, column_num, l, title_vec, ...
    lineWidth,{"k"}, LineStyles(iX,2), lineName(iX,2))
end
if figure_save == 1
    %     saveas(gcf, [home 'Figures/' char(pref) '/' char(figname_ICU_local) '.png']);
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_ICU_local) '.png']);
end

% Transitions of Hospital Rate
figure('Name', [char(figname_Hospital) num2str(xvec(iX),'%.2f')]);
    set(gcf, 'Position', [100, 100, 1200, 800])
    title_vec = ["Transitions of Hospital Rate", "入院率の推移"];
for iX = 1:nX

    plot_4Dfunction(Hospital_rate, SimHospital_rate_Path(:,:,:,2), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,{"k"}, LineStyles(iX,2), lineName(iX,2))
end
  %ここも変更
    if figure_save == 1
        saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_Hospital) '.png']);
    end
