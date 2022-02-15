% Transitions of Beta
figure('Name', char(figname_beta));
LineStyles = {"-", "-", "-","-";
              "--", "--", "--","--";
              ":", ":", ":",":"};
color_vec = ["b", "k", "r"];
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of \beta", "\betaの推移"];
yft = '%.3f';
for iX = 1:nX
    lineNameJP = strings(nY,nZ);
    lineNameEN = strings(nY,nZ);
    for iY = 1:nY
        for iZ = 1:nZ
            lineNameJP(iY,iZ) = ['相対感染力 = ', num2str(omicron_relative_infectivity_vector(iX)), ', 相対ワクチン有効性 = ', num2str(omicron_E2_vector(iY))];
            lineNameEN(iY,iZ) = ['Rel. Inf = ', num2str(omicron_relative_infectivity_vector(iX)), ', Rel. VE = ', num2str(omicron_E2_vector(iY))];
            linecolor{iY, iZ} = color_vec(iX);
        end
    end

    l = 2;
    if l == 1
        lineName = lineNameEN;
    elseif l == 2
        lineName = lineNameJP;
    end

    plot_4Dfunction(beta, betaPath(:,:,:,1), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,linecolor, LineStyles, lineName)
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
    lineNameJP = strings(nY,nZ);
    lineNameEN = strings(nY,nZ);
    for iY = 1:nY
        for iZ = 1:nZ
            lineNameJP(iY,iZ) = ['相対感染力 = ', num2str(omicron_relative_infectivity_vector(iX)), ', 相対ワクチン有効性 = ', num2str(omicron_E2_vector(iY))];
            lineNameEN(iY,iZ) = ['Rel. Inf = ', num2str(omicron_relative_infectivity_vector(iX)), ', Rel. VE = ', num2str(omicron_E2_vector(iY))];
            linecolor{iY, iZ} = color_vec(iX);
        end
    end

    if l == 1
        lineName = lineNameEN;
    elseif l == 2
        lineName = lineNameJP;
    end
    plot_4Dfunction(beta_tilde, betaTildePath(:,:,:,1), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,linecolor, LineStyles, lineName)
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

    lineNameJP = strings(nY,nZ);
    lineNameEN = strings(nY,nZ);
    for iY = 1:nY
        for iZ = 1:nZ
            lineNameJP(iY,iZ) = ['相対感染力 = ', num2str(omicron_relative_infectivity_vector(iX)), ', 相対ワクチン有効性 = ', num2str(omicron_E2_vector(iY))];
            lineNameEN(iY,iZ) = ['Rel. Inf = ', num2str(omicron_relative_infectivity_vector(iX)), ', Rel. VE = ', num2str(omicron_E2_vector(iY))];
            linecolor{iY, iZ} = color_vec(iX);
        end
    end

    if l == 1
        lineName = lineNameEN;
    elseif l == 2
        lineName = lineNameJP;
    end

    plot_4Dfunction(ERN, SimERN(:,:,:,1), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,linecolor, LineStyles, lineName)
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

    lineNameJP = strings(nY,nZ);
    lineNameEN = strings(nY,nZ);
    for iY = 1:nY
        for iZ = 1:nZ
            lineNameJP(iY,iZ) = ['相対感染力 = ', num2str(omicron_relative_infectivity_vector(iX)), ', 相対ワクチン有効性 = ', num2str(omicron_E2_vector(iY))];
            lineNameEN(iY,iZ) = ['Rel. Inf = ', num2str(omicron_relative_infectivity_vector(iX)), ', Rel. VE = ', num2str(omicron_E2_vector(iY))];
            linecolor{iY, iZ} = color_vec(iX);
        end
    end

    if l == 1
        lineName = lineNameEN;
    elseif l == 2
        lineName = lineNameJP;
    end

    plot_4Dfunction(BRNpast, SimBRN(:,:,:,1), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,linecolor, LineStyles, lineName)
end

if figure_save == 1
    %     saveas(gcf, [home 'Figures/' char(pref) '/' char(figname_ERN) '.png']);
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_BRN) '.png']);
end

% Transitions of Death Rate
figure('Name', char(figname_delta));
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of Death Rate", "死亡率の推移"];
yft = '%.3f';
iX = 1;
LineStyles = {"-", "--", ":","-.";"-", "--", ":","-.";"-", "--", ":","-."};
lineNameJP = strings(nY,nZ);
lineNameEN = strings(nY,nZ);
for iY = 1:nY
    for iZ = 1:nZ
        lineNameJP(iY,iZ) = ['相対重症化率 = ', num2str(omicron_realtive_severity_vector(iZ))];
        lineNameEN(iY,iZ) = ['Rel. Sev. = ', num2str(omicron_realtive_severity_vector(iZ))];
        linecolor{iY, iZ} = 'r';
    end
end

if l == 1
    lineName = lineNameEN;
elseif l == 2
    lineName = lineNameJP;
end

plot_4Dfunction(delta, deltaPath(:,:,:,:), iX, ...
    WeekNumber, YearMonth, xmin, xmax, ...
    fn, fs, lgdfs, axfs,yft,...
    lgdLocation, column_num, l, title_vec, ...
    lineWidth,linecolor, LineStyles, lineName)

if figure_save == 1
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_delta) '.png']);
end

% Transitions of ICU Rate
figure('Name', char(figname_ICU_nation));
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of ICU Rate (National Standard)", "重症化率の推移(国基準)"];
plot_4Dfunction(delta .* ICU_nation_inflow, ICU_nation_inflow_avg * delta_ICU_nationPath(:,:,:,:), iX, ...
    WeekNumber, YearMonth, xmin, xmax, ...
    fn, fs, lgdfs, axfs,yft,...
    lgdLocation, column_num, l, title_vec, ...
    lineWidth,linecolor, LineStyles, lineName)

if figure_save == 1
    %     saveas(gcf, [home 'Figures/' char(pref) '/' char(figname_ICU_nation) '.png']);
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_ICU_nation) '.png']);
end

% Transitions of ICU Rate

figure('Name', char(figname_ICU_local));
set(gcf, 'Position', [100, 100, 1200, 800])
title_vec = ["Transitions of ICU Rate (Local Standard)", "重症化率の推移(都基準)"];
plot_4Dfunction(delta .* ICU_pref_inflow, ICU_pref_inflow_avg * delta_ICU_prefPath(:,:,:,:), iX, ...
    WeekNumber, YearMonth, xmin, xmax, ...
    fn, fs, lgdfs, axfs,yft,...
    lgdLocation, column_num, l, title_vec, ...
    lineWidth,linecolor, LineStyles, lineName)

if figure_save == 1
    %     saveas(gcf, [home 'Figures/' char(pref) '/' char(figname_ICU_local) '.png']);
    saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_ICU_local) '.png']);
end

% Transitions of Hospital Rate
for iX = 1:nX

    figure('Name', [char(figname_Hospital) num2str(xvec(iX),'%.2f')]);
    set(gcf, 'Position', [100, 100, 1200, 800])
    title_vec = [["Transitions of Hospital Rate", ''], "入院率の推移"];
    % iX = 2;
    % iY = 2;
    plot_4Dfunction(delta .* Hospital_inflow, Hospital_inflow_avg * delta_HospitalPath(:,:,:,:), iX, ...
        WeekNumber, YearMonth, xmin, xmax, ...
        fn, fs, lgdfs, axfs,yft,...
        lgdLocation, column_num, l, title_vec, ...
        lineWidth,linecolor, LineStyles, lineName)
    
    %ここも変更
    if figure_save == 1
        saveas(gcf, [home 'Figures/' char(pref) '/' char(figfolder) '/' char(figname_Hospital) '.png']);
    end

end

