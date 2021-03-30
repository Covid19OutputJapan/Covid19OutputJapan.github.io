% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

clear variables
close all
iPC=0;
if iPC==1
    %home = '\Users\shcor\Dropbox\Website\Codes\';
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
    %home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/Codes';
    %home = '/Users/asaihiroyuki/covid/';
    home = '/Users/asaihiroyuki/Dropbox/fujii_nakata/Website/Codes/GDO_prefecture/';
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 0;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 0;      % save back data
% in the "Figure" folder
fs = 16;            % common font size for many figures
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
%======================================================================%

PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka','Aichi','Fukuoka'};%AichiとFukuokaは外す
GDOVector = [106,36,23,21,40,41,20];

GDO_preds = [];
s_1s ={};

for reg_mode = 3:5
    %1:全て使用 2:過去を使わない 3:Maltのみ 4:Fujiのみ 5:過去のみ 6:


    for pindex = 1:1%:length(PrefVector)
        close all

        %====================== Model parameter values ======================%
        pref = PrefVector{pindex};        % prefecture to be analyzed
        prefGDO = GDOVector(pindex);


        if iPC==1
            covid = importdata([home '\Covid_weekly.csv']);  % Import weekly Covid data by prefecture
        else
            covid = importdata([home 'Covid_monthly.csv']);  % Import weekly Covid data by prefecture
        end
        Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
        % Columns: 1 = date, 2 = new positive, 3 = death, 4 = mobility, 5 = GDO, 6 = population
        dateD = Data(:,1);
        N = Data(:,2);
        dD = Data(:,3);
        M = Data(:,5);
        Fujii =Data(:,7);
        GDO = Data(:,6);
        %ERN_TK = Data(:,7);
        Tdata= size(Data,1);    % Data period in month
        xtick1 = 1:Tdata;
        %xticklabel1 = [202001, 202002, 202003, 202004,202005, 202006, 202007, 202008, 202009,202010, 202011,202012,202101,202102,202103]; %1から始まってTdataで終わる、13ずつインクリメントする目盛
        xticklabel1 = [2001, 2002, 2003, 2004,2005, 2006, 2007, 2008, 2009,2010, 2011,2012,2101,2102,2103];


        M = 1+0.01*M;
        TdataGDO = Tdata-sum(isnan(GDO));
        TdataFujii = Tdata-sum(isnan(Fujii));


        %--- Impute alpha (regress alpha on M)---%

        Malt=M;

        %Malt(50)=0.5*(Malt(49)+Malt(51)); %外れ値を穴埋め？

        %1st for Janu Y(t) = Fujii(t) + Mobility(t) + Y(t-1)

        GDO_pred = GDO; 



        X1_1 = Malt(2:TdataGDO);
        X2_1 = Fujii(2:TdataGDO);
        X3_1 = GDO(1:TdataGDO-1);
        Y_1 = GDO(2:TdataGDO);

        %reg_mode = 1
        %1:全て使用 2:過去を使わない 3:Maltのみ 4:Fujiのみ 5:過去のみ

        if reg_mode == 1
            XC_1 = [ones(length(X1_1),1), X1_1,X2_1,X3_1];
        elseif reg_mode == 2
            XC_1 = [ones(length(X1_1),1), X1_1,X2_1];
        elseif reg_mode == 3
            XC_1 = [ones(length(X1_1),1), X1_1];
        elseif reg_mode == 4
            XC_1 = [ones(length(X1_1),1), X2_1];
        elseif reg_mode == 5
            XC_1 = [ones(length(X1_1),1), X3_1];
        elseif reg_mode == 6
            XC_1 = [ones(length(X1_1),1), X2_1,X3_1];
        end

        s_1 = (XC_1'*XC_1)\XC_1'*Y_1;   
        
        eval(strcat('s_1_',num2str(reg_mode) ,'= s_1;'));
 
        reg_1 = XC_1*s_1;
        r_p_1 = Y_1 - reg_1; %residual
        SSE_p_1 = sum(r_p_1.^2);
        SSY_1 = sum((Y_1 - mean(Y_1)).^2);
        R_1 = 1-SSE_p_1/SSY_1;

        for i = 1:TdataFujii-TdataGDO
            if reg_mode == 1
                GDO_pred(TdataGDO+i) = [1,Malt(TdataGDO+i),Fujii(TdataGDO+i),GDO(TdataGDO+i-1)]*s_1;
            elseif reg_mode == 2
                GDO_pred(TdataGDO+i) = [1,Malt(TdataGDO+i),Fujii(TdataGDO+i)]*s_1;
            elseif reg_mode == 3
                GDO_pred(TdataGDO+i) = [1,Malt(TdataGDO+i)]*s_1;
            elseif reg_mode == 4
                GDO_pred(TdataGDO+i) = [1,Fujii(TdataGDO+i)]*s_1;
            elseif reg_mode == 5
                GDO_pred(TdataGDO+i) = [1,GDO(TdataGDO+i-1)]*s_1;
            elseif reg_mode == 5
                GDO_pred(TdataGDO+i) = [1,Fuji(TdataGDO+i),GDO(TdataGDO+i-1)]*s_1;
            end
        end

        ['mobility Fujii pastGDO GDO']
        correl  = corrcoef([X1_1,X2_1,X3_1,Y_1])


        %eps_p_1 = zeros(Tdata-TdataGDO,1);

        %2nd for Feb Mar Y(t) = Mobility(t) + Y(t-1)
        X1_2 = Malt(2:TdataGDO);
        X2_2 = GDO(1:TdataGDO-1);
        Y_2 = GDO(2:TdataGDO);
        XC_2 = [ones(length(X1_2),1), X1_2, X2_2];
        s_2 = (XC_2'*XC_2)\XC_2'*Y_2;         
        reg_2 = XC_2*s_2;
        r_p_2 = Y_2 - reg_2; %residual
        SSE_p_2 = sum(r_p_2.^2);
        SSY_2 = sum((Y_2 - mean(Y_2)).^2);
        R_2 = 1-SSE_p_2/SSY_2;

        for i = 1:Tdata-TdataFujii
            GDO_pred(TdataFujii+i) = [1,Malt(TdataFujii+i),GDO_pred(TdataFujii+i-1)]*s_2;
        end

        %subplot(3,3,pindex)
        %figure(2)
        yyaxis left
        plot(GDO_pred,'k','LineWidth',1.5)
        ylabel('prdicted')
        yyaxis right
         %--- Plot mobility data ---%
        %figure(2)
        yyaxis left
        plot(Malt,'k','LineWidth',1.5)
        ylabel('Mobility')
        yyaxis right
        hold on
        plot(GDO_pred,'r-.','LineWidth',1.5)
        hold on
        plot(GDO,'b-.','LineWidth',1.5)
        hold on
        plot(Fujii,'k-.','LineWidth',1.5)
        ylabel('GDO')
        %xlim([1 Tdata])
        xticks(xtick1)

        xticklabels(xticklabel1)

        legend('Mobility (L axis)','Forecasted GDO (R axis)', 'Data GDO (R axis)','Data Fujii (R axis)','FontSize',16,'Location','southeast');

        xlabel('year(last two digits)+month','FontSize',10,'FontName',fn)
        title(strcat('mode=',num2str(reg_mode),' Monthly GDO of  ', pref,',R1:',num2str(R_1),',R2:',num2str(R_2)),'FontSize',10,'FontWeight','normal')

        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'b';
        ax.YAxis(1).FontSize = fs;
        ax.YAxis(2).FontSize = fs;
        ax.XAxis.FontSize = fs;
        xtickangle(45)
        saveas(gcf,strcat('mode=',num2str(reg_mode),pref),'fig');

        GDO_preds = [GDO_preds,GDO_pred];


    end
    

end

save('GDO_preds','GDO_preds') %for sharing 