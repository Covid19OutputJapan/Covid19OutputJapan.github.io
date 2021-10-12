% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

%To estimate RA(1) coefficient of error term, I used Cochrane-Orcutt estimation method.
%https://en.wikipedia.org/wiki/Cochrane%E2%80%93Orcutt_estimation
%This code is consistent until December 2021.
%All the important regression results are stored in the struct var "Result".

clear variables
%close all
iPC=0;
if iPC==1
    %home = '\Users\shcor\Dropbox\Website\Codes\';
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
    %home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/CodesGDO_prefecture/';
    %home = '/Users/asaihiroyuki/covid/';
    home = '/Users/asaihiroyuki/Dropbox/fujii_nakata/monthly_gdp_prefectures/GDO_prefecture/';
end

%cd(home);
%====================== Program parameter values ======================%
figure_save = 1;    % 0 = figures won't be saved, 1 = they will be saved
% in the "Figure" folder
fs = 16;            % common font size for many figures
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
%======================================================================%


close all



%PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka'};
PrefVector ={"Aichi", 'Akita', 'Aomori', 'Chiba', 'Ehime','Fukui', 'Fukuoka', 'Fukushima', 'Gifu', 'Gunma', 'Hiroshima','Hokkaido', 'Hyogo', 'Ibaraki', 'Ishikawa', 'Iwate','Kagawa', 'Kagoshima', 'Kanagawa', 'Kochi', 'Kumamoto', 'Kyoto','Mie', 'Miyagi', 'Miyazaki', 'Nagano', 'Nagasaki', 'Nara','Niigata', 'Oita', 'Okayama', 'Okinawa', 'Osaka', 'Saga','Saitama', 'Shiga', 'Shimane', 'Shizuoka', 'Tochigi', 'Tokushima','Tokyo', 'Tottori', 'Toyama', 'Wakayama', 'Yamagata', 'Yamaguchi','Yamanashi'};
%GDOVector = [106,36,23,21,40];

%Pref_VarVectors = {};structの形で、県名、covid_monthlyでの順番、変数名の３つを含む形にする


GDO_preds_star_simple=[]; %Strage for predicted GDO

plot_mode = 0;%If 1, some statistics are plotted in output glaph.

Result = struct();

for pindex = 41:41%:length(PrefVector)
    
    %%%%%ここを地域ごとに変える
    VarVector = {'Fujii','Malt','Cred','Credserv','Pos','Trav','Event','Rest','Job','Eventinfo','Travelinfo','Movieinfo','TDB','Lag'};
    %Pref_VarVector = Pref_VarVectors;%ここにインポートしてきた使用変数名ファイルを読み込む
    %VarVector = VarVector + Pref_VarVector;
    
    
    size_var = size(VarVector,2);

    VarVector_lag = horzcat(VarVector,{'Lag'});

    var_plot='Malt';


    VarVectorlag = {'Lag'};
    
    
    

    %====================== Model parameter values ======================%
    pref = PrefVector{pindex}      % prefecture to be analyzed
    %prefGDO = GDOVector(pindex);

    if iPC==1
        covid = importdata([home '\Covid_weekly.csv']);  % Import weekly Covid data by prefecture
    else
        covid = importdata([home 'Covid_monthly.csv']);  % Import weekly Covid data by prefecture
    end
    Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
    Tdata= size(Data,1); 
    xtick1 = 1:Tdata;


    year = 21; 
    month = Tdata-12;
    xticklabel1 = horzcat([2001:2012],[2101:year*100+month;]);


    %%データを変数に格納
    for varindex = 1:size_var-size(VarVectorlag,2) %lag,job_pro以外をimport
        eval(strcat(VarVector{varindex},'_raw= Data(:,',num2str(6+varindex),');'));
    end

    %TDB_raw = reallog(TDB_raw); %TDBはアンケートデータなのでlogをとる

    for varindex = 1:size_var-size(VarVectorlag,2)
        eval(strcat(VarVector{varindex},'= normalize(',VarVector{varindex},'_raw);'));
    end
    
    %{
    for varindex = pref_varindex %県にspecificなデータ系列を読み込む
        eval(strcat(VarVector{varindex},'_raw= Data(:,',num2str(6+varindex),');'));
        eval(strcat(VarVector{varindex},'= normalize(',VarVector{varindex},'_raw);'));
    end
    %}
    
    Depart = Data(:,29); 
    Conveni =  Data(:,30);  
    Kaden = Data(:,31);
    Drugstore = Data(:,32);
    
    Depart = normalize(Depart);
    Conveni = normalize(Conveni);
    Kaden = normalize(Kaden);
    Drugstore = normalize(Drugstore);


    GDO = Data(:,6);

    GDO_pred = GDO; 
    GDO_pred_star = GDO;
    TdataGDO = Tdata-sum(isnan(GDO));

    Lag = nan(Tdata,1);
    Lag(2:TdataGDO+1) = GDO(1:TdataGDO);

    %変数を変える時は以下を変える。(Data_cellとVarVector_useから落とす）
    Data_cell = {Fujii,Malt,Cred,Credserv,Pos,Trav,Event,Rest,Job,Eventinfo,Travelinfo,Movieinfo,TDB,Depart,Conveni,Lag};
    
    %VarVector_useはVarVectorの部分集合
    %VarVector_use = {'Fujii','Malt','Cred','Credserv','Pos','Trav','Event','Rest','Job','Eventinfo','Travelinfo','Movieinfo','TDB','Depart','Conveni','Kaden','Drugstore','Lag'};
    VarVector_use = {'Fujii','Malt','Cred','Credserv','Pos','Trav','Event','Rest','Job','Eventinfo','Travelinfo','Movieinfo','TDB','Depart','Conveni','Lag'};
    Data_struct = cell2struct(Data_cell,VarVector_use,2);
    %%%%%%%%%ここまで弄ればOK
    %以下はData_structベースで動くので、ここが変わってれば下は全部O'K

    VarVector_empty ={};

    for var_name = VarVector_use


        if sum(isnan(getfield(Data_struct,var_name{1}))) == Tdata
            VarVector_empty = horzcat(VarVector_empty,var_name);
            Data_struct = rmfield(Data_struct,var_name);
        end
    end

    VarRemain_store = setdiff(setdiff(VarVector_use,VarVector_empty),{'Lag'});%全てemptyの変数を落とす 


    VarVector_cell = {};
    VarRemain_cell = {};


    for j=1:Tdata-TdataGDO

        %tempの初期化
        VarRemain_temp = {};
        VarVector_temp = {};

        for varindex = 1:size(VarRemain_store,2)

            data_var = getfield(Data_struct,VarRemain_store{varindex});
            if isnan(data_var(12+month+1-j))
                VarRemain_temp = horzcat(VarRemain_temp,VarRemain_store{varindex});
            else
                VarVector_temp = horzcat(VarVector_temp,VarRemain_store{varindex});
            end
        end

        VarRemain_store = VarRemain_temp;

        if j == 1
            VarVector_temp = horzcat(VarVector_temp,'Lag');
        end

        VarVector_cell{j} = VarVector_temp;
        VarRemain_cell{j} = VarRemain_temp;


    end

    VarVector_cell = flip(VarVector_cell);



    str = {'Var name : coeff'};
    str_var  = '';
    Xs_name_cell ={};
    ss_star_cell = {};
    Rs_star_cell = {};
    rhos_star_cell = {};

    for j=1:Tdata-TdataGDO


        str  =horzcat(str,strcat('~+ ',num2str(j),' period data~'));

        Xs_name_j ={};
        Xs_j ={};

        for i = j:Tdata-TdataGDO
            Xs_name_j = horzcat(Xs_name_j,VarVector_cell{i});
        end

        for varindex = 1:size(Xs_name_j,2)
            Xs_j = horzcat(Xs_j,getfield(Data_struct,Xs_name_j{varindex}));
        end



        GDO_preds_j = ones(1,size(Xs_name_j,2)); 
        ss_j = ones(2,size(Xs_name_j,2));
        Rs_j = ones(1,size(Xs_name_j,2));

        GDO_preds_star_j = ones(1,size(Xs_name_j,2));
        ss_star_j = ones(2,size(Xs_name_j,2));
        Rs_star_j = ones(1,size(Xs_name_j,2));
        rhos_star_j = ones(1,size(Xs_name_j,2));

        for varindex=1:size(Xs_name_j,2)
            X_raw_j = Xs_j{varindex};
            X_j = X_raw_j(1:TdataGDO+j-1);

            if isnan(X_j(1))

                %derive rho

                X_j = X_j(2:TdataGDO+j-1);
                Y_j = GDO_pred(2:TdataGDO+j-1);%X_j,Y_jのsizeはTdataGDO-1

                XC_j = [ones(length(X_j),1), X_j];
                s_j = (XC_j'*XC_j)\XC_j'*Y_j; 
                reg_j = XC_j*s_j;
                r_p_j = Y_j - reg_j; %residual
                SSE_p_j = sum(r_p_j.^2);
                SSY_j = sum((Y_j - mean(Y_j)).^2);
                R_j = 1-SSE_p_j/SSY_j;

                GDO_pred_j =[1, X_raw_j(TdataGDO+1)]*s_j;
                GDO_preds_j(varindex) = GDO_pred_j;

                ss_j(:,varindex) = s_j;
                Rs_j(:,varindex) = R_j;

                res = r_p_j;
                res_lag = res(1:TdataGDO-2);

                s_res = (res_lag'*res_lag)\res_lag'*res(2:TdataGDO-1);
                rhos_star_j(varindex) = s_res;

                %derive coeff

                Y_star_j = Y_j(2:TdataGDO-1 +j-1) - s_res*Y_j(1:TdataGDO-2 +j-1);
                X_star_j = X_j(2:TdataGDO-1 +j-1) - s_res*X_j(1:TdataGDO-2 +j-1);

                XC_star_j = [ones(length(X_star_j),1), X_star_j];
                s_star_j = (XC_star_j'*XC_star_j)\XC_star_j'*Y_star_j; 
                reg_star_j = XC_star_j*s_star_j;
                r_p_star_j = Y_star_j - reg_star_j; %residual
                SSE_p_star_j = sum(r_p_star_j.^2);
                SSY_star_j = sum((Y_star_j - mean(Y_star_j)).^2);
                R_star_j = 1-SSE_p_star_j/SSY_star_j;

                GDO_pred_star_j =[1, X_raw_j(TdataGDO+1 +j-1)-s_res*X_raw_j(TdataGDO +j-1)]*s_star_j+s_res*Y_j(TdataGDO-1 +j-1);
                GDO_preds_star_j(varindex) = GDO_pred_star_j;

                ss_star_j(:,varindex) = s_star_j;
                Rs_star_j(:,varindex) = R_star_j;

            else

                Y_j = GDO_pred(1:TdataGDO +j-1);
                XC_j = [ones(length(X_j),1), X_j];
                s_j = (XC_j'*XC_j)\XC_j'*Y_j; 
                reg_j = XC_j*s_j;
                r_p_j = Y_j - reg_j; %residual
                SSE_p_j = sum(r_p_j.^2);
                SSY_j = sum((Y_j - mean(Y_j)).^2);
                R_j = 1-SSE_p_j/SSY_j;

                GDO_pred_j =[1, X_raw_j(TdataGDO+1 +j-1)]*s_j;
                GDO_preds_j(varindex) = GDO_pred_j;

                ss_j(:,varindex) = s_j;
                Rs_j(:,varindex) = R_j;


                res = r_p_j;
                res_lag = res(1:TdataGDO-1 +j-1);

                s_res = (res_lag'*res_lag)\res_lag'*res(2:TdataGDO +j-1);
                rhos_star_j(varindex) = s_res;

                %derive coeff

                Y_star_j = Y_j(2:TdataGDO +j-1) - s_res*Y_j(1:TdataGDO-1 +j-1);
                X_star_j = X_j(2:TdataGDO +j-1) - s_res*X_j(1:TdataGDO-1 +j-1);

                XC_star_j = [ones(length(X_star_j),1), X_star_j];
                s_star_j = (XC_star_j'*XC_star_j)\XC_star_j'*Y_star_j; 
                reg_star_j = XC_star_j*s_star_j;
                r_p_star_j = Y_star_j - reg_star_j; %residual
                SSE_p_star_j = sum(r_p_star_j.^2);
                SSY_star_j = sum((Y_star_j - mean(Y_star_j)).^2);
                R_star_j = 1-SSE_p_star_j/SSY_star_j;

                GDO_pred_star_j =[1, X_raw_j(TdataGDO+1 +j-1)-s_res*X_raw_j(TdataGDO +j-1)]*s_star_j+s_res*Y_j(TdataGDO +j-1);
                GDO_preds_star_j(varindex) = GDO_pred_star_j;

                ss_star_j(:,varindex) = s_star_j;
                Rs_star_j(:,varindex) = R_star_j;

            end

            GDO_pred_star(TdataGDO+1 +j-1) = mean(GDO_preds_star_j);
            GDO_pred = GDO_pred_star;

            %LAgの更新
            Data_struct.Lag(2:TdataGDO+2 +j-1) = GDO_pred(1:TdataGDO+1 +j-1);
            %plot用の説明文作成

            str_var = strcat(' ',Xs_name_j{varindex},':',num2str(round(ss_star_j(2,varindex),2)));
            str = horzcat(str,str_var);


        end
        str = horzcat(str,' ' );


        GDO_pred = GDO_pred_star;

        Xs_name_cell{j} = Xs_name_j;
        Rs_star_cell{j} = Rs_star_j;
        ss_star_cell{j} = ss_star_j;
        rhos_star_cell{j} = rhos_star_j;

        if j==1
            coef_table = table(ss_star_j(2,:)','RowNames',Xs_name_j);
            R_table = table(Rs_star_j','RowNames',Xs_name_j);
        end
        
        IIP = [99.8	99.5	95.8	86.4	78.7	80.2	87.2	88.1	91.5	95.2	94.7	93.8	97.8	95.6]
        ITA = [101.6	101.4	97.4	89	86.7	94.3	94.6	95.4	97.3	98.1	98.1	97.7	96.7	97]
        ITA_mob = [102.7	102.1	94.9	81.8	77.2	86.8	88.1	88.2	91.4	91.2	91.4	90.8	89.5	90.3]%運輸業・郵便業 968.8/10000
        ITA_oroshi = [95.8	95.7	95.1	88.6	82.1	86.5	88.8	89.2	91.7	94.3	92.5	94	94.8	92] %卸売り 1350/10000
        
        JCER = [100	100.6510024	96.69076979	89.69919226	88.31737659	94.65092244	94.71066259	95.35102107	96.99090547	99.1233873	98.77579007	97.23384391	96.50281264]
        RDEI_tokyo = [100	97.29449	92.67529	83.33447	79.09724	87.17791	86.63431	89.2946	90.69135	92.43358	93.75739	92.21624]
        
        Average = 1/2*(RDEI_tokyo' + Fujii_raw(1:12))



    end 


    GDO_pred_star_simple = GDO_pred_star;
    GDO_preds_star_simple = [GDO_preds_star_simple,GDO_pred_star_simple];


    %if pindex == 4 | pindex == 19 |pindex == 32 | pindex == 33 |pindex == 35 |pindex == 41
    if pindex == 41


        %subplot(3,3,pindex)
        figure(pindex)

         %--- Plot mobility data ---%
        %figure(2)
        
        
        %{
        if pindex ~= 48
            yyaxis left
            eval(strcat('plot(',var_plot,');'))%',''k'',''LineWidth'',1.5'';')) %1つ目
            ylabel('Var Index (Normalized)')
            hold on 
        end
        %}
        
        %PLotする変数をいじる時は、ここをいじる！！
        %yyaxis right
        %yyaxis left
        %plot(IIP,'k:+','LineWidth',1.5) %２つ目'g-.',
        %hold on
        %plot(ITA,'k:o','LineWidth',1.5) %3つ目'k-.',
        %hold on
        %plot(ITA_mob,'k:*','LineWidth',1.5) %４つ目'b-.',
        %hold on
        %plot(ITA_oroshi,'k:s','LineWidth',1.5) %４つ目'b-.',
        
        
        %plot(GDO_pred,'r-.','LineWidth',1.5) %２つ目
        %hold on
        plot(GDO,'k-','LineWidth',1.5) %3つ目
        hold on
        plot(GDO_pred,'k--','LineWidth',1.5) %3つ目
        
        %plot(GDO_pred_weight,'g-.','LineWidth',1.5)%表示しないかも
        hold on
        %plot(GDO,'b-.','LineWidth',1.5) %４つ目
        ylabel('Index')
        %xlim([1 Tdata])
        xticks(xtick1)

        xticklabels(xticklabel1)

        %legend('Malt(L axis)','Forecasted GDO simple (R axis)','Forecasted GDO weight (R axis)', 'Data GDO (R axis)','Data Fujii (R axis)','FontSize',13,'Location','southeast');
        %legend('Pos(L axis)','Forecasted GDO naive (R axis)','Forecasted GDO star(R axis)', 'Data Fujii (R axis)','Data GDO (R axis)','FontSize',13,'Location','southeast');
        
        
        %{
        text = '(L axis)'',''Forecasted GDO AR1 (R axis)'', ''Data Fujii (R axis)'',''Data GDO (R axis)'',''FontSize'',13,''Location'',''southeast'');';
        text2 = strcat('legend(',var_plot,text);
        eval(strcat('legend(''',var_plot,text))
        %}
        
        %%plotする変数を変える時は、ここもいじる！！！！
        %legend('IIP','ITA','ITA Transportation','ITA Oroshi','Production-side output ','FontSize',13,'Location','southeast');
        legend('Fujii-Nakata GDP','Nowcasted Fujii-Nakata GDP','FontSize',13,'Location','southeast');
        
        xlabel('year(last two digits)+month','FontSize',10,'FontName',fn)
        title(strcat(pref),'FontSize',10,'FontWeight','normal')
        ax = gca;
        %ax.YAxis(1).Color = 'k';
        %ax.YAxis(2).Color = 'b';
        %ax.YAxis(1).FontSize = fs;
        %ax.YAxis(2).FontSize = fs;
        ax.XAxis.FontSize = fs;
        xtickangle(45)



        dim = [.13 .27 .35 .65];

        %{
        VarVector1_withR = {};
        for i =1:size1
            var_with_R = strcat(VarVector1{i},':',num2str(round(Rs_star_1(i),2)),', ',num2str(round(ss_star_1(2,i),2)),', ',num2str(round(rhos_star_1(i),2)));
            VarVector1_withR = horzcat(VarVector1_withR,var_with_R);
        end

        VarVector2_withR = {};
        for i =1:size2
            var_with_R = strcat(VarVector2{i},':',num2str(round(Rs_star_2(i),2)),', ',num2str(round(ss_star_2(2,i),2)),', ',num2str(round(rhos_star_1(i),2)));
            VarVector2_withR = horzcat(VarVector2_withR,var_with_R);
        end

        VarVector3_withR = {};
        for i =1:size3
            var_with_R = strcat(VarVector3{i},':',num2str(round(Rs_star_3(i),2)),', ',num2str(round(ss_star_3(2,i),2)),', ',num2str(round(rhos_star_1(i),2)));
            VarVector3_withR = horzcat(VarVector3_withR,var_with_R);
        end
        %}
        %str = horzcat({'Var:R^2,coeff,rho',' ','~Until Jan~'},VarVector1_withR,{' '},{'~Until Feb~'},VarVector2_withR,{' '},{'~Until Mar~'},VarVector3_withR);

        plot_mode = 1; %使用した変数の説明を加えるかどうか。


        strs = {};
        for t = 1:Tdata-TdataGDO
            strs = horzcat(strs,strcat('+',num2str(t),'month data'));
        end
        
        plot_mode = 0;

        if plot_mode == 1
            %strs = {'~Until Jan~','~Until Feb~','~Until Mar~','~Until Apr~',};%RDEIが新しくなるまではこれでOK

            str_plot = {'Var names : coeff, R^2'};
            for t = 1:Tdata-TdataGDO
                str_plot = horzcat(str_plot,strs{t});
                if t ~= Tdata-TdataGDO
                    vars_t = setdiff(Xs_name_cell{1,t},Xs_name_cell{1,t+1});
                    for var_index = 1:size(vars_t,2)
                        var_coeff = coef_table{vars_t(var_index),'Var1'};
                        var_R = R_table{vars_t(var_index),'Var1'};
                        str_var = strcat(vars_t(var_index),' : ',num2str(round(var_coeff,2)),', ',num2str(round(var_R,2)));
                        str_plot = horzcat(str_plot,str_var);

                    end

                else
                    vars_t = Xs_name_cell{1,t};
                    for var_index = 1:size(vars_t,2)
                        var_coeff = coef_table{vars_t(var_index),'Var1'};
                        var_R = R_table{vars_t(var_index),'Var1'};
                        str_var = strcat(vars_t(var_index),' : ',num2str(round(var_coeff,2)),', ',num2str(round(var_R,2)));
                        str_plot = horzcat(str_plot,str_var);

                    end

                end
            str_plot = horzcat(str_plot,' ');
            end
          annotation('textbox',dim,'String',str_plot,'Edgecolor','none');%'FitBoxToText','on')
        end


        if figure_save ==1
            saveas(gcf,strcat(pref,'-AR-AUTO'),'fig');
            saveas(gcf,strcat(pref,'-AR-AUTO'),'png');
        end
    end


    Result(pindex).pref_name = pref;
    Result(pindex).var_name = Xs_name_cell;
    Result(pindex).R = Rs_star_cell;
    Result(pindex).coeff = ss_star_cell;
    Result(pindex).rho = rhos_star_cell;

    
end
%close all
%{
share = [0.071768677
0.006345257
0.00791276
0.037588689
0.009171117
0.005918902
0.035046134
0.014360385
0.013835352
0.015975174
0.020997916
0.03460255
0.037983855
0.024591009
0.008327455
0.008283249
0.006849074
0.009802724
0.063380856
0.004326541
0.01079133
0.019232711
0.014651634
0.016854026
0.006701261
0.015033527
0.008148816
0.006580397
0.016017821
0.008031656
0.013914263
0.007860925
0.07135939
0.005245057
0.041727656
0.011634847
0.004403961
0.030768923
0.016297328
0.005621999
0.189196432
0.00337771
0.008163665
0.006185557
0.007598903
0.011420981
0.00611151
];
%}

%{
Prefecture_JP = ["愛知県","秋田県","青森県","千葉県","愛媛県","福井県","福岡県","福島県","岐阜県","群馬県","広島県","北海道","兵庫県","茨城県","石川県","岩手県","香川県","鹿児島県","神奈川県","高知県","熊本県","京都府","三重県","宮城県","宮崎県","長野県","長崎県","奈良県","新潟県","大分県","岡山県","沖縄県","大阪府","佐賀県","埼玉県","滋賀県","島根県","静岡県","栃木県","徳島県","東京都","鳥取県","富山県","和歌山県","山形県","山口県","山梨県"];

share_data = readtable('gdpshare_fujii.xls');

share = [];
for pindex = 1:47
    for p =1:47
        if strcmp(Prefecture_JP(pindex),share_data{p,1})
            share(pindex) = share_data{p,2};
        end
    end
end
share = share';



for j = 1:size(GDO_preds_star_simple,1)
    GDO_preds_star_simple(j,48) = GDO_preds_star_simple(j,1:47)*share
end

%Japan GDOのplot

plot(GDO_preds_star_simple(:,48),'g-.','LineWidth',1.5) %２つ目
hold on
plot(GDO_preds_star_simple(:,41),'b-.','LineWidth',1.5) %４つ目
hold on
plot(GDO_preds_star_simple(:,32),'c-.','LineWidth',1.5) %４つ目
ylabel('GDO')
%xlim([1 Tdata])
xticks(xtick1)
xticklabels(xticklabel1)

%legend('Malt(L axis)','Forecasted GDO simple (R axis)','Forecasted GDO weight (R axis)', 'Data GDO (R axis)','Data Fujii (R axis)','FontSize',13,'Location','southeast');
legend('Japan','Tokyo','Osaka','FontSize',13,'Location','northwest');
text = '(L axis)'',''Forecasted GDO AR1 (R axis)'', ''Data Fujii (R axis)'',''Data GDO (R axis)'',''FontSize'',13,''Location'',''southeast'');';
%text2 = strcat('legend(',var_plot,text);
%eval(strcat('legend(''',var_plot,text))
xlabel('year(last two digits)+month','FontSize',10,'FontName',fn)
title('Japan & Tokyo & Osaka','FontSize',10,'FontWeight','normal')

if figure_save ==1
    saveas(gcf,strcat("Japan",'-AR-AUTO'),'png');
    saveas(gcf,strcat("Japan",'-AR-AUTO'),'fig');
end





%for share

GDO_preds_star_simple;
yearlabel = horzcat(repelem([2020],12),repelem([2021],month));
monthlabel = horzcat([1:12],[1:month]);
ValCell = num2cell([yearlabel', monthlabel', GDO_preds_star_simple]);
Test = horzcat({'year', 'month'},horzcat(PrefVector,'Japan'));
TESTtable = cell2table([Test; ValCell]);
writetable(TESTtable, 'for_share.csv', 'WriteVariableNames', false);% ヘッダーを書かない設定
    
 
if iPC==1
    gop = importdata([home 'for_share.csv']);  % Import weekly Covid data by prefecture
else
    %gop = importdata([home 'for_share.csv']);  % Import weekly Covid data by prefecture
    gop = importdata([home 'for_share.csv']);
end

Data = gop.data;

Prefecture_EN = ["Aichi", "Akita", "Aomori", "Chiba", "Ehime","Fukui", "Fukuoka", "Fukushima", "Gifu", "Gunma", "Hiroshima","Hokkaido", "Hyogo", "Ibaraki", "Ishikawa", "Iwate","Kagawa", "Kagoshima", "Kanagawa", "Kochi", "Kumamoto", "Kyoto","Mie", "Miyagi", "Miyazaki", "Nagano", "Nagasaki", "Nara","Niigata", "Oita", "Okayama", "Okinawa", "Osaka", "Saga","Saitama", "Shiga", "Shimane", "Shizuoka", "Tochigi","Tokushima","Tokyo", "Tottori", "Toyama", "Wakayama", "Yamagata", "Yamaguchi","Yamanashi","Japan"];
Prefecture_JP = ["愛知県","秋田県","青森県","千葉県","愛媛県","福井県","福岡県","福島県","岐阜県","群馬県","広島県","北海道","兵庫県","茨城県","石川県","岩手県","香川県","鹿児島県","神奈川県","高知県","熊本県","京都府","三重県","宮城県","宮崎県","長野県","長崎県","奈良県","新潟県","大分県","岡山県","沖縄県","大阪府","佐賀県","埼玉県","滋賀県","島根県","静岡県","栃木県","徳島県","東京都","鳥取県","富山県","和歌山県","山形県","山口県","山梨県","日本"];

year = [];
month = [];
GOP = [];

for i_pref = 1:1:length(Prefecture_EN)
    year = vertcat(year,Data(:,1));
    month = vertcat(month,Data(:,2));
    GOP = vertcat(GOP,Data(:,i_pref+2));
end
pref_length = length(Data(:,1));
for i_pref = 1:1:length(Prefecture_EN)
    prefecture(1+pref_length*(i_pref-1):pref_length*(i_pref),1) = Prefecture_EN(i_pref);
    prefectureJP(1+pref_length*(i_pref-1):pref_length*(i_pref),1) = Prefecture_JP(i_pref);
end

prefGOP(1,:) = ["year","month","prefecture","Prefecture_JP","gop"];
prefGOP(2:size(year,1)+1,:) = [year,month,prefecture,prefectureJP,GOP];

%writematrix(prefGOP, [home 'prefecture_gop.xls']);

    %}
   