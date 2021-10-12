%To estimate RA(1) coefficient of error term, I used Cochrane-Orcutt estimation method.
%https://en.wikipedia.org/wiki/Cochrane%E2%80%93Orcutt_estimation
%This code is consistent until December 2021.
%All the important regression results are stored in the struct var "Result".

clear variables
clear all
%close all
iPC=0;
if iPC==1
    %home = '\Users\shcor\Dropbox\Website\Codes\';
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
    %home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/CodesGDO_prefecture/';
    %home = '/Users/asaihiroyuki/covid/';
    home = '/Users/asaihiroyuki/Dropbox/fujii_nakata/Website/Codes/GDO_prefecture_update/';
end

cd(home);
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
if iPC==1
        covid = importdata([home '\Covid_weekly.csv']);  % Import weekly Covid data by prefecture
    else
        covid = importdata([home 'Covid_monthly.csv']);  % Import weekly Covid data by prefecture
    end
columns = covid.textdata(1,:);
useful_seq = 1:size(columns,2);

previous = load('previousdata.mat');



close all

PrefVector ={"Aichi", 'Akita', 'Aomori', 'Chiba', 'Ehime','Fukui', 'Fukuoka', 'Fukushima', 'Gifu', 'Gunma', 'Hiroshima','Hokkaido', 'Hyogo', 'Ibaraki', 'Ishikawa', 'Iwate','Kagawa', 'Kagoshima', 'Kanagawa', 'Kochi', 'Kumamoto', 'Kyoto','Mie', 'Miyagi', 'Miyazaki', 'Nagano', 'Nagasaki', 'Nara','Niigata', 'Oita', 'Okayama', 'Okinawa', 'Osaka', 'Saga','Saitama', 'Shiga', 'Shimane', 'Shizuoka', 'Tochigi', 'Tokushima','Tokyo', 'Tottori', 'Toyama', 'Wakayama', 'Yamagata', 'Yamaguchi','Yamanashi'};

%結果を格納する変数を作成

Datas_struct = struct();
GDO_preds_star_simple=[]; %Strage for predicted GDO
Result = struct();



Prefs_Vars = readtable('pref_vars.csv'); %使用する変数を読み込む


k_set =1 ; %2:RDEI予測しない

for pindex = 1:47%:length(PrefVector)    sagaは34
   
    
    for k = 1:2
        %k
        Data_struct = struct();%Datas_structに格納するようの変数、一時的

        pref = PrefVector{pindex}

        Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);



    
        Positive = Data(:,3);
        GDO = Data(:,8);
        Fujii_raw = Data(:,10);
        qRDEI = Data(:,15);
        qe = Data(:,6);
        
        
        if k==1
            GDO = qRDEI;
            Tdata= size(Fujii_raw,1)-sum(isnan(Fujii_raw));
            TdataGDO = size(qRDEI,1)-sum(isnan(qRDEI));
            
        else
            qFujii = Data(:,13);
            if k_set ~= 2
                GDO = (qrdei_pred+qFujii)/2;
            end
            Tdata = size(Fujii_raw,1);
            TdataGDO = Tdata-sum(isnan(GDO));
            
        end

        xtick1 = 1:Tdata;
        year = 21; 
        month = Tdata-12; %2020年分の１２ヶ月を引く
        xticklabel1 = horzcat([2001:2012],[2101:year*100+month;]);




        %GDO = RDEI;
        Datas_struct.(pref).positive = Positive;


        GDO_pred = GDO; 
        GDO_pred_star = GDO;




        Lag = nan(Tdata,1);
        Lag(2:TdataGDO+1) = GDO(1:TdataGDO);





        %使用する変数の取り込み
        Pref_Vars = (Prefs_Vars{:,pref})';
        Pref_Vars = rmmissing(Pref_Vars);
        Pref_Vars_lag = horzcat(Pref_Vars,{'Lag'});
        add_size = size(Pref_Vars,2);

        %あとでcellにデータを入れていくときのための文字列セット
        text1 = "Data_cell = {";
        text1_raw = "Data_cell_raw = {";
        text2 = "";
        text2_raw = "";

        text3 = ",Lag};";

        %pref
        for add_var_index=1:add_size
            %Pref_Vars(add_var_index);
            col_index = strcmp(columns,Pref_Vars(add_var_index)) * useful_seq' - 1; %prefectureが消えてるからそのぶん引かないといけない
            if col_index < 1
                pref
                Pref_Vars(add_var_index)
                %test = strcat('No ' ,Pref_Vars(add_var_index) ' in ',pref)

            end
            %col_index;


            if strcmp(Pref_Vars(add_var_index),"TDBall")
                %eval(strcat(Pref_Vars{add_var_index},"_raw=log(Data(:,",num2str(col_index),"));"))
                eval(strcat(Pref_Vars{add_var_index},"_raw=Data(:,",num2str(col_index),");"))
                eval(strcat(Pref_Vars{add_var_index},"=normalize(",Pref_Vars{add_var_index},"_raw);"))
                text2_raw = strcat(text2_raw,",",Pref_Vars(add_var_index),"_raw");
                text2 = strcat(text2,",",Pref_Vars(add_var_index));
            else

                eval(strcat(Pref_Vars{add_var_index},"_raw=Data(:,",num2str(col_index),");"))
                eval(strcat(Pref_Vars{add_var_index},"=normalize(",Pref_Vars{add_var_index},"_raw);"))
                text2_raw = strcat(text2_raw,",",Pref_Vars(add_var_index),"_raw");
                text2 = strcat(text2,",",Pref_Vars(add_var_index));
            end
        end

        eval(strcat(text1,text2,text3));
        eval(strcat(text1_raw,text2_raw,text3));
        Data_struct =  cell2struct(Data_cell, Pref_Vars_lag,2 );
        Data_struct =  cell2struct(Data_cell_raw, Pref_Vars_lag,2 );
        
        if k ==2
            Datas_struct.(pref).vars = Data_struct;%全体の構造体の方にはall Nanの系列も含まれている。
            Datas_struct.(pref).vars_raw = Data_struct;
        end

        VarVector_use = Pref_Vars_lag; %現状ではCSVを直接いじる構造だけど、こっちでいじる仕組みにするときのために一回受けとく

        VarVector_empty ={};

        for var_name = VarVector_use
            a = isnan(getfield(Data_struct,var_name{1}));
            %if sum(isnan(getfield(Data_struct,var_name{1}))) == Tdata %全部の変数が0だったら落とす
            if a(2) == 1 %var_name{1}になってるのはvar_nameがセルになってるからなだけ
                VarVector_empty = horzcat(VarVector_empty,var_name);
                %var_name
                Data_struct = rmfield(Data_struct,var_name);
            end
        end

        VarRemain_store = setdiff(setdiff(VarVector_use,VarVector_empty),{'Lag'});

        VarVector_cell = {};
        VarRemain_cell = {};


        for j=1:Tdata-TdataGDO

            %j

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


        VarVector_cell = flip(VarVector_cell); %flip:順番を逆にする
        
        if k==2
            Datas_struct.(pref).varname = VarVector_cell;
            Datas_struct.(pref).iip_TDBseizo = cov(Kasaiiip(2:12),TDBseizou(2:12));
        end


        str = {'Var name : coeff'};
        str_var  = '';
        Xs_name_cell ={};
        ss_star_cell = {};
        Rs_star_cell = {};
        rhos_star_cell = {};


        for j=1:Tdata-TdataGDO
            
            asai = j;


            str  =horzcat(str,strcat('~+ ',num2str(j),' period data~'));

            Xs_name_j ={};
            Xs_j ={};
            mat_residuals = []; %residual保存用

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
                varname_j = Xs_name_j{varindex};
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
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% structに保存する
                
                if k==2
                    if j==1
                        Datas_struct.(pref).predict.(varname_j) = GDO;

                        Datas_struct.(pref).predict.(varname_j)(12+j) = GDO_pred_star_j;
                    else
                        Datas_struct.(pref).predict.(varname_j)(12+j)= GDO_pred_star_j;
                    end

                    Datas_struct.(pref).Result.(varname_j).coef = s_star_j;
                    Datas_struct.(pref).Result.(varname_j).R = R_star_j;
                    Datas_struct.(pref).Result.(varname_j).rho = s_res;
                    Datas_struct.(pref).Result.(varname_j).coef_raw = s_j;
                    Datas_struct.(pref).Result.(varname_j).R_raw = R_j;

                    Data_struct.Lag(2:TdataGDO+2 +j-1) = GDO_pred(1:TdataGDO+1 +j-1);
                end
                %plot用の説明文作成

                str_var = strcat(' ',Xs_name_j{varindex},':',num2str(round(ss_star_j(2,varindex),2)));
                str = horzcat(str,str_var);


            end
            
            simple = 0;
            if simple ==1
                GDO_pred_star(TdataGDO+1 +j-1) = mean(GDO_preds_star_j);
            else
                GDO_pred_star(TdataGDO+1 +j-1) = GDO_preds_star_j*((1-Rs_star_j).^(-1))'/sum((1-Rs_star_j).^(-1));
            end

            GDO_pred = GDO_pred_star;
            %LAgの更新
            Data_struct.Lag(2:TdataGDO+2 +j-1) = GDO_pred(1:TdataGDO+1 +j-1);




        str = horzcat(str,' ' );


            GDO_pred = GDO_pred_star;

            Xs_name_cell{j} = Xs_name_j;
            Rs_star_cell{j} = Rs_star_j;
            ss_star_cell{j} = ss_star_j;
            rhos_star_cell{j} = rhos_star_j;
            
            
                
            if k==1
                qrdei_pred = GDO_pred_star;
            elseif k==2
                if j==1
                    coef_table = table(ss_star_j(2,:)','RowNames',Xs_name_j);
                    R_table = table(Rs_star_j','RowNames',Xs_name_j);
                end

                if j==1
                    Datas_struct.(pref).predict.ALL = GDO;
                    Datas_struct.(pref).predict.ALL(12+j) = GDO_pred_star_j;
                else
                    Datas_struct.(pref).predict.ALL(12+j)= GDO_pred_star_j;
                end
                
                GDO_pred_star_simple = GDO_pred_star;
                
            end
           

        end


        
        
        
        
        if k==1
            %clearvars -except k_set previous qrdei_pred pindex columns covid Datas_struct figure_save fn fs GDO_preds_star_simple home iPC Prefs_Vars PrefVector Result useful_seq;
        else
            GDO_preds_star_simple = [GDO_preds_star_simple,GDO_pred_star_simple];
        end

    end
    
    if pindex == 1|pindex == 4 | pindex == 19 |pindex == 32 | pindex == 33 |pindex == 35 |pindex == 41|pindex == 34

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
        yyaxis left
        plot(conveni,'g:','LineWidth',1.5) %２つ目'g-.',
        hold on
        plot(TDBall,'m:','LineWidth',1.5) %3つ目'k-.',
        hold on
        plot(newcar,'c:','LineWidth',1.5) %４つ目'b-.',
        
        
        
        yyaxis right
        %plot(GDO_pred_naive,'r-.','LineWidth',1.5)
        %hold on
        plot(GDO_pred,'r-.','LineWidth',1.5) %２つ目
        hold on
        plot(Fujii_raw,'k-.','LineWidth',1.5) %3つ目
        hold on
        %plot(GDO_pred_weight,'g-.','LineWidth',1.5)%表示しないかも
        hold on
        plot(GDO,'b-','LineWidth',1.5) %４つ目
        hold on
        plot(previous.GDO_preds_star_simple(:,pindex),'g-.','LineWidth',1.5) %４つ目
        hold on
        if k_set ~= 2
            plot(qrdei_pred,'c-.','LineWidth',1.5)
            legend('Conveni(L axis)','TDB(L axis)','newcar(L axis)','Forecasted GDO (R axis)', 'Data Fujii (R axis)','Data GDO (R axis)','the other','qRDEI-pred (R axis)','FontSize',13,'Location','southeast');
        else
            legend('Conveni(L axis)','TDB(L axis)','newcar(L axis)','Forecasted GDO (R axis)', 'Data Fujii (R axis)','Data GDO (R axis)','the other','FontSize',13,'Location','southeast');
        end
        
         %４つ目
        ylabel('GDO')
        %xlim([1 Tdata])
        xticks(xtick1)

        xticklabels(xticklabel1)
        
        
        
        xlabel('year(last two digits)+month','FontSize',10,'FontName',fn)
        title(strcat(pref),'FontSize',10,'FontWeight','normal')
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'b';
        ax.YAxis(1).FontSize = fs;
        ax.YAxis(2).FontSize = fs;
        ax.XAxis.FontSize = fs;
        xtickangle(45)



        dim = [.13 .27 .35 .65];
        

     
        plot_mode = 1; %使用した変数の説明を加えるかどうか。


        strs = {};
        for t = 1:Tdata-TdataGDO
            strs = horzcat(strs,strcat('+',num2str(t),'month data'));
        end

        plot_mode = 1;
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
            saveas(gcf,strcat(pref,''),'fig');
            saveas(gcf,strcat(pref,'-'),'png');
        end
        
    end


    Result(pindex).pref_name = pref;
    Result(pindex).var_name = Xs_name_cell;
    Result(pindex).R = Rs_star_cell;
    Result(pindex).coeff = ss_star_cell;
    Result(pindex).rho = rhos_star_cell;
    
    close all

    
end
    
    



%{
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Data_struct = struct();%Datas_structに格納するようの変数、一時的
    
    pref = PrefVector{pindex}; 
    
    Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
    Tdata= size(Data,1); 
    xtick1 = 1:Tdata;
    
    year = 21; 
    month = Tdata-12; %2020年分の12月を引く
    xticklabel1 = horzcat([2001:2012],[2101:year*100+month;]);
    
    GDO = Data(:,7);
    Positive = Data(:,3);
    RDEI = Data(:,10);
    
    %GDO = RDEI;
    Datas_struct.(pref).positive = Positive;
    

    GDO_pred = GDO; 
    GDO_pred_star = GDO;
    TdataGDO = Tdata-sum(isnan(GDO));
    Fujii_raw = Data(:,9);

    Lag = nan(Tdata,1);
    Lag(2:TdataGDO+1) = GDO(1:TdataGDO);
    
    
    
     
    %使用する変数の取り込み
    Pref_Vars = (Prefs_Vars{:,pref})';
    Pref_Vars = rmmissing(Pref_Vars);
    Pref_Vars_lag = horzcat(Pref_Vars,{'Lag'});
    add_size = size(Pref_Vars,2);
    
    %あとでcellにデータを入れていくときのための文字列セット
    text1 = "Data_cell = {";
    text1_raw = "Data_cell_raw = {";
    text2 = "";
    text2_raw = "";

    text3 = ",Lag};";
    
    %pref
    for add_var_index=1:add_size
        %Pref_Vars(add_var_index);
        col_index = strcmp(columns,Pref_Vars(add_var_index)) * useful_seq' - 1; %prefectureが消えてるからそのぶん引かないといけない
        if col_index < 1
            pref
            Pref_Vars(add_var_index)
            %test = strcat('No ' ,Pref_Vars(add_var_index) ' in ',pref)
 
        end
        %col_index;
        
        
        if strcmp(Pref_Vars(add_var_index),"TDBall")
            %eval(strcat(Pref_Vars{add_var_index},"_raw=log(Data(:,",num2str(col_index),"));"))
            eval(strcat(Pref_Vars{add_var_index},"_raw=Data(:,",num2str(col_index),");"))
            eval(strcat(Pref_Vars{add_var_index},"=normalize(",Pref_Vars{add_var_index},"_raw);"))
            text2_raw = strcat(text2_raw,",",Pref_Vars(add_var_index),"_raw");
            text2 = strcat(text2,",",Pref_Vars(add_var_index));
        else
        
            eval(strcat(Pref_Vars{add_var_index},"_raw=Data(:,",num2str(col_index),");"))
            eval(strcat(Pref_Vars{add_var_index},"=normalize(",Pref_Vars{add_var_index},"_raw);"))
            text2_raw = strcat(text2_raw,",",Pref_Vars(add_var_index),"_raw");
            text2 = strcat(text2,",",Pref_Vars(add_var_index));
        end
    end
    
    eval(strcat(text1,text2,text3));
    eval(strcat(text1_raw,text2_raw,text3));
    Data_struct =  cell2struct(Data_cell, Pref_Vars_lag,2 );
    Data_struct =  cell2struct(Data_cell_raw, Pref_Vars_lag,2 );
    
    Datas_struct.(pref).vars = Data_struct;%全体の構造体の方にはall Nanの系列も含まれている。
    Datas_struct.(pref).vars_raw = Data_struct;
    
    VarVector_use = Pref_Vars_lag; %現状ではCSVを直接いじる構造だけど、こっちでいじる仕組みにするときのために一回受けとく
    
    VarVector_empty ={};

    for var_name = VarVector_use
        a = isnan(getfield(Data_struct,var_name{1}));
        %if sum(isnan(getfield(Data_struct,var_name{1}))) == Tdata %全部の変数が0だったら落とす
        if a(2) == 1 %var_name{1}になってるのはvar_nameがセルになってるからなだけ
            VarVector_empty = horzcat(VarVector_empty,var_name);
            %var_name
            Data_struct = rmfield(Data_struct,var_name);
        end
    end

    VarRemain_store = setdiff(setdiff(VarVector_use,VarVector_empty),{'Lag'});
    
    VarVector_cell = {};
    VarRemain_cell = {};


    for j=1:Tdata-TdataGDO
        
        %j
        
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


    VarVector_cell = flip(VarVector_cell); %flip:順番を逆にする
    
    Datas_struct.(pref).varname = VarVector_cell;
    
    %pref
    Datas_struct.(pref).iip_TDBseizo = cov(Kasaiiip(2:12),TDBseizou(2:12));
    %cov(Kasaiiip(2:12),TDBseizou(2:12))
    

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
        mat_residuals = []; %residual保存用

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
        
        %Xs_name_j

        for varindex=1:size(Xs_name_j,2)
            varname_j = Xs_name_j{varindex};
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% structに保存する
            if j==1
                Datas_struct.(pref).predict.(varname_j) = GDO;
                
                Datas_struct.(pref).predict.(varname_j)(12+j) = GDO_pred_star_j;
            else
                Datas_struct.(pref).predict.(varname_j)(12+j)= GDO_pred_star_j;
            end
            
            Datas_struct.(pref).Result.(varname_j).coef = s_star_j;
            Datas_struct.(pref).Result.(varname_j).R = R_star_j;
            Datas_struct.(pref).Result.(varname_j).rho = s_res;
            Datas_struct.(pref).Result.(varname_j).coef_raw = s_j;
            Datas_struct.(pref).Result.(varname_j).R_raw = R_j;
            
            %%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%

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
        
        if j==1
            Datas_struct.(pref).predict.ALL = GDO;
            Datas_struct.(pref).predict.ALL(12+j) = GDO_pred_star_j;
        else
            Datas_struct.(pref).predict.ALL(12+j)= GDO_pred_star_j;
        end



    end


    GDO_pred_star_simple = GDO_pred_star;
    GDO_preds_star_simple = [GDO_preds_star_simple,GDO_pred_star_simple];
    
    if isnan(GDO_pred)'*[1:17]' > 0
        pref
    end

    
    
%}
    
   



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

JCER =[100	100.6510024	96.993739	89.9099266	88.09356445	94.45083198	94.26919577	95.54064928	96.993739	98.44682872	98.26519251	98.08355629	96.81210278	97.17537522	98.99173737];
RDEI_tokyo = [100	97.29449	92.67529	83.33447	79.09724	87.17791	86.63431	89.2946	90.69135	92.43358	93.75739	92.21624];
for j = 1:size(GDO_preds_star_simple,1)
    GDO_preds_star_simple(j,48) = GDO_preds_star_simple(j,1:47)*share;
end

gdo_for_jcer = GDO_preds_star_simple((1:size(JCER,2)),48);
coefs = JCER./gdo_for_jcer;



%Japan GDOのplot

%S = load("GDO_preds_star_simple_3.mat");

plot(GDO_preds_star_simple(:,48),'k-','LineWidth',1.5) %２つ目
hold on
plot(GDO_preds_star_simple(:,41),'b-.','LineWidth',1.5) %４つ目
hold on
plot(JCER,'r-.','LineWidth',1.5) %４つ目
hold on
plot(qe,'g-.','LineWidth',1.5) %４つ目
hold on
plot(previous.GDO_preds_star_simple(:,48),'c-.','LineWidth',1.5) %４つ目
hold on
%plot(S.GDO_preds_star_simple(:,48),'k:','LineWidth',1.5) %４つ目
hold on
%plot(S.GDO_preds_star_simple(:,41),'c:','LineWidth',1.5) %４つ目

ylabel('GDO')

%xlim([1 Tdata])
xticks(xtick1)
xticklabels(xticklabel1)

%legend('Malt(L axis)','Forecasted GDO simple (R axis)','Forecasted GDO weight (R axis)', 'Data GDO (R axis)','Data Fujii (R axis)','FontSize',13,'Location','southeast');
%legend('Japan','Tokyo','JCER','Japan-car-retail','Tokyo-car-retail','FontSize',13,'Location','northwest');
legend('Japan','Tokyo','JCER','QE','The other','FontSize',13,'Location','northwest');
text = '(L axis)'',''Forecasted GDO AR1 (R axis)'', ''Data Fujii (R axis)'',''Data GDO (R axis)'',''FontSize'',13,''Location'',''southeast'');';
%text2 = strcat('legend(',var_plot,text);
%eval(strcat('legend(''',var_plot,text))
xlabel('year(last two digits)+month','FontSize',10,'FontName',fn)
title('Japan & Tokyo','FontSize',10,'FontWeight','normal')

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

writematrix(prefGOP, [home 'prefecture_gop.xls']);


%save('previousdata');


