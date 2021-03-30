% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

%コクラン=オーカット推定量を用いた


clear variables
%close all
iPC=0;
if iPC==1
    %home = '\Users\shcor\Dropbox\Website\Codes\';
    home = '\Users\masam\Dropbox\fujii_nakata\Website\Codes\';
else
    %home = '/Users/Daisuke/Desktop/Dropbox/Research/fujii_nakata/Website/CodesGDO_prefecture/';
    %home = '/Users/asaihiroyuki/covid/';
    home = '/Users/asaihiroyuki/Dropbox/fujii_nakata/Website/Codes/GDO_prefecture/';
end
cd(home);
%====================== Program parameter values ======================%
figure_save = 1;    % 0 = figures won't be saved, 1 = they will be saved
data_save = 0;      % save back data
% in the "Figure" folder
fs = 16;            % common font size for many figures
if iPC == 1
    fn = 'Yu Gothic';     % Font style for xaxis, yaxis, title
else
    fn = 'YuGothic';
end
%======================================================================%






PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka'};
GDOVector = [106,36,23,21,40];

VarVector = {'Fujii','Malt','Cred','Credserv','Pos','Trav','Event','Rest','Job','Eventinfo','Travelinfo','Movieinfo','Lag'};

size_var = size(VarVector,2);

VarVector_lag = horzcat(VarVector,{'Lag'});

var_plot='Pos';


VarVectorlag = {'Lag'};




GDO_preds_simple=[];
GDO_preds_star_simple=[];

predmode =2;
mode  = 1;
for pindex = 1:5%:length(PrefVector)
    
    %close all

        %====================== Model parameter values ======================%
        pref = PrefVector{pindex};        % prefecture to be analyzed
        prefGDO = GDOVector(pindex);

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
 
        
    
        for varindex = 1:size_var-size(VarVectorlag,2) %lag,job_pro以外をimport
            eval(strcat(VarVector{varindex},'_raw= Data(:,',num2str(6+varindex),');'));
            eval(strcat(VarVector{varindex},'= normalize(',VarVector{varindex},'_raw);'));
        end
        
        GDO = Data(:,6);

        GDO_pred = GDO; 
        GDO_pred_star = GDO;
        TdataGDO = Tdata-sum(isnan(GDO));

        Lag = nan(Tdata,1);
        Lag(2:TdataGDO+1) = GDO(1:TdataGDO);
        
        Data_cell = {Fujii,Malt,Cred,Credserv,Pos,Trav,Event,Rest,Job,Eventinfo,Travelinfo,Movieinfo,Lag};
        Data_struct = cell2struct(Data_cell,VarVector,2);
        
        
     
        
        VarRemain_cell = {};
        
        VarRemain_store = setdiff(VarVector,{'Lag'});
        
        
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
        
        VarVector1 = VarVector_cell{1};
        VarVector2 = VarVector_cell{2};
        VarVector3 = VarVector_cell{3};
        
        size1 = size(VarVector1,2);
        size2 = size(VarVector2,2);
        size3 = size(VarVector3,2);
        
        
    
    for j=1:Tdata-TdataGDO
        
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
        rhos_star_1 = ones(1,size(Xs_name_j,2));
        
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
        
        end
    
    %{
        %1st for Janu

        Xs_1_name = horzcat(VarVector1, VarVector2, VarVector3);
        Xs_1 = [];
        for varindex = 1:size(Xs_1_name,2)
            eval(strcat('Xs_1 = [Xs_1,',Xs_1_name{varindex},'];'));
        end
        
        correl  = corrcoef([GDO(2:TdataGDO),Xs_1(2:TdataGDO,:)]);

        GDO_preds_1 = ones(1,size(Xs_1_name,2)); 
        ss_1 = ones(2,size(Xs_1_name,2));
        Rs_1 = ones(1,size(Xs_1_name,2));
        
        
        GDO_preds_star_1 = ones(1,size(Xs_1_name,2));
        ss_star_1 = ones(2,size(Xs_1_name,2));
        Rs_star_1 = ones(1,size(Xs_1_name,2));
        rhos_star_1 = ones(1,size(Xs_1_name,2));
        

        for varindex=1:size(Xs_1_name,2)

            X_raw_1 = Xs_1(:,varindex);
            X_1 = X_raw_1(1:TdataGDO);

            if isnan(X_1(1))
                
                %derive rho

                X_1 = X_1(2:TdataGDO);
                Y_1 = GDO(2:TdataGDO);%X_1,Y_1のsizeはTdataGDO-1

                XC_1 = [ones(length(X_1),1), X_1];
                s_1 = (XC_1'*XC_1)\XC_1'*Y_1; 
                reg_1 = XC_1*s_1;
                r_p_1 = Y_1 - reg_1; %residual
                SSE_p_1 = sum(r_p_1.^2);
                SSY_1 = sum((Y_1 - mean(Y_1)).^2);
                R_1 = 1-SSE_p_1/SSY_1;
                
                GDO_pred_1 =[1, X_raw_1(TdataGDO+1)]*s_1;
                GDO_preds_1(varindex) = GDO_pred_1;
                
                ss_1(:,varindex) = s_1;
                Rs_1(:,varindex) = R_1;
                
                
                res = r_p_1;
                res_lag = res(1:TdataGDO-2);

                s_res = (res_lag'*res_lag)\res_lag'*res(2:TdataGDO-1);
                rhos_star_1(varindex) = s_res;

                %derive coeff

                Y_star_1 = Y_1(2:TdataGDO-1) - s_res*Y_1(1:TdataGDO-2);
                X_star_1 = X_1(2:TdataGDO-1) - s_res*X_1(1:TdataGDO-2);

                XC_star_1 = [ones(length(X_star_1),1), X_star_1];
                s_star_1 = (XC_star_1'*XC_star_1)\XC_star_1'*Y_star_1; 
                reg_star_1 = XC_star_1*s_star_1;
                r_p_star_1 = Y_star_1 - reg_star_1; %residual
                SSE_p_star_1 = sum(r_p_star_1.^2);
                SSY_star_1 = sum((Y_star_1 - mean(Y_star_1)).^2);
                R_star_1 = 1-SSE_p_star_1/SSY_star_1;

                GDO_pred_star_1 =[1, X_raw_1(TdataGDO+1)-s_res*X_raw_1(TdataGDO)]*s_star_1+s_res*Y_1(TdataGDO-1);
                GDO_preds_star_1(varindex) = GDO_pred_star_1;

                ss_star_1(:,varindex) = s_star_1;
                Rs_star_1(:,varindex) = R_star_1;
                

                
            else

                Y_1 = GDO(1:TdataGDO);
                XC_1 = [ones(length(X_1),1), X_1];
                s_1 = (XC_1'*XC_1)\XC_1'*Y_1; 
                reg_1 = XC_1*s_1;
                r_p_1 = Y_1 - reg_1; %residual
                SSE_p_1 = sum(r_p_1.^2);
                SSY_1 = sum((Y_1 - mean(Y_1)).^2);
                R_1 = 1-SSE_p_1/SSY_1;
                
                GDO_pred_1 =[1, X_raw_1(TdataGDO+1)]*s_1;
                GDO_preds_1(varindex) = GDO_pred_1;
                
                ss_1(:,varindex) = s_1;
                Rs_1(:,varindex) = R_1;
                
                
                res = r_p_1;
                res_lag = res(1:TdataGDO-1);

                s_res = (res_lag'*res_lag)\res_lag'*res(2:TdataGDO);
                rhos_star_1(varindex) = s_res;

                %derive coeff

                Y_star_1 = Y_1(2:TdataGDO) - s_res*Y_1(1:TdataGDO-1);
                X_star_1 = X_1(2:TdataGDO) - s_res*X_1(1:TdataGDO-1);

                XC_star_1 = [ones(length(X_star_1),1), X_star_1];
                s_star_1 = (XC_star_1'*XC_star_1)\XC_star_1'*Y_star_1; 
                reg_star_1 = XC_star_1*s_star_1;
                r_p_star_1 = Y_star_1 - reg_star_1; %residual
                SSE_p_star_1 = sum(r_p_star_1.^2);
                SSY_star_1 = sum((Y_star_1 - mean(Y_star_1)).^2);
                R_star_1 = 1-SSE_p_star_1/SSY_star_1;

                GDO_pred_star_1 =[1, X_raw_1(TdataGDO+1)-s_res*X_raw_1(TdataGDO)]*s_star_1+s_res*Y_1(TdataGDO);
                GDO_preds_star_1(varindex) = GDO_pred_star_1;

                ss_star_1(:,varindex) = s_star_1;
                Rs_star_1(:,varindex) = R_star_1;

            end
            
            GDO_pred_star(TdataGDO+1) = mean(GDO_preds_star_1);
            GDO_pred = GDO_pred_star;

            
        end

  
        GDO_pred(TdataGDO+1) = mean(GDO_preds_1);
            if predmode ==2
                GDO_pred_star(TdataGDO+1) = mean(GDO_preds_star_1);
                GDO_pred = GDO_pred_star;
            end
        
        Lag(2:TdataGDO+2) = GDO_pred(1:TdataGDO+1);
 

        %2nd

        Xs_2_name = horzcat(VarVector2, VarVector3);
        Xs_2 = [];
        for varindex = 1:size(Xs_2_name,2)
            eval(strcat('Xs_2 = [Xs_2,',Xs_2_name{varindex},'];'));
        end

        GDO_preds_2 = ones(1,size(Xs_2_name,2)); 
        ss_2 = ones(2,size(Xs_2_name,2));
        Rs_2 = ones(1,size(Xs_2_name,2));
        
        if predmode ==2
            GDO_preds_star_2 = ones(1,size(Xs_2_name,2)); 
            ss_star_2 = ones(2,size(Xs_2_name,2));
            Rs_star_2 = ones(1,size(Xs_2_name,2));
            rhos_star_2 = ones(1,size(Xs_2_name,2));
        end

        for varindex=1:size(Xs_2_name,2)

            X_raw_2 = Xs_2(:,varindex);
            X_2 = X_raw_2(1:TdataGDO+1);

            if isnan(X_2(1))

                X_2 = X_2(2:TdataGDO+1);%サイズはTdataGDO
                Y_2 = GDO_pred(2:TdataGDO+1);

                XC_2 = [ones(length(X_2),1), X_2];
                s_2 = (XC_2'*XC_2)\XC_2'*Y_2; 
                reg_2 = XC_2*s_2;
                r_p_2 = Y_2 - reg_2; %residual
                SSE_p_2 = sum(r_p_2.^2);
                SSY_2 = sum((Y_2 - mean(Y_2)).^2);
                R_2 = 1-SSE_p_2/SSY_2;

                GDO_pred_2 =[1, X_raw_2(TdataGDO+2)]*s_2;

                GDO_preds_2(varindex) = GDO_pred_2;
                ss_2(:,varindex) = s_2;
                Rs_2(:,varindex) = R_2;
                
                if predmode ==2 
                
                    res = r_p_2;
                    res_lag = res(1:TdataGDO-1);

                    s_res = (res_lag'*res_lag)\res_lag'*res(2:TdataGDO);
                    rhos_star_2(varindex) = s_res;

                    %derive coeff

                    Y_star_2 = Y_2(2:TdataGDO) - s_res*Y_2(1:TdataGDO-1);
                    X_star_2 = X_2(2:TdataGDO) - s_res*X_2(1:TdataGDO-1);

                    XC_star_2 = [ones(length(X_star_2),1), X_star_2];
                    s_star_2 = (XC_star_2'*XC_star_2)\XC_star_2'*Y_star_2; 
                    reg_star_2 = XC_star_2*s_star_2;
                    r_p_star_2 = Y_star_2 - reg_star_2; %residual
                    SSE_p_star_2 = sum(r_p_star_2.^2);
                    SSY_star_2 = sum((Y_star_2 - mean(Y_star_2)).^2);
                    R_star_2 = 1-SSE_p_star_2/SSY_star_2;

                    GDO_pred_star_2 =[1, X_raw_2(TdataGDO+2)-s_res*X_raw_2(TdataGDO+1)]*s_star_2+s_res*Y_2(TdataGDO);
                    GDO_preds_star_2(varindex) = GDO_pred_star_2;


                    ss_star_2(:,varindex) = s_star_2;
                    Rs_star_2(:,varindex) = R_star_2;
                end


            else
                
                %derive rho

                Y_2 = GDO_pred(1:TdataGDO+1);
                XC_2 = [ones(length(X_2),1), X_2];
                s_2 = (XC_2'*XC_2)\XC_2'*Y_2; 
                reg_2 = XC_2*s_2;
                r_p_2 = Y_2 - reg_2; %residual
                SSE_p_2 = sum(r_p_2.^2);
                SSY_2 = sum((Y_2 - mean(Y_2)).^2);
                R_2 = 1-SSE_p_2/SSY_2;

                GDO_pred_2 =[1, X_raw_2(TdataGDO+2)]*s_2;

                GDO_preds_2(varindex) = GDO_pred_2;
                ss_2(:,varindex) = s_2;
                Rs_2(:,varindex) = R_2;
                
                if predmode == 2 
                
                    res = r_p_2;
                    res_lag = res(1:TdataGDO);

                    s_res = (res_lag'*res_lag)\res_lag'*res(2:TdataGDO+1);
                    rhos_star_2(varindex) = s_res;
                    %derive coeff

                    Y_star_2 = Y_2(2:TdataGDO+1) - s_res*Y_2(1:TdataGDO);
                    X_star_2 = X_2(2:TdataGDO+1) - s_res*X_2(1:TdataGDO);

                    XC_star_2 = [ones(length(X_star_2),1), X_star_2];
                    s_star_2 = (XC_star_2'*XC_star_2)\XC_star_2'*Y_star_2; 
                    reg_star_2 = XC_star_2*s_star_2;
                    r_p_star_2 = Y_star_2 - reg_star_2; %residual
                    SSE_p_star_2 = sum(r_p_star_2.^2);
                    SSY_star_2 = sum((Y_star_2 - mean(Y_star_2)).^2);
                    R_star_2 = 1-SSE_p_star_2/SSY_star_2;

                    GDO_pred_star_2 =[1, X_raw_2(TdataGDO+2)-s_res*X_raw_2(TdataGDO+1)]*s_star_2+s_res*Y_2(TdataGDO+1);
                    GDO_preds_star_2(varindex) = GDO_pred_star_2;


                    ss_star_2(:,varindex) = s_star_2;
                    Rs_star_2(:,varindex) = R_star_2;
                end

            end

        end

       
        
        if predmode ==2
            GDO_pred_star(TdataGDO+2) = mean(GDO_preds_star_2);
            GDO_pred = GDO_pred_star;
        else
            GDO_pred(TdataGDO+2) = mean(GDO_preds_2);
        end

        Lag(2:TdataGDO+3) = GDO_pred(1:TdataGDO+2);

        %3rd


        Xs_3_name = horzcat(VarVector3);
        Xs_3 = [];
        for varindex = 1:size(Xs_3_name,2)
            eval(strcat('Xs_3 = [Xs_3,',Xs_3_name{varindex},'];'));
        end

        GDO_preds_3 = ones(1,size(Xs_3_name,2)); 
        ss_3 = ones(2,size(Xs_3_name,2));
        Rs_3 = ones(1,size(Xs_3_name,2));
        
        if predmode ==2
            GDO_preds_star_3 = ones(1,size(Xs_3_name,2)); 
            ss_star_3 = ones(2,size(Xs_3_name,2));
            Rs_star_3 = ones(1,size(Xs_3_name,2));
            rhos_star_3 = ones(1,size(Xs_3_name,2));
        end

        for varindex=1:size(Xs_3_name,2)

            X_raw_3 = Xs_3(:,varindex);
            X_3 = X_raw_3(1:TdataGDO+2);

            if isnan(X_3(1))

                X_3 = X_3(2:TdataGDO+2);
                Y_3 = GDO_pred(2:TdataGDO+2); %Y_3のサイズはTdataGDO+1

                XC_3 = [ones(length(X_3),1), X_3];
                s_3 = (XC_3'*XC_3)\XC_3'*Y_3; 
                reg_3 = XC_3*s_3;
                r_p_3 = Y_3 - reg_3; %residual
                SSE_p_3 = sum(r_p_3.^2);
                SSY_3 = sum((Y_3 - mean(Y_3)).^2);
                R_3 = 1-SSE_p_3/SSY_3;

                GDO_pred_3 =[1, X_raw_3(TdataGDO+3)]*s_3;

                GDO_preds_3(varindex) = GDO_pred_3;
                ss_3(:,varindex) = s_3;
                Rs_3(:,varindex) = R_3;
                
                if predmode ==2
                
                    %derive coeff
                    res = r_p_3;%resのサイズはTdataGDO+1
                    res_lag = res(1:TdataGDO);

                    s_res = (res_lag'*res_lag)\res_lag'*res(2:TdataGDO+1);
                    rhos_star_3(varindex) = s_res;

                    Y_star_3 = Y_3(2:TdataGDO+1) - s_res*Y_3(1:TdataGDO); 
                    X_star_3 = X_3(2:TdataGDO+1) - s_res*X_3(1:TdataGDO);

                    XC_star_3 = [ones(length(X_star_3),1), X_star_3];
                    s_star_3 = (XC_star_3'*XC_star_3)\XC_star_3'*Y_star_3; 
                    reg_star_3 = XC_star_3*s_star_3;
                    r_p_star_3 = Y_star_3 - reg_star_3; %residual
                    SSE_p_star_3 = sum(r_p_star_3.^2);
                    SSY_star_3 = sum((Y_star_3 - mean(Y_star_3)).^2);
                    R_star_3 = 1-SSE_p_star_3/SSY_star_3;

                    GDO_pred_star_3 =[1, X_raw_3(TdataGDO+3)-s_res*X_raw_3(TdataGDO+2)]*s_star_3+s_res*Y_3(TdataGDO+1);
                    GDO_preds_star_3(varindex) = GDO_pred_star_3;

                    ss_star_3(:,varindex) = s_star_3;
                    Rs_star_3(:,varindex) = R_star_3;
                end


            else

                Y_3 = GDO_pred(1:TdataGDO+2);
                XC_3 = [ones(length(X_3),1), X_3];
                s_3 = (XC_3'*XC_3)\XC_3'*Y_3; 
                reg_3 = XC_3*s_3;
                r_p_3 = Y_3 - reg_3; %residual
                SSE_p_3 = sum(r_p_3.^2);
                SSY_3 = sum((Y_3 - mean(Y_3)).^2);
                R_3 = 1-SSE_p_3/SSY_3;

                GDO_pred_3 =[1, X_raw_3(TdataGDO+3)]*s_3;

                GDO_preds_3(varindex) = GDO_pred_3;
                ss_3(:,varindex) = s_3;
                Rs_3(:,varindex) = R_3;
                
                if predmode ==2
                
                    %derive coeff
                    res = r_p_3;
                    res_lag = res(1:TdataGDO);

                    s_res = (res_lag'*res_lag)\res_lag'*res(2:TdataGDO+1);
                    rhos_star_3(varindex) = s_res;

                    Y_star_3 = Y_3(2:TdataGDO+2) - s_res*Y_3(1:TdataGDO+1);
                    X_star_3 = X_3(2:TdataGDO+2) - s_res*X_3(1:TdataGDO+1);

                    XC_star_3 = [ones(length(X_star_3),1), X_star_3];
                    s_star_3 = (XC_star_3'*XC_star_3)\XC_star_3'*Y_star_3; 
                    reg_star_3 = XC_star_3*s_star_3;
                    r_p_star_3 = Y_star_3 - reg_star_3; %residual
                    SSE_p_star_3 = sum(r_p_star_3.^2);
                    SSY_star_3 = sum((Y_star_3 - mean(Y_star_3)).^2);
                    R_star_3 = 1-SSE_p_star_3/SSY_star_3;

                    GDO_pred_star_3 =[1, X_raw_3(TdataGDO+3)-s_res*X_raw_3(TdataGDO+2)]*s_star_3+s_res*Y_3(TdataGDO+2);
                    GDO_preds_star_3(varindex) = GDO_pred_star_3;

                    ss_star_3(:,varindex) = s_star_3;
                    Rs_star_3(:,varindex) = R_star_3;
                end

                
            end
           
        end


        
        if predmode ==2
            GDO_pred_star(TdataGDO+3) = mean(GDO_preds_star_3);
            GDO_pred = GDO_pred_star;
        else
            GDO_pred(TdataGDO+3) = mean(GDO_preds_3);
        end

     %}
        GDO_pred = GDO_pred_star;
    end
    
            
    
        
    GDO_pred_star_simple = GDO_pred_star;
    GDO_preds_star_simple = [GDO_preds_star_simple,GDO_pred_star_simple];
   
        
        

    

    %subplot(3,3,pindex)
    figure(pindex)
   
     %--- Plot mobility data ---%
    %figure(2)
    yyaxis left
    eval(strcat('plot(',var_plot,');'))%',''k'',''LineWidth'',1.5'';')) %1つ目
    ylabel('Var Index (Normalized)')
    hold on 
    yyaxis right
    %plot(GDO_pred_naive,'r-.','LineWidth',1.5)
    %hold on
    plot(GDO_pred,'g-.','LineWidth',1.5) %２つ目
    hold on
    plot(Fujii_raw,'k-.','LineWidth',1.5) %3つ目
    hold on
    %plot(GDO_pred_weight,'g-.','LineWidth',1.5)%表示しないかも
    hold on
    plot(GDO,'b-.','LineWidth',1.5) %４つ目
    ylabel('GDO')
    %xlim([1 Tdata])
    xticks(xtick1)

    xticklabels(xticklabel1)

    %legend('Malt(L axis)','Forecasted GDO simple (R axis)','Forecasted GDO weight (R axis)', 'Data GDO (R axis)','Data Fujii (R axis)','FontSize',13,'Location','southeast');
    %legend('Pos(L axis)','Forecasted GDO naive (R axis)','Forecasted GDO star(R axis)', 'Data Fujii (R axis)','Data GDO (R axis)','FontSize',13,'Location','southeast');
    text = '(L axis)'',''Forecasted GDO AR1 (R axis)'', ''Data Fujii (R axis)'',''Data GDO (R axis)'',''FontSize'',13,''Location'',''southeast'');';
    text2 = strcat('legend(',var_plot,text);
    eval(strcat('legend(''',var_plot,text))
    xlabel('year(last two digits)+month','FontSize',10,'FontName',fn)
    title(strcat(pref),'FontSize',10,'FontWeight','normal')
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'b';
    ax.YAxis(1).FontSize = fs;
    ax.YAxis(2).FontSize = fs;
    ax.XAxis.FontSize = fs;
    xtickangle(45)
    
    
    
    dim = [.13 .17 .15 .65];
    
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
    str = horzcat({'Var',' ','~Until Jan~'},VarVector1,{' '},{'~Until Feb~'},VarVector2,{' '},{'~Until Mar~'},VarVector3);
    annotation('textbox',dim,'String',str,'Edgecolor','none');%'FitBoxToText','on')
    
    saveas(gcf,strcat(pref,'-AR-AUTO'),'fig');
    
end
  
    

save('GDO_preds','GDO_preds_simple') %for sharing 
%save('GDO_preds_weight','GDO_preds_weight') 
   