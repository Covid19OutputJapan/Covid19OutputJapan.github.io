% This m-file executes simulation and generates figures for the
% analysis of the state of emergency placed in Tokyo in the paper
% "Covid-19 and Output in Japan" by Daisuke Fujii and Taisuke
% Nakata

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

PrefVector = {'Tokyo','Kanagawa','Saitama','Chiba','Osaka'};%AichiとFukuokaは外す
GDOVector = [106,36,23,21,40];

VarVector = {'Fujii','Malt','Cred','Credserv','Pos','Trav','Event','Rest','Job','Eventinfo','Travelinfo','Movieinfo','Lag','Joblag'};
size_var = size(VarVector,2);

VarVector_lag = horzcat(VarVector,{'Lag'});

VarVector1 = {'Fujii','Cred','Credserv','Trav','Event','Jobpre'};
size1 = size(VarVector1,2);
VarVector2 = {'Rest','Job'};
size2 = size(VarVector2,2);
VarVector3 = {'Malt','Pos','Eventinfo','Travelinfo','Movieinfo','Lag'};
size3 = size(VarVector3,2);
VarVectorlag = {'Lag','Jobpre'};

xticklabel1 = [2001, 2002, 2003, 2004,2005, 2006, 2007, 2008, 2009,2010, 2011,2012,2101,2102,2103];

GDO_preds_weight=[];
GDO_preds_simple=[];

for pindex = 5:5%:length(PrefVector)
    for mode = 1:2
        

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
        Tdata= size(Data,1); 
        xtick1 = 1:Tdata;
        
    
        for varindex = 1:size_var-size(VarVectorlag,2) %lag,job_lag以外をimport
            eval(strcat(VarVector{varindex},'_raw= Data(:,',num2str(6+varindex),');'));
            eval(strcat(VarVector{varindex},'= normalize(',VarVector{varindex},'_raw);'));

        end

        GDO = Data(:,6);

        GDO_pred = GDO; 
        TdataGDO = Tdata-sum(isnan(GDO));

        Lag = nan(Tdata,1);
        Lag(2:TdataGDO+1) = GDO(1:TdataGDO);
        
        Jobpre = nan(Tdata,1);
        Jobpre(1:Tdata-sum(isnan(Job))-1) = Job(2:Tdata-sum(isnan(Job)));

        %1st

        Xs_1_name = horzcat(VarVector1, VarVector2, VarVector3);
        Xs_1 = [];
        for varindex = 1:size(Xs_1_name,2)
            eval(strcat('Xs_1 = [Xs_1,',Xs_1_name{varindex},'];'));
        end
        
        
        correl  = corrcoef([GDO(2:TdataGDO),Xs_1(2:TdataGDO,:)]);

        GDO_preds_1 = ones(1,size(Xs_1_name,2)); 
        ss_1 = ones(2,size(Xs_1_name,2));
        Rs_1 = ones(1,size(Xs_1_name,2));

        for varindex=1:size(Xs_1_name,2)

            X_raw_1 = Xs_1(:,varindex);
            X_1 = X_raw_1(1:TdataGDO);

            if isnan(X_1(1))

                X_1 = X_1(2:TdataGDO);
                Y_1 = GDO(2:TdataGDO);

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

            end

        end

        if mode == 1
            GDO_pred(TdataGDO+1) = mean(GDO_preds_1);
        else
            GDO_pred(TdataGDO+1) = GDO_preds_1*((1-Rs_1).^(-1))'/sum((1-Rs_1).^(-1));
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

        for varindex=1:size(Xs_2_name,2)

            X_raw_2 = Xs_2(:,varindex);
            X_2 = X_raw_2(1:TdataGDO+1);

            if isnan(X_2(1))

                X_2 = X_2(2:TdataGDO+1);
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

            else

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

            end

        end

        if mode == 1
            GDO_pred(TdataGDO+2) = mean(GDO_preds_2);
        else
            GDO_pred(TdataGDO+2) = GDO_preds_2*((1-Rs_2).^(-1))'/sum((1-Rs_2).^(-1));
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

        for varindex=1:size(Xs_3_name,2)

            X_raw_3 = Xs_3(:,varindex);
            X_3 = X_raw_3(1:TdataGDO+2);

            if isnan(X_3(1))

                X_3 = X_3(2:TdataGDO+2);
                Y_3 = GDO_pred(2:TdataGDO+2);

                XC_3 = [ones(length(X_3),1), X_3];
                s_3 = (XC_3'*XC_3)\XC_3'*Y_3; 
                reg_2 = XC_3*s_3;
                r_p_3 = Y_3 - reg_2; %residual
                SSE_p_3 = sum(r_p_3.^2);
                SSY_3 = sum((Y_3 - mean(Y_3)).^2);
                R_3 = 1-SSE_p_3/SSY_3;

                GDO_pred_3 =[1, X_raw_3(TdataGDO+3)]*s_3;

                GDO_preds_3(varindex) = GDO_pred_3;
                ss_3(:,varindex) = s_3;
                Rs_3(:,varindex) = R_3;

            else

                Y_3 = GDO_pred(1:TdataGDO+2);
                XC_3 = [ones(length(X_3),1), X_3];
                s_3 = (XC_3'*XC_3)\XC_3'*Y_3; 
                reg_2 = XC_3*s_3;
                r_p_3 = Y_3 - reg_2; %residual
                SSE_p_3 = sum(r_p_3.^2);
                SSY_3 = sum((Y_3 - mean(Y_3)).^2);
                R_3 = 1-SSE_p_3/SSY_3;

                GDO_pred_3 =[1, X_raw_3(TdataGDO+3)]*s_3;

                GDO_preds_3(varindex) = GDO_pred_3;
                ss_3(:,varindex) = s_3;
                Rs_3(:,varindex) = R_3;

            end

        end

        if mode== 1
            GDO_pred(TdataGDO+3) = mean(GDO_preds_3);
        else
            GDO_pred(TdataGDO+3) = GDO_preds_3*((1-Rs_3).^(-1))'/sum((1-Rs_3).^(-1));
        end
            
        
        if mode ==1
            GDO_pred_simple = GDO_pred;
            GDO_preds_simple = [GDO_preds_simple,GDO_pred_simple];
        else
            GDO_pred_weight = GDO_pred;
            GDO_preds_weight = [GDO_preds_weight,GDO_pred_weight];
        end
        
        
        
    end
    

    %subplot(3,3,pindex)
    figure(pindex)
    
    plot(Fujii,'LineWidth',1.5)
    hold on 
    plot(Cred,'LineWidth',1.5)
    ylabel('Var Index')
    hold on
    plot(Credserv,'LineWidth',1.5)
    hold on
    plot(Trav,'LineWidth',1.5)
    hold on
    plot(Event,'LineWidth',1.5)
    hold on
    plot(Jobpre,'LineWidth',1.5)
    hold on
    plot(Rest,'LineWidth',1.5)
    
    hold on
    plot(Malt,'LineWidth',1.5)
    hold on
    plot(Pos,'LineWidth',1.5)
    hold on
    plot(Eventinfo,'LineWidth',1.5)
    hold on
    plot(Travelinfo,'LineWidth',1.5)
    hold on
    plot(Movieinfo,'LineWidth',1.5)
    hold on
    plot(normalize(Lag),'LineWidth',1.5)
    
    %xlim([1 Tdata])
    xticks(xtick1)

    xticklabels(xticklabel1)
    %legend('Fuji','Cred','Credserv','Trav','Event','Jobpre','Rest','Job','Malt','pos','eventinfo','travelinfo','movieinfo','Lag')
    legend('Fuji','Cred','Credserv','Trav','Event','Jobpre','Rest','Malt','pos','eventinfo','travelinfo','movieinfo','Lag')
    %legend('Jobpre(L axis)','Forecasted GDO simple (R axis)','Forecasted GDO weight (R axis)', 'Data GDO (R axis)','Data Fujii (R axis)','FontSize',13,'Location','southeast');
    %legend('Jobpre(L axis)','Forecasted GDO simple (R axis)', 'Data GDO (R axis)','Data Fujii (R axis)','FontSize',13,'Location','southeast');
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
    
    VarVector1_withR = {};
    for i =1:size1
        var_with_R = strcat(VarVector1{i},':',num2str(round(Rs_1(i),2)),', ',num2str(round(ss_1(2,i),2)));
        VarVector1_withR = horzcat(VarVector1_withR,var_with_R);
    end
    
    VarVector2_withR = {};
    for i =1:size2
        var_with_R = strcat(VarVector2{i},':',num2str(round(Rs_2(i),2)),', ',num2str(round(ss_2(2,i),2)));
        VarVector2_withR = horzcat(VarVector2_withR,var_with_R);
    end
        
    VarVector3_withR = {};
    for i =1:size3
        var_with_R = strcat(VarVector3{i},':',num2str(round(Rs_3(i),2)),', ',num2str(round(ss_3(2,i),2)));
        VarVector3_withR = horzcat(VarVector3_withR,var_with_R);
    end
    str = horzcat({'Var :R^2,coeff',' ','~Until Jan~'},VarVector1_withR,{' '},{'~Until Feb~'},VarVector2_withR,{' '},{'~Until Mar~'},VarVector3_withR);
    %str = horzcat({'Var :R^2,coeff',' ','~Until Jan~'},VarVector1,{' '},{'~Until Feb~'},VarVector2,{' '},{'~Until Mar~'},VarVector3);
    annotation('textbox',dim,'String',str,'Edgecolor','none');%'FitBoxToText','on')
    
    saveas(gcf,strcat(pref,'-allvar'),'png');
   

end


    

save('GDO_preds_simple','GDO_preds_simple') %for sharing 
save('GDO_preds_weight','GDO_preds_weight') 
   