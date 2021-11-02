%--- Import data ---%
% Covid data are recorded at weekly frequencies (Mon-Sun)
% The first week start on January 20 (Mon), 2020
Data = covid.data(strcmp(covid.textdata(2:end,1),pref),:);
date = Data(:,1) + 21916;
date1 =datetime(date,'ConvertFrom','excel','Format','MMM-yy');
date = datetime(date,'ConvertFrom','excel');

dateD = Data(:,1) + 21916;

N = Data(:,2);
dD = Data(:,3);
M = Data(:,4);
GDP = Data(:,5);
Tdata = length(N);
ICU = zeros(Tdata+1,1);
BED = zeros(Tdata,1);%

I_data = zeros(Tdata+1,1);
I_data(2:end) = Data(:,26);
R_data = zeros(Tdata+1,1);
R_data(2:end)= Data(:,27);
D_data = zeros(Tdata+1,1);
D_data(2:end) = Data(:,28);
S_data = ones(Tdata+1, 1) * POP0 - (I_data + R_data + D_data);
dI_data = I_data(2:end) - I_data(1:end-1);
dR_data = R_data(2:end) - R_data(1:end-1);
dD_data = D_data(2:end) - D_data(1:end-1);
dS_data = S_data(2:end) - S_data(1:end-1);

if data_switch == 1
    I = I_data;
    gamma_data = zeros(Tdata, 1);
    %     delta_data = zeros(Tdata,1);
    for i = 1:Tdata
        
        if I(i) > 0
            gamma_data(i) = dR_data(i) / I_data(i);
            %             delta_data(i) = dD_data(i)/I_data(i);
        end
        
    end
    
    gamma = mean(gamma_data(end - RetroPeriod + 1:end));
end

ICU(2:Tdata+1,1) = Data(:,22);
BED(1:Tdata,1) = Data(:,23);
ps = POP0/125710000;
xtick1 = 1:13:Tdata;
dateEN = datetime(date);
SimDate = date(end)+7:7:date(end)+7*SimPeriod;
SimDateEN = datetime(SimDate);
date = [date; SimDate'];

M = 1+0.01*M;
TdataGDP = Tdata-sum(isnan(GDP));
RetroH = TdataGDP-4;
Month = string(date1);

%--- Constructing the reference level of output ---%
potentialGDP = zeros(52*3,1);       % potential GDP for the next 3 years

potentialGDP(1) = (548182/(1.0122))*(1.0063^(1/12));

for i = 2:length(potentialGDP)
    if i <= 13
        potentialGDP(i) = potentialGDP(i-1)*(1.0063^(1/52));
    elseif i <= 52
        potentialGDP(i) = potentialGDP(i-1)*(1.0021^(1/52));
    elseif i <= 104
        potentialGDP(i) = potentialGDP(i-1)*(1.0024^(1/52));
    elseif i <= 156
        potentialGDP(i) = potentialGDP(i-1)*(1.0021^(1/52));
    end
end

referenceGDP = potentialGDP.*(1+0.0166);
referenceGDP(1:2) = [];


%--- Impute alpha (regress alpha on M)---%
Malt=M;
Malt(50)=0.5*(Malt(49)+Malt(51));
alpha = (1 - GDP(1:TdataGDP)./referenceGDP(1:TdataGDP));   % output loss in percentage
X = Malt(TdataGDP-17:TdataGDP);
Y = alpha(TdataGDP-17:TdataGDP);
XC = [ones(length(X),1), X];
s = (XC'*XC)\XC'*Y;         % OLS estimate of h with constant
reg = XC*s;
r = Y - reg;
SSE = sum(r.^2);
eps_p = zeros(Tdata-TdataGDP,1);
eps_p(1) = r(end);
for i = 1:Tdata-TdataGDP-1
    eps_p(i+1) = 1*eps_p(i);
end
alpha_pred = s(1)+s(2)*Malt(TdataGDP+1:Tdata)+eps_p;

alpha = [alpha;alpha_pred];