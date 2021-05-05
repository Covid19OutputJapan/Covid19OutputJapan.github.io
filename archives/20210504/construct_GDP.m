function [potentialGDP, referenceGDP, alpha] = construct_GDP(GDP,TdataGDP)
potentialGDP = zeros(52*3,1);
%potentialGDP(1) = (100/(1.0122))*(1.0063^(1/12)) %GDP in 2019/12 = 100
potentialGDP(1) = 100;
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

alpha = (1 - GDP(1:TdataGDP)./referenceGDP(1:TdataGDP));   % output loss in percentage