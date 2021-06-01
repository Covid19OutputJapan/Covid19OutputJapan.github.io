function [S,I,R,D]...
    = SIRD(Tdata,POP0,N,E1,E2,V1_elderly,V1_medical,V2_elderly,V2_medical,gamma,dD,TdataGDP,referenceGDP,alpha)
%plot Simulated New Cases
%--- Compute the history of S, I, R, D in the data period ---%
S = zeros(Tdata+1,1);
I = zeros(Tdata+1,1);
R = zeros(Tdata+1,1);
D = zeros(Tdata+1,1);
S(1)=POP0;
for i = 1:Tdata
    S(i+1)=S(i)-N(i)-E1*(V1_elderly(i)+V1_medical(i))-(E2-E1)*(V2_elderly(i)+V2_medical(i));
    I(i+1)=I(i)+N(i)-gamma*I(i)-dD(i);
    R(i+1)=R(i)+gamma*I(i)+E1*(V1_elderly(i)+V1_medical(i))+(E2-E1)*(V2_elderly(i)+V2_medical(i));
    D(i+1)=D(i)+dD(i);
    if i > TdataGDP
        GDP(i) = referenceGDP(i)*(1-alpha(i));
    end
end
