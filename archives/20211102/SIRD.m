function [S,I,R,D,GDP]...
    = SIRD(Tdata,POP0,N,E1,E2,V1_elderly,V1_medical,V1_others,V2_elderly,V2_medical,V2_others,gamma,dD,TdataGDP,referenceGDP,alpha,GDP)
%plot Simulated New Cases
%--- Compute the history of S, I, R, D in the data period ---%
S = zeros(Tdata+1,1);
I = zeros(Tdata+1,1);
R = zeros(Tdata+1,1);
D = zeros(Tdata+1,1);
V1_elderly = [0;0;V1_elderly(1:end-2)]; % Do not make these vectors as output
V1_medical = [0;0;V1_medical(1:end-2)]; % Do not make these vectors as output
V1_others = [0;0;V1_others(1:end-2)]; % Do not make these vectors as output

V2_elderly = [0;0;V2_elderly(1:end-2)]; % Do not make these vectors as output
V2_medical = [0;0;V2_medical(1:end-2)]; % Do not make these vectors as output
V2_others = [0;0;V2_others(1:end-2)]; % Do not make these vectors as output

S(1)=POP0;
for i = 1:Tdata
    S(i+1)=S(i)-N(i)-E1*(V1_elderly(i)+V1_medical(i)+V1_others(i))-(E2-E1)*(V2_elderly(i)+V2_medical(i)+V2_others(i));
    I(i+1)=I(i)+N(i)-gamma*I(i)-dD(i);
    R(i+1)=R(i)+gamma*I(i)+E1*(V1_elderly(i)+V1_medical(i)+V1_others(i))+(E2-E1)*(V2_elderly(i)+V2_medical(i)+V2_others(i));
    D(i+1)=D(i)+dD(i);
    if i > TdataGDP
        GDP(i) = referenceGDP(i)*(1-alpha(i));
    end
end
