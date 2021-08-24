elderly_total = 100;
E1 = 0.1;
E2 = 0.5;

V1_elderly = [10,10,10];
V2_elderly = [0,0,30];
SimPeriod = 10;
VT = zeros(SimPeriod,6);
VT(:,5) = ones(SimPeriod,1) * 10;
VT(:,6) = ones(SimPeriod,1) * 10;

S_elderly_prev = elderly_total - sum(E1*V1_elderly(1:end-2) + (E2 - E1)*V2_elderly(1:end-2));
V1_elderly_Sim = [V1_elderly(end-1);V1_elderly(end);VT(1:end-2,5)];
V2_elderly_Sim = [V2_elderly(end-1);V2_elderly(end);VT(1:end-2,6)];
S_elderly_path = S_elderly_prev - cumsum(E1 * V1_elderly_Sim+ (E2 - E1) * V2_elderly_Sim);
