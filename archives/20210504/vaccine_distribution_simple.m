function [V,delta,VT] = vaccine_distribution_simple(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2)
% Compute the path of vaccinations for each group

T = length(vacpath);
VT = zeros(T,4);
for t = 1:T
    CV = cumsum(VT);
    if CV(t,1) + vacpath(t) <= elderly
        VT(t,1) = vacpath(t);
    elseif CV(t,1) + vacpath(t) > elderly && CV(t,1) < elderly
        VT(t,1) = elderly - CV(t,1);
        VT(t,2) = vacpath(t) - VT(t,1);
    elseif CV(t,2) + vacpath(t) <= elderly
        VT(t,2) = vacpath(t);
    elseif CV(t,2) + vacpath(t) > elderly && CV(t,2) < elderly
        VT(t,2) = elderly - CV(t,2);
        VT(t,3) = vacpath(t) - VT(t,2);
    elseif CV(t,3) + vacpath(t) <= ordinary
        VT(t,3) = vacpath(t);
    elseif CV(t,3) + vacpath(t) > ordinary && CV(t,3) < ordinary
        VT(t,3) = ordinary - CV(t,3);
        VT(t,4) = vacpath(t) - VT(t,3);
    elseif CV(t,4) + vacpath(t) <= ordinary
        VT(t,4) = vacpath(t);
    elseif CV(t,4) + vacpath(t) > ordinary && CV(t,4) < ordinary
        VT(t,4) = ordinary - CV(t,4);
    end
end

V1 = E1*(VT(:,1)+VT(:,3))+(E2-E1)*(VT(:,2)+VT(:,4));
V2 = D1*VT(:,1)+(D2-D1)*VT(:,2);
V2 = [0;0;V2(1:end-2)];
V = [0;0;V1(1:end-2)];
V_ord = D1*(VT(:,3))+D2*(VT(:,4));
V_ord = [0;0;V_ord(1:end-2)];

delta_ss = delta_average*(0.09/1.28);

delta = (1-(cumsum(V2)/elderly_total))*(delta_average-delta_ss)+(1-(cumsum(V_ord)/(population-elderly_total)))*delta_ss;
