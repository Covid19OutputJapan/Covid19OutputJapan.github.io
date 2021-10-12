function [V,delta,VT] = vaccine_distribution(vacpath,elderly,ordinary,elderly_total,delta_average,E1,E2,D1,D2,lag)
% Compute the path of vaccinations for each group

T = length(vacpath);
VT = zeros(T,4);
for t = 1:T
    CV = cumsum(VT);
    if t <= lag
        if CV(t,1) + vacpath(t) <= elderly
            VT(t,1) = vacpath(t);
        elseif CV(t,1) + vacpath(t) > elderly && CV(t,1) < elderly
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1);
        elseif CV(t,1) >= elderly
            if CV(t,3) + vacpath(t) <= ordinary
                VT(t,3) = vacpath(t);
            elseif CV(t,3) + vacpath(t) > ordinary && CV(t,3) < ordinary
                VT(t,3) = ordinary - CV(t,3);
            end
        end
    elseif CV(t,1) + vacpath(t) <= elderly
        if VT(t-lag,1) == 0
            VT(t,1) = vacpath(t);
        elseif VT(t-lag,1) > 0
            VT(t,2) = VT(t-lag,1);
            VT(t,1) = vacpath(t) - VT(t,2);
        end
    elseif CV(t,1) + vacpath(t) > elderly && CV(t,1) + vacpath(t) - VT(t-lag,1) <= elderly
        if VT(t-lag,1) == 0
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1);
        elseif VT(t-lag,1) > 0
            VT(t,2) = VT(t-lag,1);
            VT(t,1) = vacpath(t) - VT(t,2);
        end
    elseif CV(t,1) + vacpath(t) - VT(t-lag,1) >= elderly && CV(t,1) < elderly
        if VT(t-lag,1) == 0
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1);
        elseif VT(t-lag,1) > 0
            VT(t,2) = VT(t-lag,1);
            VT(t,1) = elderly - CV(t,1);
            VT(t,3) = vacpath(t) - VT(t,1) - VT(t,2);
        end
    elseif CV(t,1) >= elderly
        if VT(t-lag,1) > 0
            if VT(t-lag,3) == 0
                VT(t,2) = VT(t-lag,1);
                VT(t,3) = vacpath(t) - VT(t,2);
            else
                VT(t,2) = VT(t-lag,1);
                VT(t,4) = VT(t-lag,3);
                VT(t,3) = vacpath(t) - VT(t,2) - VT(t,4);
            end
        elseif CV(t,3) < ordinary
            if VT(t-lag,3) > 0
                VT(t,4) = VT(t-lag,3);
                VT(t,3) = vacpath(t) - VT(t,4);
            else
                VT(t,3) = vacpath(t);
            end
        elseif CV(t,3) >= ordinary
            if VT(t-lag,3) > 0
                VT(t,4) = VT(t-lag,3);
            end
        end
    end
end

V1 = E1*(VT(:,1)+VT(:,3))+(E2-E1)*(VT(:,2)+VT(:,4));
V2 = D1*VT(:,1)+(D2-D1)*VT(:,2);
V2 = [0;0;V2(1:end-2)];
V = [0;0;V1(1:end-2)];
delta_ss = delta_average*(0.09/1.28);

delta = (1-(cumsum(V2)/elderly_total))*(delta_average-delta_ss)+delta_ss;