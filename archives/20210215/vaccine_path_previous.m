function [VM] = vaccine_path(VacTotal,medical,elderly,ordinary,elderly_start,ordinary_start,lag,a)
% Compute the path of vaccinations for each group

T = length(VacTotal);
VM = zeros(T,6);        % Record vaccinations (first and second times) for each group
for t = 1:T
    if VacTotal(t) == 0
        VM(t,:) = 0;
    else
        if VM(t-lag,1) > 0
            VM(t,2) = VM(t-lag,1);
            if sum(VM(1:t-1,1))+VacTotal(t)-VM(t,2) < medical
            VM(t,1) = VacTotal(t)-VM(t,2);
            else
            VM(t,1) = medical-sum(VM(1:t-1,1));
            VM(t,3) = VacTotal(t)-VM(t,1)-VM(t,2);
            end
        else
            if sum(VM(1:t-1,1))+VacTotal(t) < medical
            VM(t,1) = VacTotal(t);
            else
                VM(t,1) = medical - sum(VM(1:t-1,1));
                VM(t,3) = VacTotal(t)-VM(t,1);
            end
            
            
            
            
            
            
            
            
        else
            VM(t,1) = VacTotal(t);
        end
        
        
        
        
        
        
        
       
    elseif t >= elderly_start
        if VM(t-lag,1) > 0
            VM(t,2) = VM(t-lag,1);
            VM(t,1) = (1-a)*(VacTotal(t)-VM(t,2));
            VM(t,3) = 
        else
            VM(t,1) = VacTotal(t);
        end
            
            
            
            
        if sum(VM(1:t-1,1)) < medical
            if sum(VM(1:t-1,1))+VacTotal(t) < medical
                VM(t,1) = VacTotal(t);
            else
                VM(t,1) = medical - sum(VM(1:t-1,1));
                VM(t,2) = VacTotal(t) - VM(t,1);
            end
    elseif sum(VM(1:t-1,2)) < medical
        if sum(VM(1:t-1,2))+VacTotal(t) < medical
            VM(t,2) = VacTotal(t);
        else
            VM(t,2) = medical - sum(VM(1:t-1,2));
            VM(t,3) = VacTotal(t) - VM(t,2);
        end
    elseif sum(VM(1:t-1,3)) < elderly
        if sum(VM(1:t-1,3))+VacTotal(t) < elderly
            VM(t,3) = VacTotal(t);
        else
            VM(t,3) = medical - sum(VM(1:t-1,3));
            VM(t,4) = VacTotal(t) - VM(t,3);
        end
        
        
        
        
        
        
        
        