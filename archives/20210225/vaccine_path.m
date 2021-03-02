function [V,delta,VT,real_pace] = vaccine_path(pace,T,medical,elderly,ordinary,medical_start,elderly_start,lag,effectiveness,delta_average,elderly_total)
% Compute the path of vaccinations for each group

VM = zeros(10000,6);        % Record vaccinations (first and second times) for each group
pace_e = pace/2;%pace_elderly
pace_m = pace_e/4;%pace_medical

medical_vector = ones(ceil(medical/pace_m),1)*medical/ceil(medical/pace_m);
elderly_vector = ones(ceil(elderly/pace_e),1)*elderly/ceil(elderly/pace_e);
ordinary_vector = ones(ceil(ordinary/pace_e),1)*ordinary/ceil(ordinary/pace_e);
ordinary_start = elderly_start + length(elderly_vector);
real_pace = (elderly/ceil(elderly/pace_e))*2;

VM(medical_start:medical_start+length(medical_vector)-1,1) = medical_vector;
VM(medical_start+lag:medical_start+length(medical_vector)-1+lag,2) = medical_vector;
VM(elderly_start:elderly_start+length(elderly_vector)-1,3) = elderly_vector;
VM(elderly_start+lag:elderly_start+length(elderly_vector)-1+lag,4) = elderly_vector;
VM(ordinary_start:ordinary_start+length(ordinary_vector)-1,5) = ordinary_vector;
VM(ordinary_start+lag:ordinary_start+length(ordinary_vector)-1+lag,6) = ordinary_vector;

VT = VM(1:T,:);
V = effectiveness*(VT(:,2)+VT(:,4)+VT(:,6));
delta_ss = delta_average*(0.09/1.28);

delta = (1-(cumsum(VT(:,4))/elderly_total))*(delta_average-delta_ss)+delta_ss;
