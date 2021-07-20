function[V1_w,V2_w] = vaccine_daily_to_weekly(vaccine, ps, dateEN,iPC)

if iPC == 1
    %         dateV = datetime(vaccine.textdata(2:end,1),'InputFormat','yyyy/MM/dd');
    dateV = datetime(vaccine.textdata(2:end,1),'InputFormat','MM/dd/yyyy');
    vaccine_data = vaccine.data;
    totalV_d = vaccine_data(:,1);
    V1_d = vaccine_data(:,2);
    V2_d = vaccine_data(:,3);
else
    vaccine_data = vaccine.data;
    dateV = datetime(vaccine_data(:,1),'ConvertFrom','excel');
    totalV_d = vaccine_data(:,2);
    V1_d = vaccine_data(:,3);
    V2_d = vaccine_data(:,4);
end

V1_d(isnan(V1_d)) = 0;
V2_d(isnan(V2_d)) = 0;
Vdata = size(vaccine,1);
V1_w = zeros(length(dateEN),1);
V2_w = zeros(length(dateEN),1);
for i = 1:1:length(dateEN)
    index_date = (dateV >= dateEN(i)-3 & dateV <= dateEN(i)+3);
    V1_w(i) = sum(V1_d(index_date));
    V2_w(i) = sum(V2_d(index_date));
end
V1_w = ps*V1_w;
V2_w = ps*V2_w;