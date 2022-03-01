function [avg, se] = ts_avg(x, sample_period, weight)

    x_sample            = x(sample_period);
    weighted_x_sample   = x_sample.*weight;
    
    avg     = sum(weighted_x_sample);
    se      = sqrt(var(x_sample)*sum((weight).^2));
     
end