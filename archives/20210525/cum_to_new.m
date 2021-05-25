function new = cum_to_new(cum)

new = zeros(length(cum),1);
cum(isnan(cum)) = 0;
for vi = 2:1:length(cum)
    if  cum(vi)>0
        new(vi) = cum(vi) - cum(vi-1);
    end
end

