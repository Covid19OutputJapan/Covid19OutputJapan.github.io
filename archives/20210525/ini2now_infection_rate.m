function var_initial_vector = ini2now_infection_rate(var_initial_vector,var_growth,var_ss_vector,SimPeriod,share_index)


for pindex = 1:5
    var_initial = var_initial_vector(pindex);
    var_ss = var_ss_vector(pindex);
%     logit_initial = log(var_initial/(1-var_initial)); % Logit of the variant share, most recently
    logit_initial = log(var_initial/(var_ss-var_initial)); % Logit of the variant share, most recently
    var_share = exp((1:SimPeriod)'*var_growth+logit_initial).*var_ss./(1+exp((1:SimPeriod)'*var_growth+logit_initial));
    var_initial_vector(pindex) = var_share(share_index); %var_share(2)
end

disp(var_initial_vector)

end