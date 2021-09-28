function [betaT] = beta_AR1(betaT_temp_ini, beta_rho, betaT, start_beta)

betaT_temp = betaT_temp_ini;
betaT(start_beta,1) = (1+betaT_temp)*betaT(start_beta,1);
for i = start_beta+1:length(betaT)
    betaT_temp = betaT_temp * beta_rho;
    betaT(i) = betaT(i) * (1+betaT_temp);
end