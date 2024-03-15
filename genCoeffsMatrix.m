function coeffs_matrix = genCoeffsMatrix(fitParams,n_synth_prof)

[n_coeff,~] = size(fitParams);
coeffs_matrix = zeros(n_synth_prof,n_coeff);
for i_coeff=1:1:n_coeff
    coeff_mu = fitParams(i_coeff,1);
    coeff_sd = fitParams(i_coeff,2);
    coeffs_matrix(:,i_coeff) = coeff_mu + coeff_sd * randn(n_synth_prof,1); % Generation of coefficients with extracted variability
end