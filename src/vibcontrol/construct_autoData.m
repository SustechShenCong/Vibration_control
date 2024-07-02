function autData = construct_autoData(DS, modes, R_0)
% CONSTRUCT_AUTODATA This function constructs autoData that is used for
% reduced dynamics

% extract coefficients and exponents
beta  = [];
kappa = [];
for k = 2:numel(R_0)
    R = R_0{k};
    betak = R.coeffs;
    if ~isempty(betak)
        kappak = R.ind;
        % assemble terms
        beta  = [beta betak];
        kappa = [kappa; kappak];
    end
end
autData       = struct();
autData.lamd  = DS.spectrum.Lambda(modes);
autData.beta  = beta;
autData.kappa = kappa;

end