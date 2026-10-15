function prior = ttPrior(model)
%TTPRIOR  Default Gaussian prior (mean and precision per parameter) for the MAP fit of a model; the empirical-Bayes start.
%
%   prior = ttPrior(model)
%
%   Penalty 0.5 * sum(w .* (theta - mu).^2). prior.wMin / wMax bound the precisions during the
%   empirical-Bayes updates of ttFitAnimal.
n = model.nParam; mu = zeros(n, 1); w = ones(n, 1);
for i = 1:n
    nm = model.names{i};
    if startsWith(nm, 'b_z1'),       mu(i) = log(0.3); w(i) = 0.5;
    elseif startsWith(nm, 'b_z'),    mu(i) = log(0.7); w(i) = 0.5;
    elseif startsWith(nm, 'b_logr'), mu(i) = 0;        w(i) = 0.5;
    elseif startsWith(nm, 'b_ell'),  mu(i) = -1;       w(i) = 1;
    elseif startsWith(nm, 'w_lsig'), mu(i) = 0;        w(i) = 2;
    elseif strcmp(nm, 'w_z_tauR'),   mu(i) = 0;        w(i) = 2;
    elseif startsWith(nm, 'w_'),     mu(i) = 0;        w(i) = 2;
    elseif startsWith(nm, 'alpha_v'), mu(i) = 0.01;    w(i) = 10;
    elseif any(strcmp(nm, {'alpha_up', 'alpha_down', 'alpha_q'})), mu(i) = 0.05; w(i) = 5;
    elseif strcmp(nm, 'alpha_w'),    mu(i) = 0.01;     w(i) = 10;
    elseif startsWith(nm, 'alpha'),  mu(i) = 0.01;     w(i) = 10;
    elseif strcmp(nm, 'm_z'),        mu(i) = 0.2;      w(i) = 1;
    elseif strcmp(nm, 'c_ssig'),     mu(i) = 0;        w(i) = 1;
    elseif strcmp(nm, 'lambda_b'),   mu(i) = 0.9;      w(i) = 1;
    elseif strcmp(nm, 'eps_b'),      mu(i) = 0.05;     w(i) = 5;
    elseif strcmp(nm, 'lambda_s'),   mu(i) = 0.5;      w(i) = 1;
    elseif any(strcmp(nm, {'c_sz', 'c_spi'})), mu(i) = 0; w(i) = 1;
    elseif strcmp(nm, 'gamma'),      mu(i) = 0.8;      w(i) = 10;
    elseif strcmp(nm, 'sigma'),      mu(i) = 0.3;      w(i) = 5;
    elseif strcmp(nm, 's1'),         mu(i) = 1.0;      w(i) = 2;
    end
end
prior = struct('mu', mu, 'w', w, 'wMin', 0.5 * ones(n, 1), 'wMax', 100 * ones(n, 1), 'names', {model.names}, 'source', 'default');
end
