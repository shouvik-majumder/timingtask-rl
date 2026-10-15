function ok = ttGradCheck(verbose)
%TTGRADCHECK  Consistency checks of the policy code: scores against finite differences, normalisation, CDF, simulate-then-forward.
%
%   ok = ttGradCheck()
%
%   1. ttPolicy: analytic scores against central finite differences for every family; the density
%      integrates to one; the CDF handle matches the integrated density; E[log tau] and the median
%      match numerical integration.
%   2. ttForward: the gradient of the total log-likelihood with respect to each channel intercept of
%      a static model matches the summed scores.
%   3. ttSimulate then ttForward: for every learning signal, rule and the state, the channel values
%      the forward pass reconstructs on a simulated session equal the ones the simulation used.
if nargin < 1, verbose = true; end
ok = true; rng(3);
fams = {'ln', 'lnexp', 'lnmix2'};
for f = fams
    for rep = 1:20
        switch f{1}
            case 'ln',     u = log(0.3 + 2 * rand); s = 0.1 + 0.5 * rand;
            case 'lnexp',  u = [log(0.3 + 2 * rand); log(0.2 + 3 * rand); -2 + 4 * rand]; s = 0.1 + 0.5 * rand;
            case 'lnmix2', u = [log(0.5 + 2 * rand); log(0.05 + 0.6 * rand); -2 + 4 * rand]; s = [0.1 + 0.5 * rand, 0.2 + 1.2 * rand];
        end
        tau = 0.02 + 3 * rand;
        [~, g, ~, q] = ttPolicy(tau, u, s, f{1});
        gfd = zeros(size(u)); h = 1e-6;
        for k = 1:numel(u), up = u; up(k) = up(k) + h; um = u; um(k) = um(k) - h; gfd(k) = (ttPolicy(tau, up, s, f{1}) - ttPolicy(tau, um, s, f{1})) / (2 * h); end
        err = max(abs(g - gfd) ./ max(1, abs(gfd)));
        if err > 1e-5, fprintf('score mismatch %s rep %d: %.2e\n', f{1}, rep, err); ok = false; end
        if rep <= 3
            Z = integral(@(x) arrayfun(@(xx) exp(ttPolicy(xx, u, s, f{1})), x), 1e-8, 200);
            cdfNum = integral(@(x) arrayfun(@(xx) exp(ttPolicy(xx, u, s, f{1})), x), 1e-8, tau);
            if abs(Z - 1) > 1e-4 || abs(cdfNum - q.cdf(tau)) > 1e-4, fprintf('normalisation or CDF mismatch %s: Z=%.5f cdf %.5f vs %.5f\n', f{1}, Z, cdfNum, q.cdf(tau)); ok = false; end
            ml = integral(@(x) arrayfun(@(xx) log(xx) * exp(ttPolicy(xx, u, s, f{1})), x), 1e-8, 200);
            if abs(ml - q.meanLog) > 1e-3, fprintf('E[log tau] mismatch %s: %.4f vs %.4f\n', f{1}, ml, q.meanLog); ok = false; end
            if abs(q.cdf(q.median) - 0.5) > 1e-3, fprintf('median mismatch %s\n', f{1}); ok = false; end
        end
    end
end
model = ttModel('family', 'lnexp', 'timed', 'fixed', 'early', 'fixed');
B = struct('tau', 0.1 + 2 * rand(200, 1), 'R', double(rand(200, 1) > 0.5));
theta = [log(0.8); log(1); -1; 0.3];
out = ttForward(theta, B, model);
for c = 1:3
    h = 1e-6; tp = theta; tp(c) = tp(c) + h; tm = theta; tm(c) = tm(c) - h;
    gfd = (ttForward(tp, B, model).totalLogLik - ttForward(tm, B, model).totalLogLik) / (2 * h);
    if abs(gfd - sum(out.score(:, c))) > 1e-4, fprintf('forward score mismatch channel %d: %.5f vs %.5f\n', c, gfd, sum(out.score(:, c))); ok = false; end
end
Bs = struct('tau', 0.1 + 2 * rand(300, 1), 'R', double(rand(300, 1) > 0.5), 'delay', [0.5 * ones(150, 1); 1.0 * ones(150, 1)]);
cfgs = {{'signal', 'rate', 'rule', 'pg'}, {'signal', 'td', 'rule', 'pg'}, {'signal', 'q', 'rule', 'pg'}, {'signal', 'unsigned', 'rule', 'pg'}, {'signal', 'raw', 'rule', 'shift'}, {'rule', 'threshold'}, {'rule', 'threshold', 'margin', true, 'tieRates', false}, ...
    {'signal', 'q', 'rule', 'qmax'}, {'signal', 'rate', 'rule', 'pg', 'state', true}, {'rule', 'stepq'}, {'rule', 'bayes'}, {'rule', 'threshold', 'wRule', 'pg'}, {'rule', 'threshold', 'wRule', 'unsigned', 'interaction', true}, {'rule', 'threshold', 'width', 'history'}, {'rule', 'threshold', 'width', 'state', 'state', true}, {'rule', 'pg', 'nLags', 2}};
for i = 1:numel(cfgs)
    m = ttModel('family', 'lnexp', 'timed', 'historyLearn', 'early', 'historyLearn', cfgs{i}{:});
    th = (m.lb + m.ub)' / 2; th(m.idx.b_z) = log(0.8); th(m.idx.b_logr) = 0; th(m.idx.b_ell) = -1.5; th(m.idx.sigma) = 0.3;
    for nm = m.names, if startsWith(nm{1}, 'alpha'), th(m.idx.(nm{1})) = 0.05; end; if startsWith(nm{1}, 'w_'), th(m.idx.(nm{1})) = 0.1; end, end
    if isfield(m.idx, 'm_z'), th(m.idx.m_z) = 0.2; end; if isfield(m.idx, 'lambda_s'), th(m.idx.lambda_s) = 0.7; th(m.idx.c_sz) = -0.2; th(m.idx.c_spi) = 0.5; end
    if isfield(m.idx, 'lambda_b'), th(m.idx.lambda_b) = 0.9; th(m.idx.eps_b) = 0.05; end; if isfield(m.idx, 'c_ssig'), th(m.idx.c_ssig) = 0.1; end; if isfield(m.idx, 'w_z_tauR'), th(m.idx.w_z_tauR) = 0.1; end
    Ssim = ttSimulate(th, Bs, m, 5); o2 = ttForward(th, Ssim, m);
    err = max(abs(o2.u(:) - Ssim.u(:)));
    if err > 1e-10, fprintf('simulate/forward mismatch %s: %.2e\n', m.shortName, err); ok = false; elseif verbose, fprintf('  consistent: %s (%d parameters)\n', m.shortName, m.nParam); end
end
if ok && verbose, fprintf('ttGradCheck: all checks passed\n'); end
end
