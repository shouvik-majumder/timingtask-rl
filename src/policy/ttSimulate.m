function S = ttSimulate(theta, B, model, seed)
%TTSIMULATE  Generate a session from a policy model on a given delay schedule.
%
%   S = ttSimulate(theta, B, model, seed)
%
%   B.delay is the delay of each trial (s); the number of trials is numel(B.delay). A lick is
%   rewarded when it occurs at or after the delay. S.tau, S.R, S.delay, S.u (the channel values used
%   on each trial). The recursion is the one of ttForward, so ttForward on the simulated session
%   reproduces S.u exactly (ttGradCheck verifies this).
rs = RandStream('mt19937ar', 'Seed', seed);
L = model.nLags; nH = model.nH; chan = model.chan; nC = numel(chan); fam = model.family;
delay = B.delay(:); n = numel(delay);
signal = model.signal; wRule = model.wRule; width = model.width;
W = zeros(nC, 1 + nH); alpha = zeros(nC, 1); learn = false(nC, 1); hasH = false(nC, 1); rule = repmat({'pg'}, nC, 1); alpha2 = zeros(nC, 1); marg = zeros(nC, 1); wInt = 0; alphaW = 0; lamB = 0; epsB = 0;
for c = 1:nC
    W(c, 1) = theta(chan(c).iB);
    if ~isempty(chan(c).rule), rule{c} = chan(c).rule; end
    if ~isempty(chan(c).iW), W(c, 2:end) = theta(chan(c).iW); hasH(c) = true; end
    if ~isempty(chan(c).iA), alpha(c) = theta(chan(c).iA); learn(c) = true; end
    if ~isempty(chan(c).iA2), alpha2(c) = theta(chan(c).iA2); end
    if ~isempty(chan(c).iM), marg(c) = theta(chan(c).iM); end
    if ~isempty(chan(c).iX), wInt = theta(chan(c).iX); end
    if ~isempty(chan(c).iAw), alphaW = theta(chan(c).iAw); end
    if ~isempty(chan(c).iL), lamB = theta(chan(c).iL); epsB = theta(chan(c).iE); learn(c) = true; end
end
switch fam, case {'ln', 'lnexp'}, sig0 = theta(model.idx.sigma); sExtra = []; case 'lnmix2', sig0 = theta(model.idx.sigma); sExtra = theta(model.idx.s1); end
wSig = zeros(nH, 1); if strcmp(width, 'history'), for k = 1:L, wSig(k) = theta(model.idx.(sprintf('w_lsig_tau%d', k))); wSig(L + k) = theta(model.idx.(sprintf('w_lsig_R%d', k))); end, end
alpha_v = 0; gam = 0; alpha_q = 0;
switch signal
    case 'td', alpha_v = theta(model.idx.alpha_v); gam = theta(model.idx.gamma);
    case 'q',  alpha_q = theta(model.idx.alpha_q);
end
useState = model.state; lam_s = 0; c_sz = 0; c_spi = 0; c_ssig = 0; iEll = find(strcmp({chan.name}, 'ell'));
if useState, lam_s = theta(model.idx.lambda_s); c_sz = theta(model.idx.c_sz); if isfield(model.idx, 'c_spi'), c_spi = theta(model.idx.c_spi); end; if isfield(model.idx, 'c_ssig'), c_ssig = theta(model.idx.c_ssig); end, end
beta = model.beta; warm = model.warmup;
tau = zeros(n, 1); R = zeros(n, 1); U = zeros(n, nC); hLt = zeros(L, 1); hR = zeros(L, 1); rho = 0; V = 0; v = zeros(3, 1); st = 0; x0 = [1; zeros(nH, 1)]; zLo = log(0.03); zHi = log(8); kappa = 0.15;
isBayes = any(strcmp(rule, 'bayes')); if isBayes, bGrid = linspace(zLo, zHi, 120)'; logPost = zeros(120, 1); end
for t = 1:n
    if t <= warm, x = x0; else, x = [1; hLt; hR]; end
    u = W * x;
    if t > warm && wInt ~= 0, u(1) = u(1) + wInt * hLt(1) * hR(1); end
    if useState, u(1) = u(1) + c_sz * st; if ~isempty(iEll), u(iEll) = u(iEll) + c_spi * st; end, end
    if isBayes, pw = exp(logPost - max(logPost)); pw = pw / sum(pw); Bhat = pw' * bGrid; u(1) = u(1) - W(1, 1) + Bhat + marg(1); end
    sgt = sig0; if t > warm && strcmp(width, 'history'), sgt = sig0 * exp(wSig' * [hLt; hR]); end; if useState && c_ssig ~= 0, sgt = sgt * exp(c_ssig * st); end
    sg = [sgt, sExtra]; U(t, :) = u';
    switch fam
        case 'ln',     tau(t) = exp(u(1) + sgt * randn(rs));
        case 'lnexp',  pi1 = 1 / (1 + exp(-u(3))); if rand(rs) < pi1, tau(t) = -log(rand(rs)) / exp(u(2)); else, tau(t) = exp(u(1) + sgt * randn(rs)); end
        case 'lnmix2', pi1 = 1 / (1 + exp(-u(3))); if rand(rs) < pi1, tau(t) = exp(u(2) + sExtra * randn(rs)); else, tau(t) = exp(u(1) + sgt * randn(rs)); end
    end
    tau(t) = min(max(tau(t), 0.01), 30); R(t) = double(tau(t) >= delay(t)); lt = log(tau(t));
    [~, g] = ttPolicy(tau(t), u, sg, fam);
    rho = (1 - beta) * rho + beta * R(t);
    switch signal
        case 'rate',     d = R(t) - rho;
        case 'td',       d = R(t) - rho + gam * V - V; V = V + alpha_v * d;
        case 'q',        phi = [1; lt; lt^2]; d = R(t) - v' * phi; v = v + alpha_q * d * phi;
        case 'unsigned', d = abs(R(t) - rho);
        case 'raw',      d = R(t);
    end
    if t > warm
        for c = find(learn)'
            switch rule{c}
                case 'pg',        stp = alpha(c) * d * g(c); if hasH(c), W(c, :) = W(c, :) + stp * x'; else, W(c, 1) = W(c, 1) + stp; end
                case 'shift',     W(c, 1) = W(c, 1) - alpha(c) * d;
                case 'threshold', Bnd = W(c, 1) - marg(c); if R(t) == 1 && lt < Bnd, W(c, 1) = W(c, 1) + alpha2(c) * (lt - Bnd); elseif R(t) == 0 && lt > Bnd, W(c, 1) = W(c, 1) + alpha(c) * (lt - Bnd); end
                case 'qmax',      if v(3) < -1e-3, target = min(max(-v(2) / (2 * v(3)), zLo), zHi); W(c, 1) = W(c, 1) + alpha(c) * (target - W(c, 1)); end
                case 'stepq',     Bnd = W(c, 1) - marg(c); sgm = 1 / (1 + exp(-(lt - Bnd) / kappa)); W(c, 1) = W(c, 1) + alpha(c) * (R(t) - sgm) * (-sgm * (1 - sgm) / kappa);
                case 'bayes',     if R(t) == 1, lik = epsB + (1 - 2 * epsB) * (bGrid <= lt); else, lik = epsB + (1 - 2 * epsB) * (bGrid > lt); end; logPost = lamB * logPost + log(lik);
            end
        end
        if ~strcmp(wRule, 'none') && hasH(1), dw = d; if strcmp(wRule, 'unsigned'), dw = abs(d); end; W(1, 2:end) = W(1, 2:end) + alphaW * dw * g(1) * x(2:end)'; end
    end
    if useState, st = lam_s * st + R(t); end
    hLt = [lt; hLt(1:end-1)]; hR = [R(t); hR(1:end-1)];
end
S = struct('tau', tau, 'R', R, 'delay', delay, 'u', U);
end
