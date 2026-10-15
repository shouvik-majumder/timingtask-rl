function out = ttForward(theta, B, model, opt)
%TTFORWARD  Run a policy model through one session: log-likelihood and trial-by-trial latent variables.
%
%   out = ttForward(theta, B, model)
%   out = ttForward(theta, B, model, struct('llMask', mask, 'predict', true))
%
%   B.tau lick times (s), B.R rewards (0/1). Per trial t after warmup, x_t = [1; log tau_{t-1..t-L}; r_{t-1..t-L}]:
%     u_c(t) = W_c(t)' x_t (+ interaction, + state terms);  sigma_t = sigma exp(w_lsig' h_t + c_ssig s_t)
%     lp_t = log p(tau_t | u(t), sigma_t), scores g_c(t)                                 (ttPolicy)
%     delta_t = learning signal; the timed intercept is updated by model.rule, the history weights
%     by model.wRule, the early channels by the policy gradient; s_{t+1} = lambda_s s_t + r_t
%   out: loglik, totalLogLik (restricted to opt.llMask when given), u, score, step (move applied to
%   each channel's intercept after the trial), resp, delta, rho, V, s, W (final weights), sigmaT,
%   timedMedian, and with opt.predict: predMedian, predMeanLog, pit
if nargin < 4, opt = struct(); end
L = model.nLags; nH = model.nH; chan = model.chan; nC = numel(chan); fam = model.family;
tau = B.tau(:); R = double(B.R(:)); n = numel(tau); lt = log(tau);
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
switch fam
    case {'ln', 'lnexp'}, sig0 = theta(model.idx.sigma); sExtra = [];
    case 'lnmix2',        sig0 = theta(model.idx.sigma); sExtra = theta(model.idx.s1);
end
wSig = zeros(nH, 1); if strcmp(width, 'history'), for k = 1:L, wSig(k) = theta(model.idx.(sprintf('w_lsig_tau%d', k))); wSig(L + k) = theta(model.idx.(sprintf('w_lsig_R%d', k))); end, end
alpha_v = 0; gam = 0; alpha_q = 0;
switch signal
    case 'td', alpha_v = theta(model.idx.alpha_v); gam = theta(model.idx.gamma);
    case 'q',  alpha_q = theta(model.idx.alpha_q);
end
useState = model.state; lam_s = 0; c_sz = 0; c_spi = 0; c_ssig = 0; iEll = find(strcmp({chan.name}, 'ell'));
if useState, lam_s = theta(model.idx.lambda_s); c_sz = theta(model.idx.c_sz); if isfield(model.idx, 'c_spi'), c_spi = theta(model.idx.c_spi); end; if isfield(model.idx, 'c_ssig'), c_ssig = theta(model.idx.c_ssig); end, end
beta = model.beta; warm = model.warmup;
wantPred = isfield(opt, 'predict') && opt.predict;
ll = zeros(n, 1); U = zeros(n, nC); G = zeros(n, nC); S = zeros(n, nC); resp = zeros(n, 1); sigT = zeros(n, 1);
delta = zeros(n, 1); rhoS = zeros(n, 1); Vs = zeros(n, 1); sS = zeros(n, 1); predMed = nan(n, 1); predML = nan(n, 1); pit = nan(n, 1);
hLt = zeros(L, 1); hR = zeros(L, 1); rho = 0; V = 0; v = zeros(3, 1); st = 0;
x0 = [1; zeros(nH, 1)]; zLo = log(0.03); zHi = log(8); kappa = 0.15;
isBayes = any(strcmp(rule, 'bayes'));
if isBayes, bGrid = linspace(zLo, zHi, 120)'; logPost = zeros(120, 1); end
for t = 1:n
    if t <= warm, x = x0; else, x = [1; hLt; hR]; end
    u = W * x;
    if t > warm && wInt ~= 0, u(1) = u(1) + wInt * hLt(1) * hR(1); end
    if useState, u(1) = u(1) + c_sz * st; if ~isempty(iEll), u(iEll) = u(iEll) + c_spi * st; end, end
    if isBayes
        pw = exp(logPost - max(logPost)); pw = pw / sum(pw); Bhat = pw' * bGrid;
        u(1) = u(1) - W(1, 1) + Bhat + marg(1);          % the intercept is replaced by the posterior mean plus the margin
    end
    sgt = sig0; if t > warm && strcmp(width, 'history'), sgt = sig0 * exp(wSig' * [hLt; hR]); end; if useState && c_ssig ~= 0, sgt = sgt * exp(c_ssig * st); end
    sg = [sgt, sExtra]; sigT(t) = sgt;
    if wantPred
        [lp, g, rE, q] = ttPolicy(tau(t), u, sg, fam); predMed(t) = q.median; predML(t) = q.meanLog; pit(t) = q.cdf(tau(t));
    else
        [lp, g, rE] = ttPolicy(tau(t), u, sg, fam);
    end
    ll(t) = lp; U(t, :) = u'; G(t, :) = g'; resp(t) = rE; sS(t) = st;
    rho = (1 - beta) * rho + beta * R(t);
    switch signal
        case 'rate',     d = R(t) - rho;
        case 'td',       d = R(t) - rho + gam * V - V; Vs(t) = V; V = V + alpha_v * d;
        case 'q',        phi = [1; lt(t); lt(t)^2]; Qv = v' * phi; d = R(t) - Qv; Vs(t) = Qv; v = v + alpha_q * d * phi;
        case 'unsigned', d = abs(R(t) - rho);
        case 'raw',      d = R(t);
    end
    delta(t) = d; rhoS(t) = rho;
    if t > warm
        for c = find(learn)'
            switch rule{c}
                case 'pg'
                    stp = alpha(c) * d * g(c); S(t, c) = stp;
                    if hasH(c), W(c, :) = W(c, :) + stp * x'; else, W(c, 1) = W(c, 1) + stp; end
                case 'shift'
                    stp = -alpha(c) * d; S(t, c) = stp; W(c, 1) = W(c, 1) + stp;
                case 'threshold'
                    Bnd = W(c, 1) - marg(c); stp = 0;
                    if R(t) == 1 && lt(t) < Bnd, stp = alpha2(c) * (lt(t) - Bnd); elseif R(t) == 0 && lt(t) > Bnd, stp = alpha(c) * (lt(t) - Bnd); end
                    S(t, c) = stp; W(c, 1) = W(c, 1) + stp;
                case 'qmax'
                    stp = 0;
                    if v(3) < -1e-3, target = min(max(-v(2) / (2 * v(3)), zLo), zHi); stp = alpha(c) * (target - W(c, 1)); end
                    S(t, c) = stp; W(c, 1) = W(c, 1) + stp;
                case 'stepq'
                    Bnd = W(c, 1) - marg(c); sgm = 1 / (1 + exp(-(lt(t) - Bnd) / kappa)); dq = R(t) - sgm;
                    stp = alpha(c) * dq * (-sgm * (1 - sgm) / kappa); S(t, c) = stp; W(c, 1) = W(c, 1) + stp;
                case 'bayes'
                    if R(t) == 1, lik = epsB + (1 - 2 * epsB) * (bGrid <= lt(t)); else, lik = epsB + (1 - 2 * epsB) * (bGrid > lt(t)); end
                    logPost = lamB * logPost + log(lik); S(t, c) = NaN;
            end
        end
        if ~strcmp(wRule, 'none') && hasH(1)
            dw = d; if strcmp(wRule, 'unsigned'), dw = abs(d); end
            W(1, 2:end) = W(1, 2:end) + alphaW * dw * g(1) * x(2:end)';
        end
    end
    if useState, st = lam_s * st + R(t); end
    hLt = [lt(t); hLt(1:end-1)]; hR = [R(t); hR(1:end-1)];
end
if isBayes, S(:, 1) = [NaN; diff(U(:, 1))]; end
out = struct('loglik', ll, 'totalLogLik', sum(ll), 'u', U, 'score', G, 'step', S, 'resp', resp, 'delta', delta, 'rho', rhoS, 'V', Vs, 's', sS, 'W', W, 'sigmaT', sigT, ...
    'predMedian', predMed, 'predMeanLog', predML, 'pit', pit, 'timedMedian', exp(U(:, 1)));
if isfield(opt, 'llMask') && ~isempty(opt.llMask), out.totalLogLik = sum(ll(opt.llMask)); end
end
