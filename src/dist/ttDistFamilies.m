function fam = ttDistFamilies(which)
%TTDISTFAMILIES  Candidate first-lick-time distributions, all with positive support.
%
%   fam = ttDistFamilies()          all families
%   fam = ttDistFamilies({'lognormal', 'lnmix2'})
%
%   Each element has
%     .name      short name used in tables
%     .label     description for figures
%     .npar      number of free parameters
%     .fit       @(x, T) -> theta (unconstrained vector), maximum likelihood on the sample x
%                (T = upper limit of the response window, taken from the TRAINING trials only;
%                used by the uniform component)
%     .logpdf    @(x, theta, T) -> log density
%     .cdf       @(x, theta, T) -> CDF
%     .summary   @(theta, T) -> struct of interpretable quantities (component medians in seconds,
%                widths in log seconds, weights, rates)
%
%   Families (K = number of components):
%     invgauss      inverse Gaussian                                       2
%     lognormal     log-normal                                             2
%     gamma         gamma                                                  2
%     weibull       Weibull                                                2
%     ig_unif       inverse Gaussian + uniform on (0, T]                   3
%     ln_unif       log-normal + uniform on (0, T]                         3
%     ln_exp        log-normal + exponential (mixture with a weight)       4
%     ln_race       log-normal process racing a constant-rate process:     3
%                   first lick = min(timed lick, constant-rate lick)
%                   p(x) = f_LN(x) e^{-hx} + h e^{-hx} (1 - F_LN(x))
%     lnmix2        two log-normals                                        5
%     lnmix3        three log-normals                                      8
%
%   All fits are by direct maximisation of the log-likelihood with fminunc from several starts
%   (mixtures start from an EM solution on log x). Parameters live on an unconstrained scale:
%   locations in log seconds, scales/rates as logs, weights as logits.

all = {'invgauss', 'lognormal', 'gamma', 'weibull', 'ig_unif', 'ln_unif', 'ln_exp', 'ln_race', 'lnmix2', 'lnmix3'};
if nargin < 1 || isempty(which), which = all; end
if ischar(which), which = {which}; end
fam = [];
for i = 1:numel(which)
    f = struct("name", which{i}, "label", "", "npar", 0, "fit", [], "logpdf", [], "cdf", [], "summary", []);
    switch which{i}
        case 'invgauss'
            f.label = 'inverse Gaussian'; f.npar = 2;
            f.logpdf = @(x, th, T) igLogpdf(x, exp(th(1)), exp(th(2)));
            f.cdf    = @(x, th, T) igCdf(x, exp(th(1)), exp(th(2)));
            f.fit    = @(x, T) igMle(x);                                                       % closed form
            f.summary = @(th, T) struct('median_s', igMedian(exp(th(1)), exp(th(2))), 'mu_s', exp(th(1)), 'lambda', exp(th(2)));
        case 'lognormal'
            f.label = 'log-normal'; f.npar = 2;
            f.logpdf = @(x, th, T) lnLogpdf(x, th(1), exp(th(2)));
            f.cdf    = @(x, th, T) normcdf((log(x) - th(1)) / exp(th(2)));
            f.fit    = @(x, T) [mean(log(x)), log(max(std(log(x)), 1e-3))];          % closed form
            f.summary = @(th, T) struct('median_s', exp(th(1)), 'width_log', exp(th(2)));
        case 'gamma'
            f.label = 'gamma'; f.npar = 2;
            f.logpdf = @(x, th, T) gamLogpdf(x, exp(th(1)), exp(th(2)));
            f.cdf    = @(x, th, T) gamcdf(x, exp(th(1)), exp(th(2)));
            f.fit    = @(x, T) log(gamfit(x));                                                 % toolbox MLE
            f.summary = @(th, T) struct('median_s', gaminv(0.5, exp(th(1)), exp(th(2))), 'shape', exp(th(1)), 'scale_s', exp(th(2)));
        case 'weibull'
            f.label = 'Weibull'; f.npar = 2;
            f.logpdf = @(x, th, T) wblLogpdf(x, exp(th(1)), exp(th(2)));
            f.cdf    = @(x, th, T) 1 - exp(-(x / exp(th(1))).^exp(th(2)));
            f.fit    = @(x, T) log(wblfit(x));                                                 % toolbox MLE
            f.summary = @(th, T) struct('median_s', exp(th(1)) * log(2)^(1 / exp(th(2))), 'scale_s', exp(th(1)), 'shape', exp(th(2)));
        case 'ig_unif'
            f.label = 'inverse Gaussian + uniform'; f.npar = 3;
            f.logpdf = @(x, th, T) mixLogpdf(igLogpdf(x, exp(th(1)), exp(th(2))), unifLogpdf(x, T), 1 - sig(th(3)), sig(th(3)));
            f.cdf    = @(x, th, T) (1 - sig(th(3))) .* igCdf(x, exp(th(1)), exp(th(2))) + sig(th(3)) .* min(x / T, 1);
            f.fit    = @(x, T) mleFit(@(th) -sum(mixLogpdf(igLogpdf(x, exp(th(1)), exp(th(2))), unifLogpdf(x, T), 1 - sig(th(3)), sig(th(3)))), ...
                                       [repmat(igMle(x), 3, 1), logit([0.05; 0.2; 0.4])]);
            f.summary = @(th, T) struct('median_s', igMedian(exp(th(1)), exp(th(2))), 'lambda', exp(th(2)), 'weight_unif', sig(th(3)), 'T_s', T);
        case 'ln_unif'
            f.label = 'log-normal + uniform'; f.npar = 3;
            f.logpdf = @(x, th, T) mixLogpdf(lnLogpdf(x, th(1), exp(th(2))), unifLogpdf(x, T), 1 - sig(th(3)), sig(th(3)));
            f.cdf    = @(x, th, T) (1 - sig(th(3))) .* normcdf((log(x) - th(1)) / exp(th(2))) + sig(th(3)) .* min(x / T, 1);
            f.fit    = @(x, T) mleFit(@(th) -sum(mixLogpdf(lnLogpdf(x, th(1), exp(th(2))), unifLogpdf(x, T), 1 - sig(th(3)), sig(th(3)))), ...
                                       [lnStarts(x), logit([0.05; 0.2; 0.4])]);
            f.summary = @(th, T) struct('median_s', exp(th(1)), 'width_log', exp(th(2)), 'weight_unif', sig(th(3)), 'T_s', T);
        case 'ln_exp'
            f.label = 'log-normal + exponential'; f.npar = 4;
            f.logpdf = @(x, th, T) mixLogpdf(lnLogpdf(x, th(1), exp(th(2))), expLogpdf(x, exp(th(3))), 1 - sig(th(4)), sig(th(4)));
            f.cdf    = @(x, th, T) (1 - sig(th(4))) .* normcdf((log(x) - th(1)) / exp(th(2))) + sig(th(4)) .* (1 - exp(-exp(th(3)) * x));
            f.fit    = @(x, T) mleFit(@(th) -sum(mixLogpdf(lnLogpdf(x, th(1), exp(th(2))), expLogpdf(x, exp(th(3))), 1 - sig(th(4)), sig(th(4)))), ...
                                       [lnStarts(x), log(1 ./ (mean(x) * [1; 0.5; 2])), logit([0.05; 0.2; 0.4])]);
            f.summary = @(th, T) struct('median_s', exp(th(1)), 'width_log', exp(th(2)), 'rate_per_s', exp(th(3)), 'weight_exp', sig(th(4)));
        case 'ln_race'
            f.label = 'log-normal racing a constant-rate process'; f.npar = 3;
            f.logpdf = @(x, th, T) raceLogpdf(x, th(1), exp(th(2)), exp(th(3)));
            f.cdf    = @(x, th, T) 1 - (1 - normcdf((log(x) - th(1)) / exp(th(2)))) .* exp(-exp(th(3)) * x);
            f.fit    = @(x, T) mleFit(@(th) -sum(raceLogpdf(x, th(1), exp(th(2)), exp(th(3)))), raceStarts(x));
            f.summary = @(th, T) struct('median_s', exp(th(1)), 'width_log', exp(th(2)), 'rate_per_s', exp(th(3)), ...
                                        'share_rate_process', raceShare(th(1), exp(th(2)), exp(th(3))));
        case 'lnmix2'
            f.label = 'two log-normals'; f.npar = 5;
            f.logpdf = @(x, th, T) lnmixLogpdf(x, th, 2);
            f.cdf    = @(x, th, T) lnmixCdf(x, th, 2);
            f.fit    = @(x, T) lnmixFit(x, 2);
            f.summary = @(th, T) lnmixSummary(th, 2);
        case 'lnmix3'
            f.label = 'three log-normals'; f.npar = 8;
            f.logpdf = @(x, th, T) lnmixLogpdf(x, th, 3);
            f.cdf    = @(x, th, T) lnmixCdf(x, th, 3);
            f.fit    = @(x, T) lnmixFit(x, 3);
            f.summary = @(th, T) lnmixSummary(th, 3);
        otherwise
            error('ttDistFamilies:unknown', 'unknown family %s', which{i});
    end
    if isempty(fam), fam = f; else, fam(end+1) = f; end
end
end

%% ---------------------------------------------------------------- densities
function lp = igLogpdf(x, mu, lam)
lp = 0.5 * (log(lam) - log(2 * pi) - 3 * log(x)) - lam .* (x - mu).^2 ./ (2 * mu^2 .* x);
end
function c = igCdf(x, mu, lam)
a = sqrt(lam ./ x); c = normcdf(a .* (x / mu - 1)) + exp(min(2 * lam / mu, 700)) .* normcdf(-a .* (x / mu + 1)); c = min(max(c, 0), 1);
end
function m = igMedian(mu, lam)
lo = 1e-4; hi = 100;
for i = 1:60, mid = sqrt(lo * hi); if igCdf(mid, mu, lam) < 0.5, lo = mid; else, hi = mid; end, end
m = sqrt(lo * hi);
end
function lp = lnLogpdf(x, m, s)
lp = -log(x) - log(s) - 0.5 * log(2 * pi) - (log(x) - m).^2 / (2 * s^2);
end
function lp = gamLogpdf(x, k, th)
lp = (k - 1) * log(x) - x / th - k * log(th) - gammaln(k);
end
function lp = wblLogpdf(x, lam, k)
lp = log(k) - log(lam) + (k - 1) * (log(x) - log(lam)) - (x / lam).^k;
end
function lp = unifLogpdf(x, T)
lp = -log(T) * ones(size(x)); lp(x > T) = -inf;
end
function lp = expLogpdf(x, rate)
lp = log(rate) - rate * x;
end
function lp = mixLogpdf(l1, l2, w1, w2)
mx = max(l1, l2); mx(~isfinite(mx)) = 0;
lp = mx + log(w1 * exp(l1 - mx) + w2 * exp(l2 - mx));
end
function lp = raceLogpdf(x, m, s, h)
z = (log(x) - m) / s; S = normcdf(-z);                     % survival of the timed process
fT = exp(lnLogpdf(x, m, s));
lp = -h * x + log(fT + h * S);
end
function sh = raceShare(m, s, h)
% probability that the constant-rate process fires first
g = @(x) h * exp(-h * x) .* normcdf(-(log(x) - m) / s);
sh = integral(g, 0, inf);
end
function lp = lnmixLogpdf(x, th, K)
[m, s, w] = lnmixUnpack(th, K);
L = zeros(numel(x), K);
for k = 1:K, L(:, k) = lnLogpdf(x(:), m(k), s(k)) + log(w(k)); end
mx = max(L, [], 2); lp = mx + log(sum(exp(L - mx), 2));
lp = reshape(lp, size(x));
end
function c = lnmixCdf(x, th, K)
[m, s, w] = lnmixUnpack(th, K); c = zeros(size(x));
for k = 1:K, c = c + w(k) * normcdf((log(x) - m(k)) / s(k)); end
end
function [m, s, w] = lnmixUnpack(th, K)
m = th(1:K); s = exp(th(K+1:2*K)); c = [0, th(2*K+1:3*K-1)]; e = exp(c - max(c)); w = e / sum(e);
end
function S = lnmixSummary(th, K)
[m, s, w] = lnmixUnpack(th, K); [m, o] = sort(m); s = s(o); w = w(o);
S = struct();
for k = 1:K, S.(sprintf('median%d_s', k)) = exp(m(k)); S.(sprintf('width%d_log', k)) = s(k); S.(sprintf('weight%d', k)) = w(k); end
end

%% ---------------------------------------------------------------- fitting
function th = mleFit(nll, starts)
% starts: one start per row. fminunc (quasi-Newton) from each; best kept. NaN-safe objective.
obj = @(th) safe(nll(th));
opts = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off', 'MaxIterations', 2000, 'MaxFunctionEvaluations', 4000, 'OptimalityTolerance', 1e-7, 'StepTolerance', 1e-9);
best = inf; th = starts(1, :);
for i = 1:size(starts, 1)
    try
        [t, v] = fminunc(obj, starts(i, :), opts);
        if v < best, best = v; th = t; end
    catch
    end
end
end
function v = safe(v)
if ~isfinite(v), v = 1e12; end
end
function S = igStarts(x)
mu = mean(x); lam = max(mu^3 / max(var(x), 1e-6), 1e-3);
S = [log(mu), log(lam); log(median(x)), log(lam * 3); log(mu), log(lam / 3)];
end
function S = lnStarts(x)
m = mean(log(x)); s = max(std(log(x)), 0.05);
S = [m, log(s); log(median(x(x >= median(x)))), log(0.3); m + 0.2, log(s * 0.7)];
end
function S = gamStarts(x)
k = max(mean(x)^2 / max(var(x), 1e-6), 0.1); th = mean(x) / k;
S = [log(k), log(th); log(k * 2), log(th / 2); log(k / 2), log(th * 2)];
end
function S = wblStarts(x)
S = [log(mean(x)), log(1.5); log(median(x)), log(2.5); log(mean(x)), log(0.8)];
end
function S = raceStarts(x)
m = log(median(x(x >= median(x)))); r = 1 / mean(x);
S = [m, log(0.3), log(r); m, log(0.5), log(r / 4); log(median(x)), log(0.3), log(r * 2)];
end
function th = lnmixFit(x, K)
% EM on log x from 3 replicates, then a quasi-Newton polish of the same log-likelihood
lx = log(x(:)); th = [];
best = inf;
for rep = 1:3
    try
        gm = fitgmdist(lx, K, 'Options', statset('MaxIter', 500), 'RegularizationValue', 1e-6, 'Replicates', 1, 'Start', 'plus');
        [m, o] = sort(gm.mu'); s = sqrt(squeeze(gm.Sigma))'; s = s(o); w = gm.ComponentProportion(o);
        t0 = [m, log(s), log(w(2:end) / w(1))];
        v = -sum(lnmixLogpdf(x, t0, K));
        if v < best, best = v; th = t0; end
    catch
    end
end
if isempty(th)
    q = quantile(lx, linspace(0.2, 0.8, K)); th = [q, log(0.3 * ones(1, K)), zeros(1, K - 1)];
end
th = mleFit(@(t) -sum(lnmixLogpdf(x, t, K)), th);
end
function y = sig(z), y = 1 ./ (1 + exp(-z)); end
function z = logit(p), z = log(p ./ (1 - p)); end
function th = igMle(x)
mu = mean(x); lam = numel(x) / sum(1 ./ x - 1 / mu); th = [log(mu), log(max(lam, 1e-6))];
end
