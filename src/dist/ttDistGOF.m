function G = ttDistGOF(D, outDir, varargin)
%TTDISTGOF  Split-half goodness of fit of a lick-time distribution family, per session, against the agreement of the data with itself.
%
%   G = ttDistGOF(D, outDir)
%   G = ttDistGOF(D, outDir, 'family', 'ln_exp', 'nSplit', 50, 'minTrials', 80)
%
%   For every session with at least minTrials cued trials and every random half split (nSplit splits):
%     ks_model   Kolmogorov-Smirnov distance between the CDF fitted to half A and the empirical CDF of half B
%     ks_data    Kolmogorov-Smirnov distance between the empirical CDFs of half A and half B (the
%                reference: how far two samples of the same session are from each other)
%     r2_quant   R2 of the fitted-model quantiles of half A against the empirical quantiles of half B
%                (quantiles 0.05 to 0.95 in steps of 0.05, in log seconds)
%   Per session: medians over splits, the fraction of splits in which ks_model exceeds ks_data, and the
%   95th percentile of ks_data. A family that fits as well as the data agree with themselves has
%   ks_model close to ks_data and that fraction near 0.5.
%   Output: <outDir>/results/dist_gof_<family>.csv
p = inputParser; p.addParameter('family', 'ln_exp', @ischar); p.addParameter('nSplit', 50, @isnumeric); p.addParameter('minTrials', 80, @isnumeric); p.addParameter('seed', 5, @isnumeric); p.parse(varargin{:}); o = p.Results;
res = fullfile(outDir, 'results'); if ~exist(res, 'dir'), mkdir(res); end
fam = ttDistFamilies({o.family}); T = D.trials(D.trials.cued, :); S = D.sessions(D.sessions.nCued >= o.minTrials, :); rows = cell(height(S), 1); qs = 0.05:0.05:0.95;
parfor s = 1:height(S)
    x = T.lickTime(T.session == S.session(s)); n = numel(x); rs = RandStream('mt19937ar', 'Seed', o.seed + s);
    km = nan(o.nSplit, 1); kd = km; r2 = km;
    for k = 1:o.nSplit
        pm = randperm(rs, n); a = x(pm(1:floor(n / 2))); b = x(pm(floor(n / 2) + 1:end)); Tm = max(a) * 1.001;
        try, th = fam.fit(a, Tm); catch, continue; end
        km(k) = ksDist(b, fam.cdf(b, th, Tm)); kd(k) = ksTwo(a, b);
        qb = quantile(log(b), qs); qm = log(arrayfun(@(q) invCdf(@(z) fam.cdf(z, th, Tm), q), qs)); r2(k) = 1 - sum((qb - qm).^2) / sum((qb - mean(qb)).^2);
    end
    rows{s} = [S(s, {'animal', 'session', 'sessionIndex', 'stage', 'phase', 'group', 'nCued'}), table(string(o.family), median(km, 'omitnan'), median(kd, 'omitnan'), prctile(kd, 95), mean(km > kd), median(r2, 'omitnan'), ...
        'VariableNames', {'family', 'ks_model', 'ks_data', 'ks_data_p95', 'frac_model_worse', 'r2_quant'})];
end
G = vertcat(rows{:}); writetable(G, fullfile(res, ['dist_gof_', o.family, '.csv']));
end
function d = ksDist(x, Fx)
[~, o] = sort(x); Fx = Fx(o); n = numel(x); d = max(max((1:n)' / n - Fx(:)), max(Fx(:) - (0:n-1)' / n));
end
function d = ksTwo(a, b)
z = sort([a; b]); Fa = arrayfun(@(v) mean(a <= v), z); Fb = arrayfun(@(v) mean(b <= v), z); d = max(abs(Fa - Fb));
end
function m = invCdf(cdf, q)
lo = 1e-3; hi = 60; for i = 1:40, mid = sqrt(lo * hi); if cdf(mid) < q, lo = mid; else, hi = mid; end, end; m = sqrt(lo * hi);
end
