function N = ttDistNoCue(D, outDir, varargin)
%TTDISTNOCUE  Lick timing on no-cue trials next to the fitted early component of the cued distribution.
%
%   N = ttDistNoCue(D, outDir)
%   N = ttDistNoCue(D, outDir, 'window', 5, 'minLick', 8)
%
%   No-cue trials enter no fit. Per session: the number of no-cue trials and of those with a lick
%   within the response window, the constant-hazard rate of those licks (licks per second of
%   exposure, trials without a lick censored at the window), their median time, and the
%   log-likelihood per trial of a censored log-normal minus that of a censored exponential. These
%   are joined with the exponential rate and weight of the log-normal + exponential full-session
%   fit from dist_components.csv.
%   Output: <outDir>/results/dist_noCue.csv
p = inputParser; p.addParameter('window', 5.0, @isnumeric); p.addParameter('minLick', 8, @isnumeric); p.parse(varargin{:}); o = p.Results;
res = fullfile(outDir, 'results'); CP = readtable(fullfile(res, 'dist_components.csv'), 'TextType', 'string');
S = D.sessions; T = D.trials; rows = {};
for s = 1:height(S)
    q = T(T.session == S.session(s) & T.noCue, :); if isempty(q), continue; end
    t = q.lickTime; lick = ~isnan(t) & t > 0 & t <= o.window; nL = sum(lick); nC = sum(~lick);
    r = S(s, {'animal', 'sessionIndex', 'session', 'stage', 'phase', 'group'}); r.nNoCue = height(q); r.nNoCueLick = nL; r.fracLick = nL / height(q);
    r.hazard_noCue = nL / (sum(t(lick)) + o.window * nC);
    r.median_noCue = NaN; r.dLL_LN_minus_exp = NaN; r.ln_median_noCue = NaN; r.ln_width_noCue = NaN;
    if nL >= o.minLick
        r.median_noCue = median(t(lick)); x = t(lick); W = o.window; h = r.hazard_noCue; llE = sum(log(h) - h * x) - h * W * nC;
        nllLN = @(th) -(sum(-log(x) - th(2) - 0.5 * log(2 * pi) - (log(x) - th(1)).^2 / (2 * exp(2 * th(2)))) + nC * log(max(1 - normcdf((log(W) - th(1)) / exp(th(2))), 1e-300)));
        th = fminsearch(nllLN, [mean(log(x)), log(max(std(log(x)), 0.1))], optimset('Display', 'off', 'MaxFunEvals', 2000));
        r.dLL_LN_minus_exp = (-nllLN(th) - llE) / height(q); r.ln_median_noCue = exp(th(1)); r.ln_width_noCue = exp(th(2));
    end
    rows{end+1} = r; %#ok<AGROW>
end
N = vertcat(rows{:});
pick = @(fam, col) arrayfun(@(b) firstOr(CP, b, fam, col), N.session);
N.exp_rate_cued = pick("ln_exp", "rate_per_s"); N.exp_weight_cued = pick("ln_exp", "weight_exp"); N.timed_median_cued = pick("ln_exp", "median_s");
writetable(N, fullfile(res, 'dist_noCue.csv'));
end
function v = firstOr(CP, b, fam, col)
v = NaN; if ~ismember(col, CP.Properties.VariableNames), return; end
x = CP.(col)(CP.session == b & CP.family == fam & CP.unitType == "session"); if ~isempty(x), v = x(1); end
end
