function T = ttComponentResponsibility(D, outDir, varargin)
%TTCOMPONENTRESPONSIBILITY  Per-trial probability that a lick came from the early component of the session's log-normal + exponential fit.
%
%   T = ttComponentResponsibility(D, outDir)
%   T = ttComponentResponsibility(D, outDir, 'threshold', 0.5)
%
%   Reads the whole-session log-normal + exponential fits from <outDir>/results/dist_components.csv
%   (ttDistSummary) and returns D.trials (cued trials) with columns respEarly (responsibility of the
%   exponential component), early (respEarly above threshold) and prevEarly (component of the
%   previous cued trial of the session, NaN for the first).
p = inputParser; p.addParameter('threshold', 0.5, @isnumeric); p.parse(varargin{:}); o = p.Results;
CP = readtable(fullfile(outDir, 'results', 'dist_components.csv'), 'TextType', 'string'); E = CP(CP.family == "ln_exp" & CP.unitType == "session", :);
T = D.trials(D.trials.cued, :); T.respEarly = nan(height(T), 1);
for s = 1:height(E)
    m = T.session == E.session(s); if ~any(m), continue; end
    x = T.lickTime(m); w = E.weight_exp(s); r = E.rate_per_s(s); mu = log(E.median_s(s)); sg = E.width_log(s);
    fE = r * exp(-r * x); fT = exp(-log(x) - log(sg) - 0.5 * log(2 * pi) - (log(x) - mu).^2 / (2 * sg^2));
    T.respEarly(m) = w * fE ./ (w * fE + (1 - w) * fT);
end
T.early = T.respEarly > o.threshold; T.prevEarly = nan(height(T), 1);
T = sortrows(T, {'animal', 'session', 'trialNum'});
same = [false; T.session(2:end) == T.session(1:end-1)]; idx = find(same); T.prevEarly(idx) = double(T.early(idx - 1));
end
