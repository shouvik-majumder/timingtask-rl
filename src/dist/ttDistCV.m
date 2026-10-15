function U = ttDistCV(D, outDir, varargin)
%TTDISTCV  Cross-validated comparison of first-lick-time distributions, per session and per delay block.
%
%   U = ttDistCV(D, outDir)
%   U = ttDistCV(D, outDir, 'families', {...}, 'sessions', {...}, 'parallel', true)
%
%   Units
%     session      all cued trials of a session (at least minTrials)
%     delayBlock   the trials at one delay level within a session (at least minBlock trials); a
%                  single-delay session is one block
%   Schemes
%     full         fit on all trials: in-sample log-likelihood per trial, BIC, fitted summaries
%     random5      5-fold random split, seeded per unit
%     forward      fit on the first trainFrac of trials in order, score the rest
%   Scores on held-out trials: mean log-likelihood per trial (testLL), Kolmogorov-Smirnov distance
%   between the held-out trials and the fitted CDF (testKS), coverage of the 10-90 % predictive
%   interval (cover80), median absolute error of the predictive median in seconds (medAE).
%   The uniform component's upper limit is the largest training lick time times 1.001.
%
%   Output: <outDir>/results/dist_units.csv, one row per unit x scheme x family.
p = inputParser;
p.addParameter('sessions', {}, @(x) iscell(x) || isstring(x)); p.addParameter('families', {}, @iscell);
p.addParameter('minTrials', 80, @isnumeric); p.addParameter('minBlock', 80, @isnumeric); p.addParameter('trainFrac', 0.7, @isnumeric);
p.addParameter('seed', 11, @isnumeric); p.addParameter('parallel', true, @islogical);
p.parse(varargin{:}); o = p.Results;
res = fullfile(outDir, 'results'); if ~exist(res, 'dir'), mkdir(res); end
T = D.trials(D.trials.cued, :); S = D.sessions;
if ~isempty(o.sessions), S = S(ismember(S.session, string(o.sessions)), :); end
fam = ttDistFamilies(o.families);
nS = height(S); rows = cell(nS, 1);
fprintf('ttDistCV: %d sessions, %d families\n', nS, numel(fam));
if o.parallel
    parfor s = 1:nS, rows{s} = oneSession(S(s, :), T(T.session == S.session(s), :), fam, o); end
else
    for s = 1:nS, rows{s} = oneSession(S(s, :), T(T.session == S.session(s), :), fam, o); end
end
U = vertcat(rows{:});
writetable(U, fullfile(res, 'dist_units.csv'));
fprintf('ttDistCV: %d rows written\n', height(U));
end
function out = oneSession(r, Ts, fam, o)
Ts = sortrows(Ts, 'trialNum'); tau = Ts.lickTime(:); dl = round(Ts.delay(:) * 20) / 20;
key = table(r.animal, r.sessionIndex, r.session, r.stage, r.phase, r.group, 'VariableNames', {'animal', 'sessionIndex', 'session', 'stage', 'phase', 'group'});
units = {};
if numel(tau) >= o.minTrials, units{end+1} = {"session", -1, tau}; end
for d = unique(dl)'
    m = dl == d; if sum(m) >= o.minBlock, units{end+1} = {"delayBlock", d, tau(m)}; end %#ok<AGROW>
end
rows = {};
for u = 1:numel(units)
    ut = units{u}{1}; dlev = units{u}{2}; x = units{u}{3}; n = numel(x);
    for sc = ["full", "random5", "forward"]
        switch sc
            case "full",    folds = zeros(n, 1); nf = 1;
            case "random5", rs = RandStream('threefry', 'Seed', o.seed + 1000 * u); folds = mod(randperm(rs, n)', 5) + 1; nf = 5;
            case "forward", folds = zeros(n, 1); folds(round(o.trainFrac * n) + 1:end) = 2; nf = 1;
        end
        for c = 1:numel(fam)
            ll = 0; cnt = 0; ks = []; cov = []; ae = []; summ = []; bic = NaN; nPar = fam(c).npar; ok = true;
            for f = 1:nf
                if sc == "full", te = false(n, 1); elseif sc == "random5", te = folds == f; else, te = folds == 2; end
                tr = ~te; xtr = x(tr); xte = x(te);
                if sum(tr) < 40 || (any(te) && sum(te) < 10), ok = false; break; end
                Tmax = max(xtr) * 1.001;
                try, th = fam(c).fit(xtr, Tmax); catch, ok = false; break; end
                if sc == "full"
                    lp = fam(c).logpdf(xtr, th, Tmax); ll = mean(lp); cnt = 1; bic = -2 * sum(lp) + nPar * log(n);
                    summ = fam(c).summary(th, Tmax);
                else
                    lp = fam(c).logpdf(xte, th, Tmax); lp = max(lp, -700);
                    ll = ll + sum(lp); cnt = cnt + numel(xte);
                    Fte = fam(c).cdf(xte, th, Tmax); ks(end+1) = ksDist(xte, Fte); %#ok<AGROW>
                    cov(end+1) = mean(Fte >= 0.1 & Fte <= 0.9); %#ok<AGROW>
                    md = predMedian(@(z) fam(c).cdf(z, th, Tmax)); ae = [ae; abs(xte - md)]; %#ok<AGROW>
                    if f == 1, summ = fam(c).summary(th, Tmax); end
                end
            end
            if ~ok || cnt == 0, continue; end
            if sc == "full", llPer = ll; else, llPer = ll / cnt; end
            rows{end+1} = [key, table(ut, dlev, n, sc, string(fam(c).name), nPar, llPer, bic, median(ks), median(cov), median(ae), string(jsonencode(summ)), ...
                'VariableNames', {'unitType', 'delayLevel', 'nUnit', 'scheme', 'family', 'nPar', 'testLL', 'bic', 'testKS', 'cover80', 'medAE', 'summary'})]; %#ok<AGROW>
        end
    end
end
if isempty(rows), out = []; else, out = vertcat(rows{:}); end
end
function d = ksDist(x, Fx)
[~, o] = sort(x); Fx = Fx(o); n = numel(x);
d = max(max((1:n)' / n - Fx(:)), max(Fx(:) - (0:n-1)' / n));
end
function md = predMedian(cdfFun)
lo = 1e-3; hi = 60;
for i = 1:50, mid = sqrt(lo * hi); if cdfFun(mid) < 0.5, lo = mid; else, hi = mid; end, end
md = sqrt(lo * hi);
end
