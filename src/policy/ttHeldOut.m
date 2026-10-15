function H = ttHeldOut(D, outDir, varargin)
%TTHELDOUT  Held-out evaluation of fitted models: forward split, single-trial prediction, simulated ceiling, shuffled-order null.
%
%   H = ttHeldOut(D, outDir)                       every <outDir>/fits/<animal>_<variant>.mat
%   H = ttHeldOut(D, outDir, 'variants', {...}, 'animals', {...}, 'nSim', 5, 'shuffle', {...variants...})
%
%   Per session and variant:
%     training fit    parameters fitted on the first trainFrac of trials under the animal's
%                     empirical-Bayes prior; the forward pass then runs through the whole session
%                     (state observed, parameters fixed) and the remaining trials are scored
%     scores          ll: mean held-out log-likelihood per trial
%                     R2log: 1 - SSE / SST of predMeanLog against log tau on the held-out trials
%                     maeLog: median |log tau - predMeanLog|; mae: median |tau - predMedian| (s)
%                     rLog: Pearson correlation of predMeanLog with log tau
%                     cover80: fraction of held-out PIT values in [0.1, 0.9]; ks: KS distance of PIT
%     benchmarks      train (median of the training trials) and run20 (median of the previous 20
%                     licks) with a log-normal of the training width around them
%     ceiling         the same scores on nSim sessions simulated from the training fit on the
%                     session's delay schedule
%     shuffle         for the variants listed in 'shuffle': trial order permuted, refitted and scored
%   Output: <outDir>/results/heldOut_sessions.csv (rows of earlier runs are kept when skipDone).
p = inputParser;
p.addParameter('variants', {}, @iscell); p.addParameter('animals', {}, @(x) iscell(x) || isstring(x)); p.addParameter('nSim', 5, @isnumeric);
p.addParameter('shuffle', {}, @iscell); p.addParameter('trainFrac', 0.7, @isnumeric); p.addParameter('nStarts', 3, @isnumeric);
p.addParameter('sessions', {}, @iscell); p.addParameter('skipDone', true, @islogical);
p.parse(varargin{:}); o = p.Results;
fd = fullfile(outDir, 'fits'); res = fullfile(outDir, 'results'); if ~exist(res, 'dir'), mkdir(res); end
T = D.trials(D.trials.cued, :);
files = dir(fullfile(fd, '*.mat')); rows = {};
prevFile = fullfile(res, 'heldOut_sessions.csv'); prev = [];
if o.skipDone && isfile(prevFile), prev = readtable(prevFile, 'TextType', 'string'); rows{end+1} = prev; end
for f = 1:numel(files)
    L = load(fullfile(fd, files(f).name), 'Res'); Res = L.Res; animal = char(Res.animal); variant = Res.model.shortName; model = Res.model; prior = Res.prior; P = Res.params;
    if ~isempty(o.animals) && ~ismember(animal, string(o.animals)), continue; end
    if ~isempty(o.variants) && ~ismember(variant, o.variants), continue; end
    if ~isempty(prev) && any(prev.animal == string(animal) & prev.variant == string(variant)), continue; end
    doShuffle = ismember(variant, o.shuffle);
    sess = P.session; if ~isempty(o.sessions), sess = sess(ismember(sess, string(o.sessions))); end
    nS = numel(sess); out = cell(nS, 1);
    fprintf('ttHeldOut %s %s: %d sessions\n', animal, variant, nS);
    parfor s = 1:nS
        q = sortrows(T(T.session == sess(s), :), 'trialNum');
        B = struct('tau', q.lickTime, 'R', double(q.rewarded), 'delay', q.delay);
        row = P(P.session == sess(s), {'animal', 'session', 'sessionIndex', 'stage', 'phase', 'nTrials', 'group'}); row.variant = string(variant);
        row = [row, oneSplit(B, model, prior, o, false)]; %#ok<AGROW>
        if doShuffle
            sh = oneSplit(B, model, prior, o, true); sh.Properties.VariableNames = strcat('shuf_', sh.Properties.VariableNames); row = [row, sh(:, {'shuf_ll', 'shuf_R2log', 'shuf_rLog', 'shuf_run20_ll', 'shuf_train_ll'})]; %#ok<AGROW>
        else
            row.shuf_ll = NaN; row.shuf_R2log = NaN; row.shuf_rLog = NaN; row.shuf_run20_ll = NaN; row.shuf_train_ll = NaN;
        end
        out{s} = row;
    end
    rows{end+1} = vertcat(out{:}); %#ok<AGROW>
end
H = vertcat(rows{:});
writetable(H, prevFile);
fprintf('ttHeldOut: %d rows\n', height(H));
end
function r = oneSplit(B, model, prior, o, shuffle)
n = numel(B.tau);
if shuffle
    rs = RandStream('mt19937ar', 'Seed', 7 + n); perm = randperm(rs, n); B = struct('tau', B.tau(perm), 'R', B.R(perm), 'delay', B.delay(perm));
end
nTr = round(o.trainFrac * n); mask = false(n, 1); mask(1:nTr) = true; te = ~mask;
theta = ttFitSession(B, model, prior, 'nStarts', o.nStarts, 'llMask', mask);
fw = ttForward(theta, B, model, struct('predict', true));
r = metrics(B.tau(te), fw.loglik(te), fw.predMeanLog(te), fw.predMedian(te), fw.pit(te));
lt = log(B.tau); sTr = std(lt(mask)); trMed = median(B.tau(mask));
run20 = nan(n, 1); for t = 2:n, run20(t) = median(B.tau(max(1, t - 20):t - 1)); end
bl = @(pred) -log(B.tau(te)) - log(sTr) - 0.9189385332046727 - (lt(te) - log(pred)).^2 / (2 * sTr^2);
b1 = metrics(B.tau(te), bl(trMed * ones(sum(te), 1)), log(trMed) * ones(sum(te), 1), trMed * ones(sum(te), 1), normcdf((lt(te) - log(trMed)) / sTr));
b2 = metrics(B.tau(te), bl(run20(te)), log(run20(te)), run20(te), normcdf((lt(te) - log(run20(te))) / sTr));
r.train_ll = b1.ll; r.train_R2log = b1.R2log; r.train_mae = b1.mae; r.train_maeLog = b1.maeLog;
r.run20_ll = b2.ll; r.run20_R2log = b2.R2log; r.run20_mae = b2.mae; r.run20_maeLog = b2.maeLog; r.run20_rLog = b2.rLog;
c = nan(o.nSim, 4);
for k = 1:o.nSim
    S = ttSimulate(theta, B, model, 100 + k); fs = ttForward(theta, S, model, struct('predict', true));
    m = metrics(S.tau(te), fs.loglik(te), fs.predMeanLog(te), fs.predMedian(te), fs.pit(te)); c(k, :) = [m.ll, m.R2log, m.mae, m.rLog];
end
r.ceil_ll = mean(c(:, 1)); r.ceil_R2log = mean(c(:, 2)); r.ceil_mae = mean(c(:, 3)); r.ceil_rLog = mean(c(:, 4));
r.nTest = sum(te); r.theta = string(jsonencode(theta'));
end
function m = metrics(tau, ll, predML, predMed, pit)
lt = log(tau); ok = isfinite(predML) & isfinite(predMed) & isfinite(lt);
if sum(ok) < 3
    m = table(mean(ll), NaN, NaN, NaN, NaN, mean(pit >= 0.1 & pit <= 0.9), ksPit(pit), 'VariableNames', {'ll', 'R2log', 'maeLog', 'mae', 'rLog', 'cover80', 'ks'}); return;
end
m = table(mean(ll), 1 - sum((lt(ok) - predML(ok)).^2) / sum((lt(ok) - mean(lt(ok))).^2), median(abs(lt(ok) - predML(ok))), median(abs(tau(ok) - predMed(ok))), ...
    corr(predML(ok), lt(ok)), mean(pit >= 0.1 & pit <= 0.9), ksPit(pit), 'VariableNames', {'ll', 'R2log', 'maeLog', 'mae', 'rLog', 'cover80', 'ks'});
end
function d = ksPit(u)
u = sort(u(isfinite(u))); n = numel(u); if n == 0, d = NaN; return; end
d = max(max((1:n)' / n - u), max(u - (0:n-1)' / n));
end
