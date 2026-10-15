function C = ttCompare(outDir, varargin)
%TTCOMPARE  Summaries of fitted models: in-sample fit, held-out prediction and paired contrasts, by phase and group.
%
%   C = ttCompare(outDir)
%   C = ttCompare(outDir, 'reference', 'lnexp_Tfixed_Efixed', 'by', {'phase'})
%
%   Reads <outDir>/fits/*.mat and <outDir>/results/heldOut_sessions.csv and writes
%     compare_params.csv     fitted parameters per session for every variant (long format)
%     compare_insample.csv   variant x (by): n, median log-likelihood per trial, median BIC per
%                            trial, fraction of sessions where the variant has the best BIC
%     compare_heldout.csv    variant x (by) and 'all': median held-out scores and their ceilings
%                            and benchmarks
%     compare_contrasts.csv  every variant minus the reference, paired within session: median
%                            difference in held-out log-likelihood per trial and R2log, fraction of
%                            sessions where the variant is better, sign-rank p
%   C.insample, C.heldout, C.contrasts, C.params, C.heldoutSessions
p = inputParser; p.addParameter('reference', 'lnexp_Tfixed_Efixed', @ischar); p.addParameter('by', {'phase'}, @iscell); p.parse(varargin{:}); o = p.Results;
fd = fullfile(outDir, 'fits'); res = fullfile(outDir, 'results'); if ~exist(res, 'dir'), mkdir(res); end
files = dir(fullfile(fd, '*.mat')); P = {}; PL = {};
for f = 1:numel(files)
    L = load(fullfile(fd, files(f).name), 'Res'); q = L.Res.params; v = string(L.Res.model.shortName); q.variant = repmat(v, height(q), 1); q.nParam = repmat(L.Res.model.nParam, height(q), 1);
    P{end+1} = q(:, [{'animal', 'sessionIndex', 'session', 'stage', 'phase', 'group', 'nTrials', 'stageIndex', 'stageIndexFromEnd', 'variant', 'nParam', 'logLik', 'llPerTrial', 'BIC', 'AIC', 'exitflag'}]); %#ok<AGROW>
    pl = stack(q(:, [{'animal', 'session', 'sessionIndex', 'stage', 'phase', 'group', 'stageIndex', 'stageIndexFromEnd'}, L.Res.model.names]), L.Res.model.names, 'NewDataVariableName', 'value', 'IndexVariableName', 'param');
    pl.param = string(pl.param); pl.variant = repmat(v, height(pl), 1); PL{end+1} = pl; %#ok<AGROW>
end
if isempty(P), error('ttCompare:noFits', 'no fits in %s', fd); end
P = vertcat(P{:}); PL = vertcat(PL{:}); writetable(PL, fullfile(res, 'compare_params.csv'));
P.bicPerTrial = P.BIC ./ P.nTrials;
[g, K] = findgroups(P(:, {'session'})); bestBic = strings(max(g), 1);
for i = 1:max(g), q = P(g == i, :); [~, j] = min(q.BIC); bestBic(i) = q.variant(j); end
K.bestBic = bestBic; P = innerjoin(P, K, 'Keys', 'session');
[g2, I] = findgroups(P(:, [{'variant'}, o.by]));
I.n = accumarray(g2, 1); I.nParam = accumarray(g2, P.nParam, [], @(x) x(1)); I.median_llPerTrial = accumarray(g2, P.llPerTrial, [], @median); I.median_bicPerTrial = accumarray(g2, P.bicPerTrial, [], @median);
I.frac_bestBic = accumarray(g2, P.bestBic == P.variant, [], @mean);
writetable(I, fullfile(res, 'compare_insample.csv'));
C = struct('insample', I, 'params', PL);
hf = fullfile(res, 'heldOut_sessions.csv'); if ~isfile(hf), fprintf('ttCompare: no heldOut_sessions.csv\n'); return; end
H = readtable(hf, 'TextType', 'string'); H.all = repmat("all", height(H), 1);
mets = {'ll', 'R2log', 'maeLog', 'mae', 'rLog', 'cover80', 'ks', 'ceil_ll', 'ceil_R2log', 'ceil_mae', 'ceil_rLog', 'train_ll', 'train_R2log', 'train_mae', 'run20_ll', 'run20_R2log', 'run20_mae', 'run20_rLog', 'shuf_ll', 'shuf_R2log'};
rows = {};
for byCol = [o.by, {'all'}]
    [g3, Hk] = findgroups(H(:, [{'variant'}, byCol]));
    Hk.n = accumarray(g3, 1);
    for m = mets, Hk.(['median_', m{1}]) = accumarray(g3, H.(m{1}), [], @(x) median(x, 'omitnan')); end
    Hk.frac_ll_gt_run20 = accumarray(g3, H.ll > H.run20_ll, [], @mean); Hk.frac_R2log_gt_run20 = accumarray(g3, H.R2log > H.run20_R2log, [], @mean);
    Hk.Properties.VariableNames{2} = 'level'; Hk.by = repmat(string(byCol{1}), height(Hk), 1); rows{end+1} = Hk; %#ok<AGROW>
end
HO = vertcat(rows{:}); writetable(HO, fullfile(res, 'compare_heldout.csv')); C.heldout = HO; C.heldoutSessions = H;
rows = {}; ref = H(H.variant == string(o.reference), :);
for v = setdiff(unique(H.variant), string(o.reference))'
    a = H(H.variant == v, :); [tf, loc] = ismember(a.session, ref.session); a = a(tf, :); b = ref(loc(tf), :); if isempty(a), continue; end
    for byCol = [o.by, {'all'}]
        for lev = unique(a.(byCol{1}))'
            m = a.(byCol{1}) == lev; if sum(m) < 3, continue; end
            dl = a.ll(m) - b.ll(m); dr = a.R2log(m) - b.R2log(m); pv = NaN; if numel(dl) >= 5, pv = signrank(dl); end
            rows{end+1} = table(v, string(o.reference), string(byCol{1}), lev, sum(m), median(dl), mean(dl > 0), pv, median(dr, 'omitnan'), mean(dr > 0), ...
                'VariableNames', {'variant', 'reference', 'by', 'level', 'n', 'd_ll', 'frac_ll', 'p_ll', 'd_R2log', 'frac_R2log'}); %#ok<AGROW>
        end
    end
end
CT = vertcat(rows{:}); writetable(CT, fullfile(res, 'compare_contrasts.csv')); C.contrasts = CT;
end
