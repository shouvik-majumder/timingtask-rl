function S = ttDistSummary(outDir, varargin)
%TTDISTSUMMARY  Tables from ttDistCV: which distribution describes first-lick times, by phase and group.
%
%   S = ttDistSummary(outDir)
%   S = ttDistSummary(outDir, 'ref', 'lognormal', 'by', {'phase'})
%
%   Reads <outDir>/results/dist_units.csv and writes
%     dist_byPhase.csv     unitType x scheme x (by) x family: n, median held-out log-likelihood
%                          per trial, median paired difference to the reference family, fraction
%                          of units where the family is best, median KS, coverage, median AE
%     dist_winners.csv     per unit x scheme: best family, best among families with at most 3 parameters
%     dist_components.csv  fitted quantities of every family from the full fits, one row per unit
%   S.units, S.byPhase, S.winners, S.components.
p = inputParser; p.addParameter('ref', 'lognormal', @ischar); p.addParameter('by', {'phase'}, @iscell); p.parse(varargin{:}); o = p.Results;
res = fullfile(outDir, 'results');
U = readtable(fullfile(res, 'dist_units.csv'), 'TextType', 'string');
U.unitKey = U.session + "|" + U.unitType + "|" + string(U.delayLevel);
fams = unique(U.family, 'stable'); nParOf = containers.Map(cellstr(fams), num2cell(arrayfun(@(f) U.nPar(find(U.family == f, 1)), fams)));
U.dLL = nan(height(U), 1);
[g, K] = findgroups(U(:, {'unitKey', 'scheme'}));
for i = 1:max(g)
    idx = find(g == i); r = U.testLL(idx(U.family(idx) == o.ref)); if isempty(r), continue; end
    U.dLL(idx) = U.testLL(idx) - r;
end
win = strings(max(g), 1); win3 = strings(max(g), 1);
for i = 1:max(g)
    q = U(g == i, :); [~, j] = max(q.testLL); win(i) = q.family(j);
    q3 = q(q.nPar <= 3, :); if ~isempty(q3), [~, j] = max(q3.testLL); win3(i) = q3.family(j); end
end
K.winner = win; K.winner3par = win3;
[tf, loc] = ismember(K.unitKey, U.unitKey); for c = {'phase', 'stage', 'group', 'unitType', 'session', 'animal'}, K.(c{1}) = U.(c{1})(loc(tf)); end
writetable(K, fullfile(res, 'dist_winners.csv'));
U = innerjoin(U, K(:, {'unitKey', 'scheme', 'winner'}), 'Keys', {'unitKey', 'scheme'});
[g2, B] = findgroups(U(:, [{'unitType', 'scheme'}, o.by, {'family'}]));
B.n = accumarray(g2, 1);
B.median_testLL = accumarray(g2, U.testLL, [], @(x) median(x, 'omitnan'));
B.median_dLL = accumarray(g2, U.dLL, [], @(x) median(x, 'omitnan'));
B.frac_better = accumarray(g2, U.dLL > 0, [], @(x) mean(x, 'omitnan'));
B.frac_best = accumarray(g2, U.winner == U.family, [], @mean);
B.median_KS = accumarray(g2, U.testKS, [], @(x) median(x, 'omitnan'));
B.median_cover80 = accumarray(g2, U.cover80, [], @(x) median(x, 'omitnan'));
B.median_medAE = accumarray(g2, U.medAE, [], @(x) median(x, 'omitnan'));
B.nPar = arrayfun(@(f) nParOf(char(f)), B.family);
writetable(B, fullfile(res, 'dist_byPhase.csv'));
C = U(U.scheme == "full", :); rows = cell(height(C), 1);
for i = 1:height(C)
    s = jsondecode(C.summary(i)); fn = fieldnames(s);
    r = C(i, {'animal', 'sessionIndex', 'session', 'stage', 'phase', 'group', 'unitType', 'delayLevel', 'nUnit', 'family'});
    for k = 1:numel(fn), v = s.(fn{k}); if isnumeric(v) && isscalar(v), r.(fn{k}) = v; end, end
    rows{i} = r;
end
allNames = {}; for i = 1:numel(rows), allNames = union(allNames, rows{i}.Properties.VariableNames, 'stable'); end
for i = 1:numel(rows), for nm = reshape(setdiff(allNames, rows{i}.Properties.VariableNames), 1, []), rows{i}.(nm{1}) = NaN; end, rows{i} = rows{i}(:, reshape(allNames, 1, [])); end
CP = vertcat(rows{:}); writetable(CP, fullfile(res, 'dist_components.csv'));
S = struct('units', U, 'byPhase', B, 'winners', K, 'components', CP);
end
