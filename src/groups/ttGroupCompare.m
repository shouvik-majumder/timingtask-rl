function G = ttGroupCompare(P, measure, groupA, groupB, varargin)
%TTGROUPCOMPARE  Difference of a per-session measure between two groups of animals, with animal as the unit.
%
%   G = ttGroupCompare(P, 'alpha_up', 'Het', 'WT')
%   G = ttGroupCompare(P, measure, groupA, groupB, 'by', 'phase', 'nBoot', 2000, 'transform', @log10)
%
%   P is a table with columns animal, group, the measure, and the column named by 'by' (default
%   phase). For every level of 'by': the per-animal median of the measure, the group medians of
%   those per-animal values, their difference (A minus B), a hierarchical bootstrap interval
%   (resample animals within group, then sessions within animal; nBoot replicates) and a
%   permutation p-value (animal labels shuffled between groups, nPerm replicates). Sessions with a
%   non-finite measure are dropped.
p = inputParser; p.addParameter('by', 'phase', @ischar); p.addParameter('nBoot', 2000, @isnumeric); p.addParameter('nPerm', 2000, @isnumeric); p.addParameter('transform', @(x) x, @(x) isa(x, 'function_handle')); p.addParameter('seed', 1, @isnumeric); p.parse(varargin{:}); o = p.Results;
rs = RandStream('mt19937ar', 'Seed', o.seed); rows = {};
P = P(isfinite(P.(measure)), :); P.value = o.transform(P.(measure));
for lev = unique(P.(o.by), 'stable')'
    q = P(P.(o.by) == lev, :); A = q(q.group == string(groupA), :); B = q(q.group == string(groupB), :);
    [va, ia] = perAnimal(A); [vb, ib] = perAnimal(B); if numel(va) < 2 || numel(vb) < 2, continue; end
    d0 = median(va) - median(vb);
    db = nan(o.nBoot, 1);
    for k = 1:o.nBoot, db(k) = median(bootAnimals(ia, rs)) - median(bootAnimals(ib, rs)); end
    allV = [va; vb]; nA = numel(va); dp = nan(o.nPerm, 1);
    for k = 1:o.nPerm, pr = randperm(rs, numel(allV)); dp(k) = median(allV(pr(1:nA))) - median(allV(pr(nA+1:end))); end
    pPerm = mean(abs(dp) >= abs(d0));
    rows{end+1} = table(string(measure), string(lev), string(groupA), string(groupB), numel(va), numel(vb), height(A), height(B), median(va), median(vb), d0, prctile(db, 2.5), prctile(db, 97.5), pPerm, ...
        'VariableNames', {'measure', 'level', 'groupA', 'groupB', 'nAnimalsA', 'nAnimalsB', 'nSessionsA', 'nSessionsB', 'medianA', 'medianB', 'diff', 'ci_lo', 'ci_hi', 'p_perm'}); %#ok<AGROW>
end
if isempty(rows)
    G = table('Size', [0 14], 'VariableTypes', [repmat({'string'}, 1, 4), repmat({'double'}, 1, 10)], 'VariableNames', {'measure', 'level', 'groupA', 'groupB', 'nAnimalsA', 'nAnimalsB', 'nSessionsA', 'nSessionsB', 'medianA', 'medianB', 'diff', 'ci_lo', 'ci_hi', 'p_perm'});
else
    G = vertcat(rows{:});
end
end
function [v, cells] = perAnimal(Q)
an = unique(Q.animal); v = nan(numel(an), 1); cells = cell(numel(an), 1);
for i = 1:numel(an), x = Q.value(Q.animal == an(i)); v(i) = median(x); cells{i} = x; end
end
function v = bootAnimals(cells, rs)
n = numel(cells); pick = randi(rs, n, n, 1); v = nan(n, 1);
for i = 1:n, x = cells{pick(i)}; v(i) = median(x(randi(rs, numel(x), numel(x), 1))); end
end
