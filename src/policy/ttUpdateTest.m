function [U, S] = ttUpdateTest(outDir, varargin)
%TTUPDATETEST  Does a model's update after trial t predict where the next lick falls?
%
%   [U, S] = ttUpdateTest(outDir)                          every fitted learning model in <outDir>/fits
%   [U, S] = ttUpdateTest(outDir, 'variants', {...})
%
%   For every session of a fitted model (whole-session fit), with u_z(t) the timed log median used on
%   trial t and step_z(t) the learning step taken after trial t:
%     target(t)  = log tau_{t+1} - u_z(t)      where the next lick falls relative to the current location
%     shift(t)   = u_z(t+1) - u_z(t)           the model's predicted shift (learning step plus the history term)
%   restricted to trials whose next lick is attributed to the timed component (early-component
%   responsibility below 0.5). Per session: Spearman correlation of target with the shift, with the
%   learning step alone and with delta alone; regression slope of target on the shift; the 97.5th
%   percentile of the correlation after permuting the shifts within the session (null).
%   U: one row per session x variant. S: summary per variant x phase.
%   Output: <outDir>/results/updateTest_sessions.csv, updateTest_byPhase.csv
p = inputParser; p.addParameter('variants', {}, @iscell); p.addParameter('minTrials', 60, @isnumeric); p.addParameter('nNull', 50, @isnumeric); p.parse(varargin{:}); o = p.Results;
fd = fullfile(outDir, 'fits'); res = fullfile(outDir, 'results'); if ~exist(res, 'dir'), mkdir(res); end
files = dir(fullfile(fd, '*.mat')); rows = {};
for f = 1:numel(files)
    L = load(fullfile(fd, files(f).name), 'Res'); Tr = L.Res.trials; P = L.Res.params; v = L.Res.model.shortName;
    if ~isempty(o.variants) && ~ismember(v, o.variants), continue; end
    if ~any(contains(Tr.Properties.VariableNames, 'step_z')), continue; end
    for s = 1:height(P)
        q = Tr(Tr.session == P.session(s), :); n = height(q); if n < o.minTrials, continue; end
        uz = q.u_z; lt = log(q.tau); st = q.step_z; dl = q.delta; rs = q.resp;
        t = (L.Res.model.warmup + 1:n - 1)'; target = lt(t + 1) - uz(t); pred = uz(t + 1) - uz(t); step = st(t); del = dl(t);
        keep = rs(t + 1) < 0.5 & isfinite(target) & isfinite(pred); target = target(keep); pred = pred(keep); step = step(keep); del = del(keep);
        if numel(target) < 30 || std(pred) == 0, continue; end
        rPred = corr(pred, target, 'Type', 'Spearman'); rDelta = corr(del, target, 'Type', 'Spearman');
        if std(step, 'omitnan') > 0, rStep = corr(step, target, 'Type', 'Spearman', 'Rows', 'complete'); else, rStep = NaN; end
        b = [ones(size(pred)), pred] \ target; slope = b(2);
        rsm = RandStream('mt19937ar', 'Seed', s); nullR = nan(o.nNull, 1);
        for k = 1:o.nNull, nullR(k) = corr(pred(randperm(rsm, numel(pred))), target, 'Type', 'Spearman'); end
        rows{end+1} = table(string(v), P.animal(s), P.sessionIndex(s), P.session(s), P.phase(s), P.group(s), numel(target), rPred, rStep, rDelta, slope, std(pred), mean(abs(step) > 0, 'omitnan'), prctile(nullR, 97.5), ...
            'VariableNames', {'variant', 'animal', 'sessionIndex', 'session', 'phase', 'group', 'nTrials', 'r_shift', 'r_step', 'r_delta', 'slope', 'sd_shift', 'frac_nonzero_step', 'null97'}); %#ok<AGROW>
    end
end
U = vertcat(rows{:}); writetable(U, fullfile(res, 'updateTest_sessions.csv'));
[g, S] = findgroups(U(:, {'variant', 'phase'}));
S.n = accumarray(g, 1); S.median_r_shift = accumarray(g, U.r_shift, [], @median); S.frac_r_shift_pos = accumarray(g, U.r_shift > 0, [], @mean); S.frac_above_null = accumarray(g, U.r_shift > U.null97, [], @mean);
S.median_r_step = accumarray(g, U.r_step, [], @(x) median(x, 'omitnan')); S.frac_learning_active = accumarray(g, U.frac_nonzero_step > 0, [], @mean); S.median_r_delta = accumarray(g, U.r_delta, [], @median); S.median_slope = accumarray(g, U.slope, [], @median);
writetable(S, fullfile(res, 'updateTest_byPhase.csv'));
end
