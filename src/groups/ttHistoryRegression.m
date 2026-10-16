function [R, N] = ttHistoryRegression(T, outDir, varargin)
%TTHISTORYREGRESSION  Per-session linear regression of the next lick time on the previous trial's lick time and reward, with a shuffle null.
%
%   [R, N] = ttHistoryRegression(T, outDir)
%   [R, N] = ttHistoryRegression(T, outDir, 'nShuffle', 200, 'minTrials', 60)
%
%   T is the cued-trial table with columns session, animal, trialNum, lickTime, rewarded, early and
%   prevEarly (ttComponentResponsibility). Per session and subset the model
%     log tau_t = b0 + b1 log tau_{t-1} + b2 r_{t-1} + b3 log tau_{t-1} r_{t-1}
%   is fitted by least squares on consecutive cued trials (log lick times z-scored within session, so
%   coefficients are in SD units; r_{t-1} is 0/1). Subsets: all trials, trials whose lick is from the
%   timed component, from the early component, and the same split by the previous trial's component.
%   Null: the previous-trial predictors permuted across trials nShuffle times (their joint
%   distribution kept, their link to the next lick broken); the 2.5th and 97.5th percentiles of each
%   coefficient are stored per session (null_lo_*, null_hi_*).
%   R: one row per session x subset with b1, b2, b3, their null bands, R2 and n.
%   Output: <outDir>/results/historyRegression.csv
p = inputParser; p.addParameter('nShuffle', 200, @isnumeric); p.addParameter('minTrials', 60, @isnumeric); p.parse(varargin{:}); o = p.Results;
res = fullfile(outDir, 'results'); if ~exist(res, 'dir'), mkdir(res); end
T = sortrows(T, {'animal', 'session', 'trialNum'}); sess = unique(T.session, 'stable'); rows = cell(numel(sess), 1);
subsets = ["all", "timedNext", "earlyNext", "timedPrev", "earlyPrev"];
parfor s = 1:numel(sess)
    t = T(T.session == sess(s), :); n = height(t); if n < o.minTrials, continue; end
    lt = log(t.lickTime); z = (lt - mean(lt)) / std(lt); r = double(t.rewarded); e = double(t.early);
    y = z(2:end); x1 = z(1:end-1); x2 = r(1:end-1); eNext = e(2:end); ePrev = e(1:end-1); rs = RandStream('mt19937ar', 'Seed', s); out = {};
    for sb = subsets
        switch sb, case "all", m = true(n - 1, 1); case "timedNext", m = eNext == 0; case "earlyNext", m = eNext == 1; case "timedPrev", m = ePrev == 0; case "earlyPrev", m = ePrev == 1; end
        if sum(m) < 30, continue; end
        ym = y(m); x1m = x1(m); x2m = x2(m); b = fitOls(ym, x1m, x2m); nb = nan(o.nShuffle, 3);
        for k = 1:o.nShuffle, pm = randperm(rs, numel(ym)); bp = fitOls(ym, x1m(pm), x2m(pm)); nb(k, :) = bp(2:4)'; end   % the previous-trial predictors are permuted together, breaking their link to the next lick
        yh = [ones(sum(m), 1), x1(m), x2(m), x1(m) .* x2(m)] * b; r2 = 1 - sum((y(m) - yh).^2) / sum((y(m) - mean(y(m))).^2);
        out{end+1} = table(t.animal(1), sess(s), sb, sum(m), b(2), b(3), b(4), prctile(nb(:, 1), 2.5), prctile(nb(:, 1), 97.5), prctile(nb(:, 2), 2.5), prctile(nb(:, 2), 97.5), prctile(nb(:, 3), 2.5), prctile(nb(:, 3), 97.5), r2, ...
            'VariableNames', {'animal', 'session', 'subset', 'n', 'b_prevLick', 'b_prevReward', 'b_interaction', 'null_lo_prevLick', 'null_hi_prevLick', 'null_lo_prevReward', 'null_hi_prevReward', 'null_lo_interaction', 'null_hi_interaction', 'r2'}); %#ok<AGROW>
    end
    if ~isempty(out), rows{s} = vertcat(out{:}); end
end
R = vertcat(rows{:}); N = [];
writetable(R, fullfile(res, 'historyRegression.csv'));
fprintf('ttHistoryRegression: %d rows\n', height(R));
end
function b = fitOls(y, x1, x2)
X = [ones(numel(y), 1), x1(:), x2(:), x1(:) .* x2(:)]; b = X \ y(:);
end
