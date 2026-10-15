function [R, S] = ttRecovery(D, outDir, varargin)
%TTRECOVERY  Parameter recovery: simulate sessions from known parameters on real delay schedules, refit, compare.
%
%   [R, S] = ttRecovery(D, outDir)
%   [R, S] = ttRecovery(D, outDir, 'variants', {...}, 'nSess', 24, 'seed', 3, 'draws', drawStruct)
%
%   For each model, nSess sessions are simulated with ttSimulate on the delay schedules of real
%   sessions drawn from D (half from delay-training sessions, half from the other stages, each with
%   at least 300 cued trials). Static parameters come from the session's lick times (moment
%   values); rule-specific parameters are drawn from the ranges in 'draws' (a struct of function
%   handles returning one value; defaults cover every parameter name). Each simulated session is
%   refitted with ttFitSession (3 starts, default prior).
%   R: one row per simulated session with true_<param> and fit_<param> columns and the fitted
%   log-likelihood minus the log-likelihood at the true parameters. S: per model and parameter,
%   Spearman correlation between true and fitted values and the median relative error.
%   Output: <outDir>/results/recovery_<variant>.csv, recovery_summary.csv
p = inputParser; p.addParameter('variants', ttModelSet('core'), @iscell); p.addParameter('nSess', 24, @isnumeric); p.addParameter('seed', 3, @isnumeric); p.addParameter('draws', struct(), @isstruct); p.addParameter('minTrials', 300, @isnumeric); p.parse(varargin{:}); o = p.Results;
res = fullfile(outDir, 'results'); if ~exist(res, 'dir'), mkdir(res); end
T = D.trials(D.trials.cued, :); Sx = D.sessions(D.sessions.nCued >= o.minTrials, :);
rs = RandStream('mt19937ar', 'Seed', o.seed);
dl = Sx(Sx.stage == "delayTraining", :); ex = Sx(Sx.stage ~= "delayTraining", :);
nHalf = floor(o.nSess / 2); pick = [dl(randperm(rs, height(dl), min(nHalf, height(dl))), :); ex(randperm(rs, height(ex), min(o.nSess - nHalf, height(ex))), :)];
draw = struct('alpha_z', @() 10^(-3.5 + 1.5 * rand(rs)), 'alpha_early', @() 10^(-3.5 + 1.5 * rand(rs)), 'alpha_up', @() 0.02 + 0.3 * rand(rs), 'alpha_down', @() 0.02 + 0.3 * rand(rs), 'm_z', @() -0.1 + 0.5 * rand(rs), ...
    'alpha_q', @() 0.02 + 0.2 * rand(rs), 'alpha_v', @() 0.02 + 0.2 * rand(rs), 'gamma', @() 0.5 + 0.4 * rand(rs), 'lambda_s', @() 0.2 + 0.7 * rand(rs), 'c_sz', @() -0.4 + 0.4 * rand(rs), 'c_spi', @() 0.8 * rand(rs), ...
    'w_z_tau1', @() 0.4 * rand(rs), 'w_z_R1', @() -0.4 * rand(rs), 'w_z_tau2', @() 0.2 * rand(rs), 'w_z_R2', @() -0.2 * rand(rs), 'alpha_w', @() 10^(-3 + 1.5 * rand(rs)), 'w_lsig_tau1', @() -0.3 + 0.6 * rand(rs), 'w_lsig_R1', @() -0.3 + 0.6 * rand(rs), ...
    'w_z_tauR', @() -0.4 + 0.8 * rand(rs), 'lambda_b', @() 0.7 + 0.29 * rand(rs), 'eps_b', @() 0.02 + 0.2 * rand(rs));
for f = fieldnames(o.draws)', draw.(f{1}) = o.draws.(f{1}); end
rows = {};
for v = 1:numel(o.variants)
    model = ttModel(o.variants{v}{:}); prior = ttPrior(model); fprintf('ttRecovery: %s\n', model.shortName);
    thAll = cell(height(pick), 1);
    for s = 1:height(pick)   % draws are made serially so that the sequence is reproducible
        q = sortrows(T(T.session == pick.session(s), :), 'trialNum'); lt = log(q.lickTime);
        th = prior.mu; th(model.idx.b_z) = median(lt(lt >= median(lt))); th(model.idx.sigma) = max(0.7 * std(lt), 0.1);
        if isfield(model.idx, 'b_logr'), th(model.idx.b_logr) = log(1 / mean(q.lickTime)); end; if isfield(model.idx, 'b_ell'), th(model.idx.b_ell) = log(0.15 / 0.85); end
        if isfield(model.idx, 'b_z1'), th(model.idx.b_z1) = log(max(prctile(q.lickTime, 15), 0.03)); end
        if pick.stage(s) == "delayTraining", th(model.idx.b_z) = log(0.5 * q.delay(1) + 0.15); end
        for ip = 1:model.nParam, pn = model.names{ip}; if isfield(draw, pn), th(ip) = draw.(pn)(); end, end
        thAll{s} = min(max(th, model.lb(:)), model.ub(:));
    end
    out = cell(height(pick), 1);
    parfor s = 1:height(pick)
        q = sortrows(T(T.session == pick.session(s), :), 'trialNum'); B = struct('tau', q.lickTime, 'R', double(q.rewarded), 'delay', q.delay); th = thAll{s};
        Ssim = ttSimulate(th, B, model, 100 + s); Bsim = struct('tau', Ssim.tau, 'R', Ssim.R);
        [thF, fit] = ttFitSession(Bsim, model, prior, 'nStarts', 3);
        llTrue = ttForward(th, Bsim, model).totalLogLik;
        r = table(string(model.shortName), pick.session(s), pick.stage(s), numel(Bsim.tau), fit.logLik - llTrue, 'VariableNames', {'variant', 'session', 'stage', 'n', 'llFit_minus_llTrue'});
        for i = 1:model.nParam, r.(['true_', model.names{i}]) = th(i); r.(['fit_', model.names{i}]) = thF(i); end
        out{s} = r;
    end
    rows{end+1} = vertcat(out{:}); %#ok<AGROW>
end
S = {};
for v = 1:numel(rows)
    q = rows{v}; nm = q.Properties.VariableNames; tn = nm(startsWith(nm, 'true_'));
    for k = 1:numel(tn)
        pn = extractAfter(tn{k}, 'true_'); x = q.(tn{k}); y = q.(['fit_', pn]); if std(x) == 0, continue; end
        S{end+1} = table(q.variant(1), string(pn), corr(x, y, 'Type', 'Spearman'), median(abs(y - x) ./ max(abs(x), 1e-3)), median(q.llFit_minus_llTrue), 'VariableNames', {'variant', 'param', 'spearman_true_fit', 'median_rel_error', 'median_llFit_minus_llTrue'}); %#ok<AGROW>
    end
end
S = vertcat(S{:}); R = rows;
for v = 1:numel(rows), writetable(rows{v}, fullfile(res, ['recovery_', char(rows{v}.variant(1)), '.csv'])); end; writetable(S, fullfile(res, 'recovery_summary.csv'));
end
