function Res = ttFitAnimal(D, animal, model, varargin)
%TTFITANIMAL  Fit one model to every session of one animal: empirical-Bayes prior, then a MAP fit per session.
%
%   Res = ttFitAnimal(D, 'RD50', ttModel('family', 'lnexp', 'timed', 'historyLearn', 'early', 'fixed'))
%   Res = ttFitAnimal(D, animal, model, 'ebIter', 8, 'nStarts', 3, 'sessions', {...}, 'outDir', ..., 'save', true)
%
%   Empirical Bayes: repeat { MAP-fit every session with one start under the current prior; move
%   the prior mean towards the across-session mean and the precision towards shrink / variance
%   (clipped to prior.wMin..wMax) with damping eta }. Then a final MAP fit with nStarts starts under
%   the final prior and a forward pass for the trial-by-trial latent variables. Sessions run in
%   parallel when a parallel pool is open.
%
%   Res.params  one row per session: parameters, logLik, llPerTrial, BIC, AIC, nTrials, exitflag,
%               plus the session columns of D.sessions
%   Res.trials  one row per fitted trial: session, trialNum, tau, R, delay, loglik, timedMedian,
%               u_<channel>, step_<channel>, resp, delta, rho, predMedian, predMeanLog, pit
%   Res.prior, Res.model, Res.eb (prior history), Res.fitOptions
%   Saved to <outDir>/<animal>_<model.shortName>.mat when 'save' is true.
p = inputParser;
p.addParameter('ebIter', 8, @isnumeric); p.addParameter('nStarts', 3, @isnumeric); p.addParameter('sessions', {}, @(x) iscell(x) || isstring(x));
p.addParameter('outDir', '', @ischar); p.addParameter('save', true, @islogical); p.addParameter('minTrials', 50, @isnumeric);
p.addParameter('eta', 0.3, @isnumeric); p.addParameter('shrink', 3, @isnumeric); p.addParameter('verbose', true, @islogical);
p.parse(varargin{:}); o = p.Results;
if o.save && isempty(o.outDir), error('ttFitAnimal:outDir', 'outDir is required when save is true'); end
if o.save && ~exist(o.outDir, 'dir'), mkdir(o.outDir); end
[Bs, S] = ttSessionBlocks(D, animal, 'sessions', o.sessions, 'minTrials', o.minTrials);
nS = numel(Bs); if nS < 2, error('ttFitAnimal:sessions', 'fewer than 2 sessions for %s', animal); end
prior = ttPrior(model); nP = model.nParam; eb = struct('mu', {{}}, 'w', {{}});
if o.verbose, fprintf('ttFitAnimal %s %s: %d sessions, %d parameters\n', animal, model.shortName, nS, nP); end
for it = 1:o.ebIter
    th = nan(nS, nP);
    parfor s = 1:nS, th(s, :) = ttFitSession(Bs{s}, model, prior, 'nStarts', 1)'; end
    ok = all(isfinite(th), 2); if sum(ok) < 2, break; end
    empMu = mean(th(ok, :), 1)'; empVar = max(var(th(ok, :), 0, 1)', 1e-4);
    rawW = min(max(o.shrink ./ empVar, prior.wMin), prior.wMax);
    newMu = (1 - o.eta) * prior.mu + o.eta * empMu; newW = (1 - o.eta) * prior.w + o.eta * rawW;
    newMu = min(max(newMu, model.lb(:)), model.ub(:));
    dW = norm(newW - prior.w); prior.mu = newMu; prior.w = newW; eb.mu{end+1} = newMu; eb.w{end+1} = newW;
    if o.verbose, fprintf('  EB %2d: |dw| %.3f\n', it, dW); end
    if dW < 1e-3, break; end
end
prior.source = sprintf('empirical Bayes, %d iterations, %d sessions', numel(eb.mu), nS);
th = nan(nS, nP); ll = nan(nS, 1); ef = nan(nS, 1); trs = cell(nS, 1);
parfor s = 1:nS
    [t, fit] = ttFitSession(Bs{s}, model, prior, 'nStarts', o.nStarts); th(s, :) = t'; ll(s) = fit.logLik; ef(s) = fit.exitflag;
    trs{s} = trialTable(t, Bs{s}, model);
end
n = cellfun(@(b) numel(b.tau), Bs)';
P = [S, array2table(th, 'VariableNames', model.names)];
P.logLik = ll; P.llPerTrial = ll ./ n; P.BIC = -2 * ll + nP * log(n); P.AIC = -2 * ll + 2 * nP; P.exitflag = ef; P.nFit = n;
Res = struct('animal', string(animal), 'params', P, 'trials', vertcat(trs{:}), 'prior', prior, 'model', model, 'eb', eb, 'fitOptions', o);
if o.save, save(fullfile(o.outDir, sprintf('%s_%s.mat', animal, model.shortName)), 'Res'); end
if o.verbose, fprintf('ttFitAnimal %s %s done: median log-likelihood per trial %.4f\n', animal, model.shortName, median(P.llPerTrial)); end
end
function tr = trialTable(theta, B, model)
out = ttForward(theta, B, model, struct('predict', true)); n = numel(B.tau);
tr = table(repmat(string(B.session), n, 1), repmat(B.sessionIndex, n, 1), B.trialNum(:), B.tau(:), B.R(:), B.delay(:), out.loglik, out.timedMedian, out.resp, out.delta, out.rho, out.predMedian, out.predMeanLog, out.pit, ...
    'VariableNames', {'session', 'sessionIndex', 'trialNum', 'tau', 'R', 'delay', 'loglik', 'timedMedian', 'resp', 'delta', 'rho', 'predMedian', 'predMeanLog', 'pit'});
for c = 1:numel(model.chan), tr.(['u_', model.chan(c).name]) = out.u(:, c); tr.(['step_', model.chan(c).name]) = out.step(:, c); end
end
