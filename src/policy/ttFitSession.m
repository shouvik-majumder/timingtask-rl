function [theta, fit] = ttFitSession(B, model, prior, varargin)
%TTFITSESSION  Maximum a posteriori fit of a model to one session (fmincon, several starts).
%
%   [theta, fit] = ttFitSession(B, model, prior, 'nStarts', 3, 'llMask', mask, 'start', theta0)
%
%   Minimises -sum(log p(tau_t)) over the trials in llMask (all trials by default) plus
%   0.5 * sum(prior.w .* (theta - prior.mu).^2), within model.lb / model.ub. Starts: a moment-based
%   start from the training lick times, the prior mean, then random points; an optional warm start
%   goes first. fit: negLogPost, logLik (masked), exitflag, allStarts (table), nTrials.
p = inputParser;
p.addParameter('nStarts', 3, @isnumeric);
p.addParameter('llMask', [], @(x) islogical(x) || isempty(x));
p.addParameter('start', [], @isnumeric);
p.addParameter('display', 'none', @ischar);
p.addParameter('maxFunEvals', 6000, @isnumeric);
p.addParameter('seed', 0, @isnumeric);
p.parse(varargin{:}); o = p.Results;
lb = model.lb(:); ub = model.ub(:); nP = model.nParam;
opt = struct('llMask', o.llMask);
obj = @(th) negLogPost(th, B, model, prior, opt);
opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', o.display, 'ScaleProblem', true, ...
    'MaxFunctionEvaluations', o.maxFunEvals, 'OptimalityTolerance', 1e-4, 'StepTolerance', 1e-6);
tau = B.tau(:); if ~isempty(o.llMask), tau = tau(o.llMask); end; lt = log(tau);
thM = prior.mu(:);
for i = 1:nP
    nm = model.names{i};
    switch nm
        case 'b_z',    thM(i) = median(lt(lt >= median(lt)));
        case 'b_z1',   thM(i) = log(max(prctile(tau, 15), 0.03));
        case 'b_logr', thM(i) = log(1 / mean(tau));
        case 'b_ell',  thM(i) = log(0.2 / 0.8);
        case 'sigma',  thM(i) = max(0.7 * std(lt), 0.05);
        case 's1',     thM(i) = 1.0;
        case {'alpha_up', 'alpha_down', 'alpha_q'}, thM(i) = 0.05;
        case 'm_z', thM(i) = 0.2;
        case 'lambda_s', thM(i) = 0.5;
        case {'c_sz', 'c_spi', 'c_ssig', 'w_z_tauR'}, thM(i) = 0;
        case 'alpha_w', thM(i) = 0.005;
        case 'lambda_b', thM(i) = 0.9;
        case 'eps_b', thM(i) = 0.05;
        otherwise, if startsWith(nm, 'alpha'), thM(i) = 0.005; elseif startsWith(nm, 'w_'), thM(i) = 0; end
    end
end
thM = min(max(thM, lb), ub); thP = min(max(prior.mu(:), lb), ub);
starts = [thM'; thP'];
if o.nStarts > 2
    rs = RandStream('mt19937ar', 'Seed', o.seed); rnd = lb' + rand(rs, o.nStarts - 2, nP) .* (ub' - lb');
    for i = 1:nP, if startsWith(model.names{i}, 'alpha'), rnd(:, i) = 0.05 * rand(rs, o.nStarts - 2, 1); end, end
    starts = [starts; rnd];
end
if o.nStarts == 1, starts = thM'; end
if ~isempty(o.start), st = min(max(o.start(:), lb), ub); starts = [st'; starts(1:max(0, o.nStarts - 1), :)]; end
res = nan(size(starts, 1), nP + 2);
for s = 1:size(starts, 1)
    try
        [th, fval, ef] = fmincon(obj, starts(s, :)', [], [], [], [], lb, ub, [], opts); res(s, :) = [th', fval, ef];
    catch ME
        warning('ttFitSession:start', 'start %d failed: %s', s, ME.message);
    end
end
[~, best] = min(res(:, end - 1)); theta = res(best, 1:nP)';
fwd = ttForward(theta, B, model, opt);
fit = struct('negLogPost', res(best, end - 1), 'logLik', fwd.totalLogLik, 'exitflag', res(best, end), 'nTrials', numel(B.tau), ...
    'allStarts', array2table(res, 'VariableNames', [model.names, {'negLogPost', 'exitflag'}]));
end
function nlp = negLogPost(theta, B, model, prior, opt)
fwd = ttForward(theta, B, model, opt);
d = theta(:) - prior.mu(:);
nlp = -fwd.totalLogLik + 0.5 * sum(prior.w(:) .* d.^2);
if ~isfinite(nlp), nlp = 1e12; end
end
