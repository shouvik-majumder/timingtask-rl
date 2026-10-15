%EXAMPLE_SESSION  Fit the core models to the single recorded session shipped with the repository.
%
%   The session (data/example_session.mat, variable S) is a fixed-delay session. Each model is
%   fitted by maximum a posteriori under the default prior (no empirical Bayes, since there is one
%   session), and scored on the last 30 % of trials after fitting on the first 70 %.
root = ttSetup();
L = load(fullfile(root, 'data', 'example_session.mat')); T = ttLoadSession(L.S, 'animal', 'M01', 'session', 'M01_s01');
q = T(T.cued, :); B = struct('tau', q.lickTime, 'R', double(q.rewarded), 'delay', q.delay);
fprintf('%d cued trials, delay %.2f s, reward rate %.2f\n', height(q), median(q.delay), mean(q.rewarded));
n = numel(B.tau); mask = false(n, 1); mask(1:round(0.7 * n)) = true;
V = ttModelSet('core'); rows = {};
for v = 1:numel(V)
    model = ttModel(V{v}{:}); prior = ttPrior(model);
    theta = ttFitSession(B, model, prior, 'nStarts', 3, 'llMask', mask);
    fw = ttForward(theta, B, model, struct('predict', true));
    rows{end+1} = table(string(model.shortName), model.nParam, mean(fw.loglik(mask)), mean(fw.loglik(~mask)), 'VariableNames', {'model', 'nParam', 'll_train', 'll_heldout'}); %#ok<AGROW>
end
disp(vertcat(rows{:}));
