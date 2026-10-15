function ok = run_tests()
%RUN_TESTS  Consistency and round-trip tests of the library.
%
%   ok = run_tests()
%
%   1. ttGradCheck: scores, normalisation, CDF, simulate-then-forward consistency.
%   2. ttClassifyStage on constructed delay sequences.
%   3. ttLoadSession on a constructed session structure.
%   4. Round trip: a session simulated by the threshold-tracking actor is refitted and the learning
%      rate and location are recovered within tolerance.
%   5. ttDistFamilies: every family fits a log-normal + exponential sample and integrates to one.
root = ttSetup(); ok = true;
ok = check(ttGradCheck(false), 'ttGradCheck') && ok;
ok = check(ttClassifyStage(0.1 * ones(200, 1), 200) == "cueAssociation", 'stage: cue association') && ok;
ok = check(ttClassifyStage(1.8 * ones(800, 1), 800) == "fixedDelay", 'stage: fixed delay') && ok;
ok = check(ttClassifyStage([0.5 * ones(300, 1); 0.6 * ones(300, 1); 0.7 * ones(300, 1)], 900) == "delayTraining", 'stage: delay training') && ok;
ok = check(ttClassifyStage([ones(400, 1); 3 * ones(400, 1)], 800) == "switching", 'stage: switching') && ok;
S = struct('numTrial', 5, 'Trial', struct('firstLickTime', [0.5 NaN 0.2 1.1 NaN], 'ResponseTypeVector', [1 2 0 1 3]), 'TrialSetting', struct('DelayPeriod', [0.4 0.4 0.4 0.4 0]));
T = ttLoadSession(S, 'animal', 'X', 'session', 'X_s1');
ok = check(isequal(T.cued', [true false true true false]) && isequal(T.noCue', [false false false false true]) && isequal(T.rewarded', [true false false true false]), 'ttLoadSession') && ok;
model = ttModel('family', 'lnexp', 'timed', 'historyLearn', 'early', 'fixed', 'rule', 'threshold'); prior = ttPrior(model);
th = prior.mu; th(model.idx.b_z) = log(0.4); th(model.idx.sigma) = 0.25; th(model.idx.b_logr) = log(2); th(model.idx.b_ell) = log(0.1 / 0.9); th(model.idx.w_z_tau1) = 0.2; th(model.idx.w_z_R1) = -0.15; th(model.idx.alpha_up) = 0.08;
delay = reshape(repmat(linspace(0.4, 1.2, 8), 100, 1), [], 1); Ssim = ttSimulate(th, struct('delay', delay), model, 11);
thF = ttFitSession(struct('tau', Ssim.tau, 'R', Ssim.R), model, prior, 'nStarts', 3);
ok = check(abs(log(thF(model.idx.alpha_up)) - log(0.08)) < log(2.5) && abs(thF(model.idx.sigma) - 0.25) < 0.08, sprintf('round trip: alpha_up %.3f (true 0.08), sigma %.3f (true 0.25)', thF(model.idx.alpha_up), thF(model.idx.sigma))) && ok;
rs = RandStream('mt19937ar', 'Seed', 2); x = [exp(log(0.8) + 0.3 * randn(rs, 400, 1)); -log(rand(rs, 100, 1)) / 3];
fam = ttDistFamilies(); good = true;
for k = 1:numel(fam)
    try, thk = fam(k).fit(x, max(x) * 1.001); Z = integral(@(z) exp(fam(k).logpdf(z, thk, max(x) * 1.001)), 1e-6, 500); good = good && abs(Z - 1) < 0.02; catch, good = false; end
end
ok = check(good, 'ttDistFamilies: fits and normalisation') && ok;
if ok, fprintf('run_tests: all passed\n'); else, fprintf('run_tests: FAILURES\n'); end
end
function ok = check(cond, name)
ok = logical(cond); if ok, fprintf('  pass  %s\n', name); else, fprintf('  FAIL  %s\n', name); end
end
