# timingtask-rl

Reinforcement-learning models of lick timing in the self-timed licking task, with the tools to
choose a lick-time distribution, fit trial-by-trial policy models to recorded sessions, evaluate
them out of sample and compare fitted parameters between groups of animals. MATLAB.

The task (Yang et al., 2025; Majumder et al., 2026): a head-fixed mouse hears a cue and must
withhold licking for an unsignalled delay. The first lick after the cue is rewarded when it comes at
or after the delay. Across training the delay rises from 0.1 s to a target and may later switch
between values in blocks. The models here describe the first-lick time of every trial as a draw
from a policy whose location changes from trial to trial through the previous outcome and through a
learning rule. A complementary gymnasium implementation of the task for training artificial agents
is at [shouvik-majumder/timingtask](https://github.com/shouvik-majumder/timingtask).

## Contents

- `src/data`: session loading, protocol-stage classification, phase assignment, synthetic datasets.
- `src/dist`: candidate first-lick-time distributions and their cross-validated comparison.
- `src/policy`: policy models, forward pass, simulation, maximum-a-posteriori fitting with an
  empirical-Bayes prior, held-out evaluation, single-trial update test, parameter recovery, model
  comparison tables.
- `src/groups`: per-session parameter tables, group trajectories, group comparison with a
  hierarchical bootstrap.
- `src/plots`: session and distribution figures.
- `examples`: `example_simulated.m` (full pipeline on synthetic animals), `example_dataset.m`
  (template for recorded sessions).
- `tests/run_tests.m`: consistency and round-trip tests.
- `docs/models.md`: model definitions.

## Installation

MATLAB R2022a or later with the Optimization Toolbox (`fmincon`, `fminunc`) and the Statistics and
Machine Learning Toolbox (`gamfit`, `wblfit`, `fitgmdist`, `signrank`). The Parallel Computing
Toolbox is used when a pool is open.

```matlab
cd timingtask-rl
ttSetup            % adds src/* to the path
run_tests          % from the tests folder
```

## Data

A dataset is a structure with two tables.

| table | one row per | columns |
|---|---|---|
| `D.trials` | trial | `animal`, `session`, `sessionIndex`, `trialNum`, `lickTime` (s, first lick after the cue), `outcome` (0 early lick, 1 reward, 2 no response, 3 no-cue), `rewarded`, `delay` (s), `cued`, `noCue`, `stage`, `phase` |
| `D.sessions` | session | `animal`, `session`, `sessionIndex`, `date`, `nTrials`, `nCued`, `nRewarded`, `nNoCue`, `delayLevels`, `delayMean`, `stage`, `phase`, `stageIndex`, `stageIndexFromEnd`, `group` |

`ttBuildDataset(root, animals, 'groups', G)` builds it from the laboratory's processed session
files (`<root>/<animal>/<animal>_<yymmdd>_<HHMMSS>.mat`, variable `S`; `ttLoadSession` lists the
fields used). Models are fitted to the cued trials: outcome 0 or 1 with a first lick after the cue.
`ttSimulateDataset` builds a synthetic dataset from a generating model and a training schedule.

Stages are read from the per-trial delays (`ttClassifyStage`): `cueAssociation`, `delayTraining`,
`fixedDelay`, `switching`. Phases (`ttAssignPhases`) divide the stages within each animal:
`cueAssoc`, `delayEarly`, `delayLate`, `fixedPre`, `switching`, `fixedPost`. `stageIndex` and
`stageIndexFromEnd` align sessions to the start or the end of a stage.

## Lick-time distributions

`ttDistFamilies` defines ten candidate distributions with positive support: inverse Gaussian,
log-normal, gamma, Weibull, inverse Gaussian + uniform, log-normal + uniform, log-normal +
exponential (mixture), log-normal racing a constant-rate process, two log-normals, three
log-normals. `ttDistCV(D, outDir)` fits each to every session and every delay block within a
session, in sample, under 5-fold random cross-validation and under a forward split (first 70 % of
trials fitted, last 30 % scored). Scores: held-out log-likelihood per trial, Kolmogorov-Smirnov
distance, coverage of the 10-90 % interval, median absolute error of the predictive median.
`ttDistSummary(outDir)` tabulates the scores relative to a reference family by phase and writes
the fitted component quantities (medians, widths, weights, rates) per unit. `ttDistNoCue`
compares the early component with licking on no-cue trials.

## Policy models

The policy is a density over the first-lick time tau of trial t. Families (`ttPolicy`):

| family | density | channels |
|---|---|---|
| `ln` | LN(z, sigma) | z |
| `lnexp` | (1 - pi) LN(z, sigma) + pi Exp(r) | z, log r, logit pi |
| `lnmix2` | (1 - pi) LN(z, sigma) + pi LN(z1, s1) | z, z1, logit pi |

z is the log median of the timed component. Each channel u_c(t) may be fixed, depend on the
previous lick and reward (`history`: u_c(t) = b_c + w_c' [log tau_{t-1}; r_{t-1}]), be updated by
a learning rule after every trial (`learn`), or both (`historyLearn`). The learning signal is the
reward prediction error delta_t = r_t - rho_t with rho the running reward rate, its absolute value,
a TD error or an action-value error (`signal`). Rules for the timed location (`rule`):

| rule | update of the timed location z after trial t | parameters |
|---|---|---|
| `pg` policy-gradient actor | z <- z + alpha delta_t g_t, g_t = d log p(tau_t) / dz | alpha |
| `pg` with `signal` `unsigned`, unsigned-error actor | z <- z + alpha abs(delta_t) g_t | alpha |
| `threshold` threshold-tracking actor | B = z - m; rewarded and log tau_t < B: B <- B + alpha (log tau_t - B); unrewarded and log tau_t > B: B <- B + alpha (log tau_t - B) | alpha, optional m |
| `stepq` | step-shaped action value in log tau with a delta rule on its location | alpha |
| `bayes` | grid posterior over the delay from the censored outcomes, with forgetting | lambda_b, eps_b |

`ttModel` assembles a model from these options and returns the parameter names and bounds;
`ttModelSet('core')` lists the reference set (static, history, and the three actors with history).
`docs/models.md` gives the full definitions, including the post-reward state, the width channel
and the two-lag history.

## Fitting and evaluation

| function | purpose | output |
|---|---|---|
| `ttFitSession(B, model, prior)` | MAP fit of one session from several starts (`fmincon`) | parameters, log-likelihood |
| `ttFitAnimal(D, animal, model, 'outDir', d)` | empirical-Bayes prior over the animal's sessions, then a MAP fit per session and a forward pass | `<d>/<animal>_<variant>.mat` with per-session parameters and per-trial latent variables |
| `ttFitAll(D, variants, outDir)` | every model of a set for every animal; existing fits are skipped | `<outDir>/fits/` |
| `ttHeldOut(D, outDir)` | forward split per session: parameters from the first 70 % of trials, scores on the rest, with benchmarks (training median, running median of 20 trials), a simulated ceiling and an optional shuffled-order null | `results/heldOut_sessions.csv` |
| `ttCompare(outDir)` | in-sample and held-out summaries by phase, paired contrasts against a reference model | `results/compare_*.csv` |
| `ttUpdateTest(outDir)` | correlation between the model's predicted shift of the timed location after trial t and where the next lick falls, against a within-session permutation null | `results/updateTest_*.csv` |
| `ttRecovery(D, outDir)` | simulate from known parameters on real delay schedules, refit, report true-versus-fitted correlations | `results/recovery_*.csv` |
| `ttGradCheck` | analytic scores against finite differences, normalisation, simulate-then-forward consistency | pass/fail |

Held-out log-likelihood per trial is the primary score: a difference of 0.1 between two models
means every held-out trial is on average 10 % more probable under one than the other. Point
prediction of single lick times (R2 in log seconds, median absolute error) is reported with the
model's own simulated ceiling and the running-median benchmark.

## Group comparison

`ttParamTable(outDir, variant)` collects the fitted parameters of one model per session with the
session's stage, phase, group and alignment indices. `ttGroupTrajectories` gives the per-group
trajectory of a parameter against `stageIndex` or `stageIndexFromEnd` (median over animals with
the interquartile range). `ttGroupCompare(P, measure, groupA, groupB, 'by', 'phase')` reports the
difference of per-animal medians between two groups with a hierarchical bootstrap interval
(animals, then sessions) and a permutation p-value over animal labels.

## Example

```matlab
root = ttSetup(); outDir = fullfile(root, 'output', 'simulated');
D = ttSimulateDataset('nAnimals', 4, 'seed', 1);            % two groups, B learns at half the rate
ttDistCV(D, outDir, 'families', {'lognormal', 'ln_exp', 'lnmix2'}); ttDistSummary(outDir);
ttFitAll(D, ttModelSet('core'), outDir);
ttHeldOut(D, outDir); C = ttCompare(outDir); ttUpdateTest(outDir);
P = ttParamTable(outDir, 'lnexp_ThistoryLearn_Efixed_Rthreshold');
G = ttGroupCompare(P, 'alpha_up', 'A', 'B', 'by', 'stage', 'transform', @log10);
```

`examples/example_simulated.m` runs this sequence; `examples/example_dataset.m` is the template
for a folder of recorded sessions. `data/example_session.mat` is one de-identified recorded
session in the laboratory's format.

## References

Yang, Z., Inagaki, M., Gerfen, C. R., Fontolan, L. and Inagaki, H. K. (2025). Integrator dynamics
in the cortico-basal ganglia loop for flexible motor timing. *Nature* 649, 1244-1253.
https://doi.org/10.1038/s41586-025-09778-2

Majumder, S., Hirokawa, K., Yang, Z., Jain, A., Paletzki, R., Gerfen, C. R., Fontolan, L.,
Romani, S., Yasuda, R. and Inagaki, H. K. (2026). Complementary roles of cell-type-specific
plasticity in shaping neocortical dynamics for learning action timing. *Nature Communications* 17,
8353. https://doi.org/10.1038/s41467-026-74869-1

Roshgadol, J., Chouinard, J. A., Majumder, S., Scott, E. L., et al. Sensitive dLight3 for imaging
broad-spectrum dopamine events across brain regions. Preprint, PMC12393496.

Williams, R. J. (1992). Simple statistical gradient-following algorithms for connectionist
reinforcement learning. *Machine Learning* 8, 229-256. https://doi.org/10.1007/BF00992696

Pearce, J. M. and Hall, G. (1980). A model for Pavlovian learning: variations in the effectiveness
of conditioned but not of unconditioned stimuli. *Psychological Review* 87, 532-552.
https://doi.org/10.1037/0033-295X.87.6.532

## License

MIT. See `LICENSE`.
