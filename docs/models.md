# Model definitions

Notation: trial t, first-lick time tau_t (s), reward r_t in {0, 1}, delay d_t (s). A lick is
rewarded when tau_t >= d_t. Logarithms are natural.

## Policy families

| family | p(tau) | channels u |
|---|---|---|
| `ln` | LN(tau; z, sigma) | [z] |
| `lnexp` | (1 - pi) LN(tau; z, sigma) + pi r e^{-r tau} | [z; log r; ell], pi = 1 / (1 + e^{-ell}) |
| `lnmix2` | (1 - pi) LN(tau; z, sigma) + pi LN(tau; z1, s1) | [z; z1; ell] |

LN(tau; z, sigma) = (tau sigma sqrt(2 pi))^{-1} exp(-(log tau - z)^2 / (2 sigma^2)). z is the log
median of the timed component. The score g = d log p / du is, for a mixture, the responsibility of
the component times the component's own score, and for ell the responsibility minus pi
(`ttPolicy`; verified against finite differences by `ttGradCheck`).

## Trial-to-trial dependence of a channel

With h_t = [log tau_{t-1}; r_{t-1}] (or two lags) and x_t = [1; h_t]:

| mode | u_c(t) | updated after trial t |
|---|---|---|
| `fixed` | b_c | nothing |
| `history` | b_c + w_c' h_t | nothing |
| `learn` | b_c(t) | b_c by the rule |
| `historyLearn` | W_c(t)' x_t | the intercept (rules other than `pg`) or all of W_c (`pg`) |

The first `warmup` trials (10) use x_t = [1; 0] and no update. The early channels (log r, ell)
share one learning rate and always use the policy gradient.

## Learning signal

| `signal` | delta_t | state |
|---|---|---|
| `rate` | r_t - rho_t | rho_t = (1 - beta) rho_{t-1} + beta r_t, beta = 0.02 |
| `unsigned` | abs(r_t - rho_t) | as above |
| `td` | r_t - rho_t + gamma V - V | V <- V + alpha_v delta_t |
| `q` | r_t - Q(tau_t) | Q(tau) = v' [1; log tau; log tau^2], v <- v + alpha_q delta_t phi_t |
| `raw` | r_t | |

## Learning rule for the timed location

| `rule` | update | parameters |
|---|---|---|
| `pg` policy-gradient actor | W_z <- W_z + alpha delta_t g_z(t) x_t | alpha_z |
| `pg` with `unsigned` unsigned-error actor | the same with abs(delta_t) | alpha_z |
| `threshold` threshold-tracking actor | B_t = b_z(t) - m. If r_t = 1 and log tau_t < B_t: b_z <- b_z + alpha_down (log tau_t - B_t). If r_t = 0 and log tau_t > B_t: b_z <- b_z + alpha_up (log tau_t - B_t). | alpha_up, alpha_down (equal when `tieRates`), m (when `margin`) |
| `stepq` | Q(tau) = s((log tau - B) / kappa), s the logistic function, kappa = 0.15; delta = r_t - Q(tau_t); B <- B + alpha delta dQ/dB | alpha_z, m |
| `bayes` | log posterior over B on a grid in [log 0.03, log 8]: lambda_b log post + log lik, lik = eps_b + (1 - 2 eps_b) [B <= log tau_t] when rewarded, [B > log tau_t] when not; z = posterior mean + m | lambda_b, eps_b, m |
| `shift` | b_z <- b_z - alpha delta_t | alpha_z (either sign) |
| `qmax` | b_z <- b_z + alpha (argmax Q - b_z) with `signal` `q` | alpha_z |

The policy-gradient actor is REINFORCE with a reward-rate baseline (Williams, 1992). The
unsigned-error actor weights the same gradient by the magnitude of the prediction error, in the
spirit of the Pearce-Hall associability. The threshold-tracking actor keeps a point estimate of the
delay and moves it only on the two trial types that bound it: a reward obtained before the estimate
lowers it, an omission after the estimate raises it.

## Optional terms

- `wRule`: the history weights of the timed channel learn, W_h <- W_h + alpha_w delta_t g_z(t) h_t
  (`pg`) or with abs(delta_t) (`unsigned`).
- `interaction`: w_z_tauR log tau_{t-1} r_{t-1} added to the timed location.
- `width`: log sigma_t = log sigma + w_lsig' h_t (`history`) or log sigma + c_ssig s_t (`state`).
- `state`: s_t = lambda_s s_{t-1} + r_{t-1}, added to the timed location (c_sz s_t) and the early
  logit (c_spi s_t).

## Fitting

Per session, the negative log posterior -sum_t log p(tau_t) + 0.5 sum_i w_i (theta_i - mu_i)^2 is
minimised with `fmincon` (SQP) within the parameter bounds from several starts. Per animal, the
prior (mu, w) is set by empirical Bayes: repeated single-start fits of every session, the prior mean
moved towards the across-session mean and the precision towards shrink / variance with damping,
precisions clipped to [0.5, 100]. The final fit uses three starts under the final prior.

## Held-out evaluation

Parameters are fitted on the first 70 % of a session's trials. The forward pass then runs through
the whole session with those parameters (the model's state is driven by the observed trials) and
the last 30 % are scored: log-likelihood per trial, R2 of the predicted mean log lick time, median
absolute error of the predicted median, coverage of the 10-90 % interval and the KS distance of the
probability integral transform. Benchmarks: the training median and the running median of the
previous 20 licks, each with a log-normal of the training width. Ceiling: the same scores on
sessions simulated from the training fit on the session's delay schedule.
