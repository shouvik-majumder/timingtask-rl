function model = ttModel(varargin)
%TTMODEL  Define a lick-time policy model: its components, how each changes from trial to trial, and the learning rule.
%
%   model = ttModel('family', 'lnexp', 'timed', 'historyLearn', 'early', 'fixed', 'signal', 'rate', 'rule', 'pg')
%
%   Policy families (ttPolicy): 'ln' (timed log-normal), 'lnexp' (timed log-normal + exponential
%   component), 'lnmix2' (timed log-normal + early log-normal).
%
%   Every quantity of the policy is a channel with an unconstrained value u_c(t):
%     timed  : z    log median of the timed component
%     early  : lnexp: log r (rate of the exponential), ell (logit of its weight)
%              lnmix2: z1 (log median of the early log-normal), ell
%   The timed width sigma is a constant of the session unless 'width' says otherwise.
%
%   Trial-to-trial dependence per component group ('timed', 'early'):
%     'fixed'         u_c = b_c
%     'history'       u_c(t) = b_c + w_c' h_t,  h_t = [log tau_{t-1}; r_{t-1}] (nLags lags)
%     'learn'         u_c(t) = b_c(t), with b_c updated after every trial by the learning rule
%     'historyLearn'  u_c(t) = W_c(t)' [1; h_t], with W_c updated by the rule (policy gradient) or
%                     its intercept updated by the rule (other rules)
%
%   Learning signal delta_t ('signal'):
%     'rate'      delta = r - rho,  rho = (1 - beta) rho + beta r     reward-rate baseline
%     'unsigned'  delta = |r - rho|                                   unsigned prediction error
%     'td'        delta = r - rho + gamma V - V,  V <- V + alpha_v delta   scalar TD critic
%     'q'         delta = r - Q(tau_t),  Q(tau) = v' [1, log tau, log tau^2],  v <- v + alpha_q delta phi
%     'raw'       delta = r
%
%   Learning rule for the timed location ('rule'); the early channels always use the policy
%   gradient with the same signal:
%     'pg'         z <- z + alpha delta g,  g = d log p / d z                       policy gradient
%     'threshold'  B = z - m is an estimate of the log delay:                      threshold tracking
%                  rewarded and log tau < B:   B <- B + alpha_down (log tau - B)
%                  unrewarded and log tau > B: B <- B + alpha_up (log tau - B);  z = B + m
%                  'tieRates' true: alpha_down = alpha_up;  'margin' true adds the free offset m
%     'stepq'      step-shaped action value Q(tau) = sigmoid((log tau - B) / kappa), kappa = 0.15,
%                  delta = r - Q(tau_t),  B <- B + alpha delta dQ/dB,  z = B + m
%     'bayes'      grid posterior over B from the censored outcomes (unrewarded at tau: B > log tau;
%                  rewarded: B <= log tau; lapse eps_b), forgetting lambda_b; z = posterior mean + m
%     'shift'      z <- z - alpha delta, alpha of either sign
%     'qmax'       z <- z + alpha (argmax Q - z), requires signal 'q'
%
%   The three featured models of the library are the policy-gradient actor ('pg' with 'rate'), the
%   unsigned-error actor ('pg' with 'unsigned') and the threshold-tracking actor ('threshold').
%
%   Further options: 'wRule' ('none' | 'pg' | 'unsigned') lets the history weights of the timed
%   channel learn; 'interaction' adds w_z_tauR log tau_{t-1} r_{t-1} to the timed location;
%   'width' ('fixed' | 'history' | 'state') lets sigma vary; 'state' adds a post-reward state
%   s_t = lambda_s s_{t-1} + r_{t-1} to the timed location (c_sz) and the early logit (c_spi).
%   Nothing learns during the first 'warmup' trials.
%
%   model.names, lb, ub   parameter names and bounds;  model.idx  name -> index;  model.chan  channels;
%   model.shortName       a file-name-safe identifier built from the options
p = inputParser;
p.addParameter('family', 'lnexp', @(x) any(strcmp(x, {'ln', 'lnexp', 'lnmix2'})));
p.addParameter('timed', 'historyLearn', @(x) any(strcmp(x, {'fixed', 'history', 'learn', 'historyLearn'})));
p.addParameter('early', 'fixed', @(x) any(strcmp(x, {'fixed', 'history', 'learn', 'historyLearn'})));
p.addParameter('signal', 'rate', @(x) any(strcmp(x, {'rate', 'unsigned', 'td', 'q', 'raw', 'abs'})));
p.addParameter('rule', 'pg', @(x) any(strcmp(x, {'pg', 'threshold', 'stepq', 'bayes', 'shift', 'qmax', 'boundary'})));
p.addParameter('tieRates', true, @islogical);
p.addParameter('margin', false, @islogical);
p.addParameter('wRule', 'none', @(x) any(strcmp(x, {'none', 'pg', 'unsigned', 'abs'})));
p.addParameter('interaction', false, @islogical);
p.addParameter('width', 'fixed', @(x) any(strcmp(x, {'fixed', 'history', 'state'})));
p.addParameter('state', false, @islogical);
p.addParameter('nLags', 1, @(x) ismember(x, [1 2]));
p.addParameter('warmup', 10, @isnumeric);
p.addParameter('beta', 0.02, @isnumeric);
p.addParameter('earlyShareAlpha', true, @islogical);
p.parse(varargin{:}); o = p.Results;
if strcmp(o.signal, 'abs'), o.signal = 'unsigned'; end
if strcmp(o.rule, 'boundary'), o.rule = 'threshold'; end
if strcmp(o.wRule, 'abs'), o.wRule = 'unsigned'; end
if strcmp(o.rule, 'qmax') && ~strcmp(o.signal, 'q'), error('ttModel:qmax', 'rule qmax requires signal q'); end
if strcmp(o.width, 'state') && ~o.state, error('ttModel:width', 'width state requires state true'); end
model = o; L = o.nLags; nH = 2 * L;
switch o.family
    case 'ln',     cn = {'z'};                 cg = {'timed'};
    case 'lnexp',  cn = {'z', 'logr', 'ell'};  cg = {'timed', 'early', 'early'};
    case 'lnmix2', cn = {'z', 'z1', 'ell'};    cg = {'timed', 'early', 'early'};
end
names = {}; lb = []; ub = []; chan = struct('name', {}, 'group', {}, 'mode', {}, 'iB', {}, 'iW', {}, 'iA', {}, 'rule', {}, 'iA2', {}, 'iM', {}, 'iX', {}, 'iAw', {}, 'iL', {}, 'iE', {});
sharedA = [];
for c = 1:numel(cn)
    mode = o.(cg{c}); ch = struct('name', cn{c}, 'group', cg{c}, 'mode', mode, 'iB', [], 'iW', [], 'iA', [], 'rule', 'pg', 'iA2', [], 'iM', [], 'iX', [], 'iAw', [], 'iL', [], 'iE', []);
    if strcmp(cg{c}, 'timed'), ch.rule = o.rule; end
    switch cn{c}
        case 'z',    [names, lb, ub] = addP(names, lb, ub, 'b_z', log(0.03), log(8));
        case 'z1',   [names, lb, ub] = addP(names, lb, ub, 'b_z1', log(0.02), log(8));
        case 'logr', [names, lb, ub] = addP(names, lb, ub, 'b_logr', log(0.02), log(50));
        case 'ell',  [names, lb, ub] = addP(names, lb, ub, 'b_ell', -7, 4);
    end
    ch.iB = numel(names);
    if any(strcmp(mode, {'history', 'historyLearn'}))
        for k = 1:L, [names, lb, ub] = addP(names, lb, ub, sprintf('w_%s_tau%d', cn{c}, k), -2, 2); end
        for k = 1:L, [names, lb, ub] = addP(names, lb, ub, sprintf('w_%s_R%d', cn{c}, k), -3, 3); end
        ch.iW = numel(names) - nH + 1:numel(names);
        if strcmp(cg{c}, 'timed') && o.interaction, [names, lb, ub] = addP(names, lb, ub, 'w_z_tauR', -3, 3); ch.iX = numel(names); end
        if strcmp(cg{c}, 'timed') && ~strcmp(o.wRule, 'none'), [names, lb, ub] = addP(names, lb, ub, 'alpha_w', 0, 1); ch.iAw = numel(names); end
    end
    if any(strcmp(mode, {'learn', 'historyLearn'}))
        if strcmp(cg{c}, 'early') && o.earlyShareAlpha && ~isempty(sharedA)
            ch.iA = sharedA;
        elseif strcmp(ch.rule, 'threshold')
            [names, lb, ub] = addP(names, lb, ub, 'alpha_up', 0, 1); ch.iA = numel(names);
            if o.tieRates, ch.iA2 = ch.iA; else, [names, lb, ub] = addP(names, lb, ub, 'alpha_down', 0, 1); ch.iA2 = numel(names); end
            if o.margin, [names, lb, ub] = addP(names, lb, ub, 'm_z', -1, 1.5); ch.iM = numel(names); end
        elseif strcmp(ch.rule, 'stepq')
            [names, lb, ub] = addP(names, lb, ub, 'alpha_z', 0, 1); ch.iA = numel(names);
            if o.margin, [names, lb, ub] = addP(names, lb, ub, 'm_z', -1, 1.5); ch.iM = numel(names); end
        elseif strcmp(ch.rule, 'bayes')
            [names, lb, ub] = addP(names, lb, ub, 'lambda_b', 0.5, 0.999); ch.iL = numel(names);
            [names, lb, ub] = addP(names, lb, ub, 'eps_b', 0.005, 0.4); ch.iE = numel(names);
            if o.margin, [names, lb, ub] = addP(names, lb, ub, 'm_z', -1, 1.5); ch.iM = numel(names); end
        else
            [names, lb, ub] = addP(names, lb, ub, sprintf('alpha_%s', ternary(strcmp(cg{c}, 'early') && o.earlyShareAlpha, 'early', cn{c})), 0, 1);
            ch.iA = numel(names); if strcmp(cg{c}, 'early') && o.earlyShareAlpha, sharedA = ch.iA; end
            if strcmp(ch.rule, 'shift'), lb(end) = -1; end
        end
    end
    chan(end+1) = ch; %#ok<AGROW>
end
switch o.family
    case {'ln', 'lnexp'}, [names, lb, ub] = addP(names, lb, ub, 'sigma', 0.02, 3);
    case 'lnmix2',        [names, lb, ub] = addP(names, lb, ub, 'sigma', 0.02, 3); [names, lb, ub] = addP(names, lb, ub, 's1', 0.02, 3);
end
if strcmp(o.width, 'history'), for k = 1:L, [names, lb, ub] = addP(names, lb, ub, sprintf('w_lsig_tau%d', k), -1, 1); end; for k = 1:L, [names, lb, ub] = addP(names, lb, ub, sprintf('w_lsig_R%d', k), -2, 2); end, end
switch o.signal
    case 'td', [names, lb, ub] = addP(names, lb, ub, 'alpha_v', 0, 1); [names, lb, ub] = addP(names, lb, ub, 'gamma', 0, 1);
    case 'q',  [names, lb, ub] = addP(names, lb, ub, 'alpha_q', 0, 1);
end
if o.state
    [names, lb, ub] = addP(names, lb, ub, 'lambda_s', 0, 0.98); [names, lb, ub] = addP(names, lb, ub, 'c_sz', -2, 2);
    if any(strcmp(cn, 'ell')), [names, lb, ub] = addP(names, lb, ub, 'c_spi', -3, 3); end
    if strcmp(o.width, 'state'), [names, lb, ub] = addP(names, lb, ub, 'c_ssig', -2, 2); end
end
model.names = names; model.lb = lb; model.ub = ub; model.nParam = numel(names);
model.idx = struct(); for i = 1:numel(names), model.idx.(names{i}) = i; end
model.chan = chan; model.nH = nH;
sn = sprintf('%s_T%s_E%s', o.family, o.timed, o.early); if strcmp(o.family, 'ln'), sn = sprintf('ln_T%s', o.timed); end
if ~strcmp(o.signal, 'rate'), sn = [sn, '_S', o.signal]; end
if ~strcmp(o.rule, 'pg'), sn = [sn, '_R', o.rule]; end
if strcmp(o.rule, 'threshold') && ~o.tieRates, sn = [sn, '_tworates']; end
if o.margin, sn = [sn, '_margin']; end
if ~strcmp(o.wRule, 'none'), sn = [sn, '_W', o.wRule]; end
if o.interaction, sn = [sn, '_int']; end
if ~strcmp(o.width, 'fixed'), sn = [sn, '_sig', o.width]; end
if o.state, sn = [sn, '_state']; end
if L == 2, sn = [sn, '_lag2']; end
model.shortName = sn;
model.description = sprintf('%s: timed %s, early %s, signal %s, rule %s, %d lag(s), %d parameters', o.family, o.timed, o.early, o.signal, o.rule, L, model.nParam);
end
function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
function [names, lb, ub] = addP(names, lb, ub, nm, l, u), names = [names, {nm}]; lb = [lb, l]; ub = [ub, u]; end
