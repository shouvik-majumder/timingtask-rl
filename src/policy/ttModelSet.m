function V = ttModelSet(which)
%TTMODELSET  Named sets of model specifications (cell arrays of ttModel argument lists).
%
%   V = ttModelSet('core')        the reference models: static, history, and the three actors with history
%   V = ttModelSet('families')    the core set for each policy family (ln, lnexp, lnmix2)
%   V = ttModelSet('early')       the three actors with the early component also moving
%   V = ttModelSet('signals')     learning signal x rule grid on the timed location (no history)
%   V = ttModelSet('extended')    core plus the post-reward state, margin, two-rate and two-lag variants
%
%   Core set (family lnexp, early component fixed):
%     static            timed fixed
%     history           timed location depends on the previous lick and reward
%     policyGradient    history + policy-gradient actor, reward-rate baseline
%     unsignedError     history + policy-gradient actor driven by the unsigned prediction error
%     threshold         history + threshold-tracking actor (one rate, no margin)
if nargin < 1, which = 'core'; end
core = @(fam) { ...
    {'family', fam, 'timed', 'fixed', 'early', 'fixed'}, ...
    {'family', fam, 'timed', 'history', 'early', 'fixed'}, ...
    {'family', fam, 'timed', 'historyLearn', 'early', 'fixed', 'signal', 'rate', 'rule', 'pg'}, ...
    {'family', fam, 'timed', 'historyLearn', 'early', 'fixed', 'signal', 'unsigned', 'rule', 'pg'}, ...
    {'family', fam, 'timed', 'historyLearn', 'early', 'fixed', 'rule', 'threshold'}};
switch which
    case 'core'
        V = core('lnexp');
    case 'families'
        V = [core('lnexp'), core('lnmix2'), {{'family', 'ln', 'timed', 'fixed'}, {'family', 'ln', 'timed', 'history'}, {'family', 'ln', 'timed', 'historyLearn', 'rule', 'pg'}, {'family', 'ln', 'timed', 'historyLearn', 'signal', 'unsigned'}, {'family', 'ln', 'timed', 'historyLearn', 'rule', 'threshold'}}];
    case 'early'
        V = {{'family', 'lnexp', 'timed', 'historyLearn', 'early', 'historyLearn', 'rule', 'pg'}, ...
             {'family', 'lnexp', 'timed', 'historyLearn', 'early', 'historyLearn', 'signal', 'unsigned'}, ...
             {'family', 'lnexp', 'timed', 'historyLearn', 'early', 'historyLearn', 'rule', 'threshold'}};
    case 'signals'
        V = {};
        for sg = {'rate', 'unsigned', 'td', 'q', 'raw'}, V{end+1} = {'family', 'lnexp', 'timed', 'learn', 'early', 'fixed', 'signal', sg{1}, 'rule', 'pg'}; end %#ok<AGROW>
        V{end+1} = {'family', 'lnexp', 'timed', 'learn', 'early', 'fixed', 'rule', 'threshold'};
        V{end+1} = {'family', 'lnexp', 'timed', 'learn', 'early', 'fixed', 'rule', 'stepq'};
        V{end+1} = {'family', 'lnexp', 'timed', 'learn', 'early', 'fixed', 'rule', 'bayes'};
        V{end+1} = {'family', 'lnexp', 'timed', 'learn', 'early', 'fixed', 'signal', 'q', 'rule', 'qmax'};
    case 'extended'
        V = core('lnexp');
        V{end+1} = {'family', 'lnexp', 'timed', 'historyLearn', 'early', 'fixed', 'rule', 'threshold', 'margin', true};
        V{end+1} = {'family', 'lnexp', 'timed', 'historyLearn', 'early', 'fixed', 'rule', 'threshold', 'tieRates', false};
        V{end+1} = {'family', 'lnexp', 'timed', 'historyLearn', 'early', 'fixed', 'rule', 'pg', 'state', true};
        V{end+1} = {'family', 'lnexp', 'timed', 'historyLearn', 'early', 'fixed', 'rule', 'threshold', 'state', true};
        V{end+1} = {'family', 'lnexp', 'timed', 'historyLearn', 'early', 'fixed', 'rule', 'pg', 'nLags', 2};
        V{end+1} = {'family', 'lnexp', 'timed', 'historyLearn', 'early', 'fixed', 'rule', 'threshold', 'nLags', 2};
    otherwise
        error('ttModelSet:unknown', 'unknown set %s', which);
end
end
