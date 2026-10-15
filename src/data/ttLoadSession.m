function T = ttLoadSession(S, varargin)
%TTLOADSESSION  Trial table of one session from the laboratory's processed session structure.
%
%   T = ttLoadSession(S)
%   T = ttLoadSession(S, 'animal', 'RD50', 'session', 'RD50_250221_155134')
%
%   S is the structure saved by the laboratory's session processing (one .mat file per session,
%   variable S). The fields used:
%     S.Trial.firstLickTime        first lick after the cue, s (NaN when none)
%     S.Trial.ResponseTypeVector   0 early lick, 1 reward, 2 no response, 3 no-cue trial
%     S.TrialSetting.DelayPeriod   delay of each trial, s
%     S.numTrial                   number of trials
%   Output columns: animal, session, trialNum, lickTime, outcome, rewarded, delay, cued, noCue.
%   cued = outcome in {0, 1} with a first lick after the cue and a defined delay; these are the
%   trials every model is fitted to. noCue = outcome 3; lickTime is kept for those trials when a
%   lick occurred.
p = inputParser; p.addParameter('animal', '', @(x) ischar(x) || isstring(x)); p.addParameter('session', '', @(x) ischar(x) || isstring(x)); p.parse(varargin{:}); o = p.Results;
fl = S.Trial.firstLickTime(:); oc = S.Trial.ResponseTypeVector(:); dl = S.TrialSetting.DelayPeriod(:);
n = min([numel(fl), numel(oc), numel(dl)]); fl = fl(1:n); oc = double(oc(1:n)); dl = dl(1:n);
T = table(repmat(string(o.animal), n, 1), repmat(string(o.session), n, 1), (1:n)', fl, oc, oc == 1, dl, ...
    'VariableNames', {'animal', 'session', 'trialNum', 'lickTime', 'outcome', 'rewarded', 'delay'});
T.cued = ismember(oc, [0 1]) & ~isnan(fl) & fl > 0 & ~isnan(dl);
T.noCue = oc == 3;
end
