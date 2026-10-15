function [Bs, S] = ttSessionBlocks(D, animal, varargin)
%TTSESSIONBLOCKS  Per-session fitting structures of one animal from a dataset.
%
%   [Bs, S] = ttSessionBlocks(D, animal)
%   [Bs, S] = ttSessionBlocks(D, animal, 'sessions', {...}, 'minTrials', 50)
%
%   Bs{s} holds the cued trials of one session in trial order: tau (lick time, s), R (reward 0/1),
%   delay (s), trialNum, session, sessionIndex, phase. S is the matching subset of D.sessions.
p = inputParser; p.addParameter('sessions', {}, @(x) iscell(x) || isstring(x)); p.addParameter('minTrials', 50, @isnumeric); p.parse(varargin{:}); o = p.Results;
S = sortrows(D.sessions(D.sessions.animal == string(animal), :), 'sessionIndex');
if ~isempty(o.sessions), S = S(ismember(S.session, string(o.sessions)), :); end
T = D.trials(D.trials.cued & D.trials.animal == string(animal), :);
Bs = {}; keep = false(height(S), 1);
for s = 1:height(S)
    q = sortrows(T(T.session == S.session(s), :), 'trialNum');
    if height(q) < o.minTrials, continue; end
    Bs{end+1} = struct('tau', q.lickTime, 'R', double(q.rewarded), 'delay', q.delay, 'trialNum', q.trialNum, 'session', char(S.session(s)), 'sessionIndex', S.sessionIndex(s), 'phase', char(S.phase(s))); %#ok<AGROW>
    keep(s) = true;
end
S = S(keep, :);
end
