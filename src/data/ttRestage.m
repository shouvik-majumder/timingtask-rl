function D = ttRestage(D, varargin)
%TTRESTAGE  Recompute the protocol stage of every session from the trials already in a dataset, then reassign phases.
%
%   D = ttRestage(D)
%   D = ttRestage(D, 'stageOptions', {...}, 'phaseOptions', {...})
%
%   Avoids reloading the session files when the classification rules change.
p = inputParser; p.addParameter('stageOptions', {}, @iscell); p.addParameter('phaseOptions', {}, @iscell); p.parse(varargin{:}); o = p.Results;
S = D.sessions; T = D.trials;
for s = 1:height(S)
    m = T.session == S.session(s); dc = T.delay(m & T.cued); [S.stage(s), lv] = ttClassifyStage(dc, sum(m), o.stageOptions{:}); S.delayLevels(s) = string(mat2str(lv, 3));
end
D.sessions = S; [tf, loc] = ismember(T.animal + "|" + T.session, S.animal + "|" + S.session); D.trials.stage(tf) = S.stage(loc(tf));
D = ttAssignPhases(D, o.phaseOptions{:});
end
