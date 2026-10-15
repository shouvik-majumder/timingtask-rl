function D = ttAssignPhases(D, varargin)
%TTASSIGNPHASES  Training phases within each animal from the protocol stages of its sessions.
%
%   D = ttAssignPhases(D)
%   D = ttAssignPhases(D, 'delaySplit', 2, 'fixedSplit', 1)
%
%   Phase names, per animal and in session order:
%     "cueAssoc"                    cue-association sessions
%     "delayEarly" / "delayLate"    delay-training sessions split into delaySplit consecutive parts
%                                   ("delay1".."delayK" for other K)
%     "fixedPre" / "fixedPost"      fixed-delay sessions before / after the first switching session
%                                   (just "fixed" when the animal has no switching sessions); each
%                                   further split into fixedSplit parts ("fixedPreEarly" ...)
%     "switching"                   switching-delay sessions
%     "other"                       unclassified sessions
%   D.sessions gains phase, stageIndex (position within the animal's sessions of that stage, from 1)
%   and stageIndexFromEnd (-1 = last session of the stage); D.trials gains phase.
p = inputParser; p.addParameter('delaySplit', 2, @isnumeric); p.addParameter('fixedSplit', 1, @isnumeric); p.parse(varargin{:}); o = p.Results;
S = D.sessions; ph = repmat("other", height(S), 1); si = nan(height(S), 1); se = nan(height(S), 1);
for a = unique(S.animal)'
    ia = find(S.animal == a); [~, ord] = sort(S.sessionIndex(ia)); ia = ia(ord); st = S.stage(ia);
    for stage = ["cueAssociation", "delayTraining", "fixedDelay", "switching"]
        k = find(st == stage); n = numel(k); if n == 0, continue; end
        si(ia(k)) = 1:n; se(ia(k)) = (1:n) - n - 1;
    end
    ph(ia(st == "cueAssociation")) = "cueAssoc";
    ph(ia(st == "switching")) = "switching";
    k = find(st == "delayTraining"); ph(ia(k)) = splitNames("delay", numel(k), o.delaySplit);
    k = find(st == "fixedDelay"); firstSw = find(st == "switching", 1);
    if isempty(firstSw)
        ph(ia(k)) = splitNames("fixed", numel(k), o.fixedSplit);
    else
        pre = k(k < firstSw); post = k(k > firstSw);
        ph(ia(pre)) = splitNames("fixedPre", numel(pre), o.fixedSplit); ph(ia(post)) = splitNames("fixedPost", numel(post), o.fixedSplit);
    end
end
S.phase = ph; S.stageIndex = si; S.stageIndexFromEnd = se; D.sessions = S;
[tf, loc] = ismember(D.trials.animal + "|" + D.trials.session, S.animal + "|" + S.session);
D.trials.phase = repmat("other", height(D.trials), 1); D.trials.phase(tf) = S.phase(loc(tf));
end
function nm = splitNames(base, n, k)
nm = strings(n, 1); if n == 0, return; end
k = max(1, min(k, n)); part = ceil((1:n)' * k / n);
if k == 1, nm(:) = base; return; end
if k == 2, lab = ["Early", "Late"]; elseif k == 3, lab = ["Early", "Mid", "Late"]; else, lab = string(1:k); end
for i = 1:n, nm(i) = base + lab(part(i)); end
end
