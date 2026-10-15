function stage = ttClassifyStage(delays, nTrials, varargin)
%TTCLASSIFYSTAGE  Protocol stage of a session from its per-trial delays.
%
%   stage = ttClassifyStage(delays, nTrials)
%   stage = ttClassifyStage(delays, nTrials, 'cueMax', 0.15, 'cueTrials', 360, 'switchRatio', 2, 'switchMinFrac', 0.15)
%
%   delays   delay of every cued trial, s
%   nTrials  number of trials in the session
%   Stages (string):
%     "cueAssociation"  one delay level at or below cueMax and fewer than cueTrials trials
%     "fixedDelay"      one delay level above cueMax (or at or below cueMax with cueTrials or more trials)
%     "switching"       two or three delay levels, each holding at least switchMinFrac of the trials,
%                       with the largest at least switchRatio times the smallest
%     "delayTraining"   any other set of delay levels (a ramp of neighbouring delays)
%     "unclassified"    no delay information
%   Delay levels are the distinct values after rounding to 0.01 s; levels at or below 0.05 s are
%   ignored.
p = inputParser; p.addParameter('cueMax', 0.15, @isnumeric); p.addParameter('cueTrials', 360, @isnumeric);
p.addParameter('switchRatio', 2, @isnumeric); p.addParameter('switchMinFrac', 0.15, @isnumeric); p.parse(varargin{:}); o = p.Results;
stage = "unclassified";
d = round(delays(isfinite(delays)), 2); d = d(d > 0.05); if isempty(d), return; end
u = unique(d);
if numel(u) == 1
    if u <= o.cueMax && nTrials < o.cueTrials, stage = "cueAssociation"; else, stage = "fixedDelay"; end
    return
end
frac = arrayfun(@(x) mean(d == x), u);
if numel(u) <= 3 && all(frac >= o.switchMinFrac) && max(u) / min(u) >= o.switchRatio, stage = "switching"; return; end
stage = "delayTraining";
end
