function fig = ttPlotHistoryRegression(R, S, levels, varargin)
%TTPLOTHISTORYREGRESSION  Coefficients of the trial-history regression (ttHistoryRegression) by level and series, with the shuffle-null band.
%
%   fig = ttPlotHistoryRegression(R, S, levels)
%   fig = ttPlotHistoryRegression(R, S, levels, 'series', 'group', 'level', 'phase', 'subsets', {...}, 'colors', C)
%
%   R is the output of ttHistoryRegression; S is D.sessions (joined on session to supply the level
%   and series columns). One row of panels per subset, one column per coefficient (previous lick
%   time, previous reward, interaction).
p = inputParser; p.addParameter('series', 'group', @ischar); p.addParameter('level', 'phase', @ischar); p.addParameter('subsets', {'all', 'timedNext', 'earlyNext'}, @iscell); p.addParameter('colors', [], @isnumeric); p.parse(varargin{:}); o = p.Results;
[tf, loc] = ismember(R.session, S.session); R = R(tf, :); R.(o.level) = S.(o.level)(loc(tf)); if ~strcmp(o.series, 'none'), R.(o.series) = S.(o.series)(loc(tf)); end
coefs = {'prevLick', 'prevReward', 'interaction'}; coefL = {'previous log lick time (SD per SD)', 'previous reward (SD per reward)', 'previous log lick time x previous reward'}; subL = struct('all', 'all trials', 'timedNext', 'next lick from the timed component', 'earlyNext', 'next lick from the early component', 'timedPrev', 'previous lick from the timed component', 'earlyPrev', 'previous lick from the early component');
nS = numel(o.subsets); fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 320 * nS]);
for i = 1:nS
    for c = 1:3
        subplot(nS, 3, (i - 1) * 3 + c); q = R(R.subset == string(o.subsets{i}), :); if isempty(q), axis off; continue; end
        ser = o.series; if strcmp(ser, 'none'), ser = ''; end
        ttDotWhisker(q, levels, 'series', ser, 'value', ['b_', coefs{c}], 'level', o.level, 'nullLo', ['null_lo_', coefs{c}], 'nullHi', ['null_hi_', coefs{c}], 'colors', o.colors);
        title(sprintf('%s: %s', subL.(o.subsets{i}), coefL{c}), 'FontWeight', 'normal', 'FontSize', 9); if c == 1, ylabel('coefficient (median over sessions, IQR)'); end
    end
end
end
