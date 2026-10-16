function fig = ttPlotDistParams(CP, S, levels, varargin)
%TTPLOTDISTPARAMS  Fitted quantities of the log-normal + exponential distribution by level and series (dot and whiskers over sessions).
%
%   fig = ttPlotDistParams(CP, S, levels)
%   fig = ttPlotDistParams(CP, S, levels, 'series', 'group', 'level', 'phase', 'unitType', 'session', 'x', 'delayLevel')
%
%   CP is dist_components.csv (ttDistSummary); S is D.sessions. Four panels: timed median (s), timed
%   width (log s), early-component weight, early-component rate (per s). With 'x' set to
%   'delayLevel' the x axis is the delay of the block (unitType delayBlock) and the levels are
%   ignored.
p = inputParser; p.addParameter('series', 'group', @ischar); p.addParameter('level', 'phase', @ischar); p.addParameter('unitType', 'session', @ischar); p.addParameter('x', '', @ischar); p.addParameter('colors', [], @isnumeric); p.addParameter('tickAngle', 30, @isnumeric); p.parse(varargin{:}); o = p.Results;
E = CP(CP.family == "ln_exp" & CP.unitType == string(o.unitType), :); [tf, loc] = ismember(E.session, S.session); E = E(tf, :);
E.(o.level) = S.(o.level)(loc(tf)); if ~strcmp(o.series, 'none'), E.(o.series) = S.(o.series)(loc(tf)); end
vars = {'median_s', 'width_log', 'weight_exp', 'rate_per_s'}; vl = {'timed median (s)', 'timed width (log s)', 'early-component weight', 'early-component rate (per s)'};
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1600 380]);
for v = 1:4
    subplot(1, 4, v); ser = o.series; if strcmp(ser, 'none'), ser = ''; end
    if isempty(o.x)
        ttDotWhisker(E, levels, 'series', ser, 'value', vars{v}, 'level', o.level, 'colors', o.colors, 'tickAngle', o.tickAngle);
    else
        E.xl = string(round(E.(o.x), 1)); lv = unique(E.(o.x)); lv = lv(arrayfun(@(d) sum(E.(o.x) == d) >= 5, lv)); ttDotWhisker(E, string(round(lv, 1)), 'series', ser, 'value', vars{v}, 'level', 'xl', 'colors', o.colors, 'showN', false); xlabel('delay of the block (s)');
        if v == 1, hold on; plot(1:numel(lv), lv, 'k:'); end
    end
    if v == 4, set(gca, 'YScale', 'log'); end; title(vl{v}, 'FontWeight', 'normal');
end
end
