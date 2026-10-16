function [M, fig] = ttGroupTrajectories(P, measure, varargin)
%TTGROUPTRAJECTORIES  Per-group trajectory of a per-session measure against an alignment index.
%
%   M = ttGroupTrajectories(P, 'alpha_up')
%   [M, fig] = ttGroupTrajectories(P, measure, 'x', 'stageIndex', 'stage', 'delayTraining', 'groups', {...}, 'transform', @log10, 'plot', true)
%
%   P is a table with columns animal, group, stage, stageIndex, stageIndexFromEnd, sessionIndex,
%   delayMean and the measure. For each group and value of the alignment variable 'x' (restricted
%   to 'stage' when given): the number of animals, the median over animals of the per-animal
%   median, and the 25th and 75th percentiles over animals. With 'plot' true a figure is drawn,
%   one line per group with the interquartile band.
p = inputParser; p.addParameter('x', 'stageIndex', @ischar); p.addParameter('stage', '', @ischar); p.addParameter('groups', {}, @iscell); p.addParameter('transform', @(x) x, @(x) isa(x, 'function_handle'));
p.addParameter('plot', false, @islogical); p.addParameter('colors', [], @isnumeric); p.addParameter('label', '', @ischar); p.parse(varargin{:}); o = p.Results;
P = P(isfinite(P.(measure)), :); if ~isempty(o.stage), P = P(P.stage == string(o.stage), :); end
P.value = o.transform(P.(measure)); groups = string(o.groups); if isempty(groups), groups = unique(P.group, 'stable')'; end
rows = {};
for g = groups
    q = P(P.group == g, :);
    for x = unique(q.(o.x))'
        qq = q(q.(o.x) == x, :); an = unique(qq.animal); v = arrayfun(@(a) median(qq.value(qq.animal == a)), an);
        rows{end+1} = table(g, x, numel(an), median(v), prctile(v, 25), prctile(v, 75), 'VariableNames', {'group', 'x', 'nAnimals', 'median', 'q25', 'q75'}); %#ok<AGROW>
    end
end
if isempty(rows), M = table('Size', [0 6], 'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'group', 'x', 'nAnimals', 'median', 'q25', 'q75'}); else, M = vertcat(rows{:}); end; fig = [];
if o.plot
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 700 420]); hold on; cols = o.colors; if isempty(cols), cols = lines(numel(groups)); end; hh = gobjects(numel(groups), 1);
    for i = 1:numel(groups)
        q = sortrows(M(M.group == groups(i), :), 'x'); if isempty(q), continue; end
        fill([q.x; flipud(q.x)], [q.q25; flipud(q.q75)], cols(i, :), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        hh(i) = plot(q.x, q.median, 'o-', 'Color', cols(i, :), 'MarkerFaceColor', cols(i, :), 'MarkerSize', 4, 'LineWidth', 1.3);
    end
    keep = arrayfun(@isgraphics, hh); legend(hh(keep), cellstr(groups(keep)), 'Box', 'off', 'Location', 'best');
    xlabel(strrep(o.x, 'stageIndexFromEnd', 'session relative to the end of the stage')); lab = o.label; if isempty(lab), lab = strrep(measure, '_', '\_'); end; ylabel(lab);
    if ~isempty(o.stage), title(sprintf('%s sessions: median over animals with the interquartile range', o.stage), 'FontWeight', 'normal'); end
end
end
