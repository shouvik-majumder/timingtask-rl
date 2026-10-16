function ax = ttDotWhisker(X, levels, varargin)
%TTDOTWHISKER  Median with interquartile whiskers of a per-session measure, by level and series, with an optional null band.
%
%   ax = ttDotWhisker(X, levels)
%   ax = ttDotWhisker(X, levels, 'series', 'region', 'value', 'r', 'level', 'phase', 'nullLo', 'null_lo', 'nullHi', 'null_hi', 'colors', C, 'labels', {...})
%
%   X is a table. For every level of X.(level) (in the order of 'levels') and every value of
%   X.(series), one dot at the median of X.(value) with whiskers at the 25th and 75th percentiles, and
%   the number of sessions printed above. When nullLo / nullHi are given, a grey band spans the
%   median over sessions of the per-session null bounds, drawn per series next to its dot.
p = inputParser; p.addParameter('series', '', @ischar); p.addParameter('value', 'r', @ischar); p.addParameter('level', 'phase', @ischar); p.addParameter('nullLo', '', @ischar); p.addParameter('nullHi', '', @ischar);
p.addParameter('colors', [], @isnumeric); p.addParameter('labels', {}, @iscell); p.addParameter('minN', 3, @isnumeric); p.addParameter('showN', true, @islogical); p.addParameter('tickAngle', 30, @isnumeric); p.parse(varargin{:}); o = p.Results;
ax = gca; hold(ax, 'on'); levels = string(levels);
if isempty(o.series), ser = ""; X.one = repmat("", height(X), 1); sc = 'one'; else, ser = unique(string(X.(o.series)), 'stable')'; sc = o.series; end
nS = numel(ser); cols = o.colors; if isempty(cols), cols = [0.85 0.2 0.1; 0.2 0.2 0.2; 0.1 0.4 0.8; 0.2 0.6 0.2]; end; hh = gobjects(nS, 1);
for j = 1:numel(levels)
    for i = 1:nS
        m = string(X.(o.level)) == levels(j) & string(X.(sc)) == ser(i); v = X.(o.value)(m); v = v(isfinite(v)); if numel(v) < o.minN, continue; end
        xo = j + (i - (nS + 1) / 2) * 0.22; med = median(v); q1 = prctile(v, 25); q3 = prctile(v, 75);
        if ~isempty(o.nullLo), lo = median(X.(o.nullLo)(m), 'omitnan'); hi = median(X.(o.nullHi)(m), 'omitnan'); fill(xo + [-0.09 0.09 0.09 -0.09], [lo lo hi hi], [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.7); end
        hh(i) = errorbar(xo, med, med - q1, q3 - med, 'o', 'Color', cols(i, :), 'MarkerFaceColor', cols(i, :), 'MarkerSize', 5, 'LineWidth', 1.3, 'CapSize', 0);
        if o.showN, text(xo, q3, sprintf(' %d', numel(v)), 'FontSize', 6, 'Color', cols(i, :), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); end
    end
end
h0 = yline(0, 'k:'); h0.Annotation.LegendInformation.IconDisplayStyle = 'off'; set(ax, 'XTick', 1:numel(levels), 'XTickLabel', cellstr(levels)); xtickangle(o.tickAngle); xlim([0.4, numel(levels) + 0.6]);
keep = arrayfun(@isgraphics, hh); lab = cellstr(ser); if ~isempty(o.labels), lab = o.labels; end; if nS > 1 && any(keep), legend(hh(keep), lab(keep), 'Box', 'off', 'FontSize', 8, 'Location', 'best'); end
end
