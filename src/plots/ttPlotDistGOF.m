function fig = ttPlotDistGOF(G, S, levels, varargin)
%TTPLOTDISTGOF  Split-half goodness of fit (ttDistGOF) by level and series: model-to-data KS distance against the data-to-data reference, and quantile R2.
%
%   fig = ttPlotDistGOF(G, S, levels, 'series', 'group', 'level', 'phase')
p = inputParser; p.addParameter('series', 'group', @ischar); p.addParameter('level', 'phase', @ischar); p.addParameter('colors', [], @isnumeric); p.parse(varargin{:}); o = p.Results;
[tf, loc] = ismember(G.session, S.session); G = G(tf, :); G.(o.level) = S.(o.level)(loc(tf)); if ~strcmp(o.series, 'none'), G.(o.series) = S.(o.series)(loc(tf)); end
G.ks_data_lo = zeros(height(G), 1); ser = o.series; if strcmp(ser, 'none'), ser = ''; end
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 380]);
subplot(1, 3, 1); ttDotWhisker(G, levels, 'series', ser, 'value', 'ks_model', 'level', o.level, 'nullLo', 'ks_data_lo', 'nullHi', 'ks_data_p95', 'colors', o.colors); ylabel('KS distance'); title('fitted CDF (half A) vs empirical CDF (half B); grey: up to the 95th percentile of half A vs half B', 'FontWeight', 'normal', 'FontSize', 8); ylim([0 0.3]);
subplot(1, 3, 2); ttDotWhisker(G, levels, 'series', ser, 'value', 'frac_model_worse', 'level', o.level, 'colors', o.colors); h = yline(0.5, 'k--'); h.Annotation.LegendInformation.IconDisplayStyle = 'off'; ylabel('fraction of splits'); title('fraction of splits where the model is farther from half B than half A is', 'FontWeight', 'normal', 'FontSize', 8); ylim([0 1]);
subplot(1, 3, 3); ttDotWhisker(G, levels, 'series', ser, 'value', 'r2_quant', 'level', o.level, 'colors', o.colors); ylabel('R2'); title('R2 of fitted quantiles (half A) against empirical quantiles (half B), log s', 'FontWeight', 'normal', 'FontSize', 8); ylim([0.5 1]);
end
