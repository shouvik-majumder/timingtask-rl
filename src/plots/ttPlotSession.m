function fig = ttPlotSession(Res, session, varargin)
%TTPLOTSESSION  Lick times of one session with the fitted model's timed median and latent variables.
%
%   fig = ttPlotSession(Res, session)
%   fig = ttPlotSession(Res, session, 'compare', Res2)
%
%   Res is the output of ttFitAnimal (or the Res variable of a saved fit); session is a session
%   identifier in Res.trials. Panels: lick times with the delay and the model's timed median (and
%   that of 'compare' when given); probability that the lick belongs to the early component; the
%   learning signal and the step applied to the timed location; the timed location in log seconds.
p = inputParser; p.addParameter('compare', [], @(x) isempty(x) || isstruct(x)); p.parse(varargin{:}); o = p.Results;
t = Res.trials(Res.trials.session == string(session), :); n = height(t); x = 1:n;
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 1400 700]);
subplot(4, 1, 1:2); hold on
scatter(x, t.tau, 8, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.5); plot(x, t.delay, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
lg = {'lick time', 'delay'};
if ~isempty(o.compare), t2 = o.compare.trials(o.compare.trials.session == string(session), :); plot(x, t2.timedMedian, 'k-', 'LineWidth', 1); lg{end+1} = ['timed median, ', strrep(o.compare.model.shortName, '_', ' ')]; end
plot(x, t.timedMedian, '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.3); lg{end+1} = ['timed median, ', strrep(Res.model.shortName, '_', ' ')];
ylim([0 prctile(t.tau, 98) * 1.2]); xlim([0 n + 1]); ylabel('lick time (s)'); legend(lg, 'Box', 'off', 'Location', 'northwest', 'FontSize', 8, 'NumColumns', 4);
title(sprintf('%s, session %s (%d trials)', Res.animal, session, n), 'FontWeight', 'normal', 'Interpreter', 'none');
subplot(4, 1, 3); hold on; stem(x, t.resp, 'Marker', 'none', 'Color', [0.4 0.4 0.4]); ylim([0 1]); xlim([0 n + 1]); ylabel('P(early component)');
subplot(4, 1, 4); hold on; plot(x, t.delta, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
if ismember('step_z', t.Properties.VariableNames), stem(x, t.step_z, 'Marker', 'none', 'Color', [0.2 0.6 0.2]); end
plot(x, t.u_z - t.u_z(1), 'r-', 'LineWidth', 1); xlim([0 n + 1]); ylabel({'learning signal (grey)', 'step (green), location change (red)'}); xlabel('trial');
end
