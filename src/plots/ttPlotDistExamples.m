function fig = ttPlotDistExamples(D, animal, varargin)
%TTPLOTDISTEXAMPLES  First-lick-time histograms of one session per phase with fitted densities.
%
%   fig = ttPlotDistExamples(D, 'RD50')
%   fig = ttPlotDistExamples(D, animal, 'families', {'lognormal', 'ln_exp', 'lnmix2'}, 'phases', {...})
%
%   For each phase the session with the most cued trials is shown. Densities are whole-session
%   maximum-likelihood fits (ttDistFamilies). Delay levels are marked.
p = inputParser; p.addParameter('families', {'lognormal', 'ln_exp', 'lnmix2'}, @iscell); p.addParameter('phases', {}, @iscell); p.parse(varargin{:}); o = p.Results;
S = D.sessions(D.sessions.animal == string(animal), :); T = D.trials(D.trials.cued & D.trials.animal == string(animal), :);
phases = string(o.phases); if isempty(phases), phases = unique(S.phase, 'stable')'; end
fam = ttDistFamilies(o.families); cols = lines(numel(fam)); nP = numel(phases); nc = min(nP, 4); nr = ceil(nP / nc);
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [50 50 350 * nc, 300 * nr]);
for j = 1:nP
    c = S(S.phase == phases(j), :); if isempty(c), continue; end
    [~, i] = max(c.nCued); sid = c.session(i); x = T.lickTime(T.session == sid); dl = T.delay(T.session == sid); Tmax = max(x) * 1.001;
    subplot(nr, nc, j); hold on
    edges = 0:0.05:min(max(4, prctile(x, 99)), 8); histogram(x, edges, 'Normalization', 'pdf', 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
    z = linspace(0.005, edges(end), 600);
    for k = 1:numel(fam), try, th = fam(k).fit(x, Tmax); plot(z, exp(fam(k).logpdf(z, th, Tmax)), '-', 'Color', cols(k, :), 'LineWidth', 1.2); catch, end, end
    for d = unique(round(dl * 20) / 20)', xline(d, ':', 'Color', [0.3 0.3 0.3]); end
    xlabel('first lick time (s)'); ylabel('density'); title(sprintf('%s: session %d, %d trials', phases(j), c.sessionIndex(i), numel(x)), 'FontWeight', 'normal');
    if j == 1, legend([{'data'}, {fam.label}], 'Box', 'off', 'FontSize', 8); end
end
end
