function ttFitAll(D, variants, outDir, varargin)
%TTFITALL  Fit a set of models to every animal of a dataset; fits already on disk are skipped.
%
%   ttFitAll(D, ttModelSet('core'), outDir)
%   ttFitAll(D, variants, outDir, 'animals', {...}, 'ebIter', 8, 'nStarts', 3, 'redo', false)
%
%   variants is a cell array of ttModel argument lists. Fits go to <outDir>/fits, a log of timings
%   and failures to <outDir>/fits/ttFitAll.log.
p = inputParser; p.addParameter('animals', {}, @(x) iscell(x) || isstring(x)); p.addParameter('ebIter', 8, @isnumeric); p.addParameter('nStarts', 3, @isnumeric); p.addParameter('redo', false, @islogical); p.parse(varargin{:}); o = p.Results;
fd = fullfile(outDir, 'fits'); if ~exist(fd, 'dir'), mkdir(fd); end
animals = string(o.animals(:))'; if isempty(animals), animals = unique(D.sessions.animal, 'stable')'; end
logf = fullfile(fd, 'ttFitAll.log'); t0 = tic;
for v = 1:numel(variants)
    model = ttModel(variants{v}{:});
    for a = animals
        fn = fullfile(fd, sprintf('%s_%s.mat', a, model.shortName));
        if isfile(fn) && ~o.redo, continue; end
        t1 = tic;
        try
            ttFitAnimal(D, a, model, 'ebIter', o.ebIter, 'nStarts', o.nStarts, 'outDir', fd, 'verbose', false);
            msg = sprintf('%s  %-44s %-8s %.1f min (elapsed %.1f h)', datestr(now, 'HH:MM'), model.shortName, a, toc(t1) / 60, toc(t0) / 3600);
        catch ME
            msg = sprintf('%s  %-44s %-8s FAILED: %s', datestr(now, 'HH:MM'), model.shortName, a, ME.message);
        end
        fprintf('%s\n', msg); fid = fopen(logf, 'a'); fprintf(fid, '%s\n', msg); fclose(fid);
    end
end
fprintf('ttFitAll done in %.1f h\n', toc(t0) / 3600);
end
