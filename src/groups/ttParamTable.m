function P = ttParamTable(outDir, variant)
%TTPARAMTABLE  Fitted parameters of one model, one row per session, with the session's stage, phase, group and alignment indices.
%
%   P = ttParamTable(outDir, 'lnexp_ThistoryLearn_Efixed_Rthreshold')
%
%   Reads <outDir>/fits/<animal>_<variant>.mat. Columns: animal, session, sessionIndex, stage,
%   phase, group, stageIndex, stageIndexFromEnd, delayMean, nTrials, one column per parameter,
%   logLik, llPerTrial, BIC.
fd = fullfile(outDir, 'fits'); files = dir(fullfile(fd, ['*_', variant, '.mat'])); rows = {};
for f = 1:numel(files)
    L = load(fullfile(fd, files(f).name), 'Res'); if ~strcmp(L.Res.model.shortName, variant), continue; end
    q = L.Res.params; rows{end+1} = q(:, [{'animal', 'session', 'sessionIndex', 'stage', 'phase', 'group', 'stageIndex', 'stageIndexFromEnd', 'delayMean', 'nTrials'}, L.Res.model.names, {'logLik', 'llPerTrial', 'BIC'}]); %#ok<AGROW>
end
if isempty(rows), error('ttParamTable:none', 'no fits of %s in %s', variant, fd); end
P = vertcat(rows{:});
end
