%EXAMPLE_DATASET  Fit the core models to a folder of recorded sessions and compare two groups of animals.
%
%   Expected layout: <dataRoot>/<animal>/<animal>_<yymmdd>_<HHMMSS>.mat, each file holding the
%   laboratory's session structure S (ttLoadSession). The group table assigns animals to groups.
root = ttSetup(); outDir = fullfile(root, 'output', 'dataset'); if ~exist(outDir, 'dir'), mkdir(outDir); end
dataRoot = fullfile(root, 'data', 'sessions');              % replace with the data folder
groups = table(["M01"; "M02"; "M03"; "M04"], ["mutant"; "control"; "mutant"; "control"], 'VariableNames', {'animal', 'group'});
if isempty(gcp('nocreate')), parpool; end

D = ttBuildDataset(dataRoot, groups.animal, 'groups', groups);
disp(groupsummary(D.sessions, {'animal', 'stage'}));

ttDistCV(D, outDir); ttDistSummary(outDir);
ttFitAll(D, ttModelSet('core'), outDir);
ttHeldOut(D, outDir);
C = ttCompare(outDir);
ttUpdateTest(outDir);

P = ttParamTable(outDir, 'lnexp_ThistoryLearn_Efixed_Rthreshold');
G = ttGroupCompare(P, 'alpha_up', 'mutant', 'control', 'by', 'phase', 'transform', @log10);
disp(G);
