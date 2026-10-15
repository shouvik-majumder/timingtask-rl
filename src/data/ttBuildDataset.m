function D = ttBuildDataset(root, animals, varargin)
%TTBUILDDATASET  Dataset structure from a folder of processed session files.
%
%   D = ttBuildDataset(root, animals)
%   D = ttBuildDataset(root, animals, 'groups', G, 'minTrials', 50, 'pattern', '*.mat', 'stageOptions', {...})
%
%   root      folder with one subfolder per animal, each holding one .mat file per session
%             (variable S, see ttLoadSession); the folder name may differ from the animal name in
%             letter case
%   animals   cell array or string array of animal names
%   groups    optional table with columns animal and group (and any others, kept in D.sessions)
%   Sessions are ordered by file name (the laboratory names files <animal>_<yymmdd>_<HHMMSS>).
%   Sessions with fewer than minTrials cued trials are dropped.
%
%   D.trials    one row per trial (ttLoadSession columns plus sessionIndex, stage, phase)
%   D.sessions  one row per session: animal, session, sessionIndex, date, nTrials, nCued, nRewarded,
%               nNoCue, delayLevels, delayMean, stage, phase, stageIndex (position within the
%               animal's sessions of that stage, from 1), stageIndexFromEnd (from the last session
%               of the stage, -1 = last), group
%   Phases are assigned by ttAssignPhases with its default split.
p = inputParser; p.addParameter('groups', table(), @istable); p.addParameter('minTrials', 50, @isnumeric); p.addParameter('pattern', '*.mat', @ischar);
p.addParameter('stageOptions', {}, @iscell); p.addParameter('verbose', true, @islogical); p.parse(varargin{:}); o = p.Results;
animals = string(animals(:))'; trials = {}; sessions = {};
for a = animals
    folder = fullfile(root, char(a));
    if ~isfolder(folder)
        cand = dir(root); cand = cand([cand.isdir]); hit = cand(strcmpi({cand.name}, char(a)));
        if isempty(hit), if o.verbose, fprintf('ttBuildDataset: no folder for %s\n', a); end; continue; end
        folder = fullfile(root, hit(1).name);
    end
    files = dir(fullfile(folder, o.pattern)); [~, ord] = sort({files.name}); files = files(ord); k = 0;
    for f = 1:numel(files)
        try, L = load(fullfile(files(f).folder, files(f).name), 'S'); catch, if o.verbose, fprintf('ttBuildDataset: cannot load %s\n', files(f).name); end; continue; end
        [~, sid] = fileparts(files(f).name); T = ttLoadSession(L.S, 'animal', a, 'session', sid);
        if sum(T.cued) < o.minTrials, continue; end
        k = k + 1; T.sessionIndex = repmat(k, height(T), 1);
        dc = T.delay(T.cued); st = ttClassifyStage(dc, height(T), o.stageOptions{:});
        tok = regexp(sid, '_(\d{6})_(\d{6})$', 'tokens', 'once'); dt = NaT; if ~isempty(tok), dt = datetime([tok{1}, tok{2}], 'InputFormat', 'yyMMddHHmmss'); end
        lv = unique(round(dc, 2)); lv = lv(lv > 0.05);
        sessions{end+1} = table(a, string(sid), k, dt, height(T), sum(T.cued), sum(T.cued & T.rewarded), sum(T.noCue), string(mat2str(lv', 3)), mean(dc), st, ...
            'VariableNames', {'animal', 'session', 'sessionIndex', 'date', 'nTrials', 'nCued', 'nRewarded', 'nNoCue', 'delayLevels', 'delayMean', 'stage'}); %#ok<AGROW>
        T.stage = repmat(st, height(T), 1); trials{end+1} = T; %#ok<AGROW>
    end
    if o.verbose, fprintf('ttBuildDataset: %s, %d sessions\n', a, k); end
end
D = struct('trials', vertcat(trials{:}), 'sessions', vertcat(sessions{:}));
D.sessions.group = repmat("", height(D.sessions), 1);
if ~isempty(o.groups)
    G = o.groups; G.animal = string(G.animal); [tf, loc] = ismember(lower(D.sessions.animal), lower(G.animal));
    for c = G.Properties.VariableNames, if strcmp(c{1}, 'animal'), continue; end; v = G.(c{1}); if isnumeric(v), col = nan(height(D.sessions), 1); else, col = repmat("", height(D.sessions), 1); v = string(v); end; col(tf) = v(loc(tf)); D.sessions.(c{1}) = col; end
end
D = ttAssignPhases(D);
end
