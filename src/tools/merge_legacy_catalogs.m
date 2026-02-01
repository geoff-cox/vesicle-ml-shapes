% File: tools/merge_legacy_catalogs.m
% Merge all legacy sim-results/**/catalog.mat into sim-results/catalog.mat.

function merge_legacy_catalogs()
    root  = pwd;
    simDir = fullfile(root,'sim-results');
    solDir = fullfile(simDir,'solutions');
    if ~exist(simDir,'dir'); error('Missing sim-results'); end
    if ~exist(solDir,'dir'); mkdir(solDir); end

    addpath(fullfile(root,'src','utils')); %#ok<ADD> for catalog_mat.m

    % 1) Start from existing unified catalog (if present)
    Tall = catalog_load(simDir);

    % 2) Find all legacy catalog.mat files (excluding the unified one)
    cats = dir(fullfile(simDir,'**','catalog.mat'));
    cats = cats(~strcmp(fullfile({cats.folder},{cats.name})', fullfile(simDir,'catalog.mat')));

    fprintf('[merge] Found %d legacy catalog.mat files.\n', numel(cats));

    for i = 1:numel(cats)
        f = fullfile(cats(i).folder, cats(i).name);
        try
            S = load(f);
        catch ME
            warning('[merge] Skipping unreadable %s (%s)', f, ME.message);
            continue;
        end
        if ~isfield(S,'T')
            warning('[merge] %s has no variable T; skipping.', f);
            continue;
        end
        Tlegacy = S.T;
        % Normalize shape to (hash,timestamp,entry)
        Tlegacy = normalize_to_entry_table(Tlegacy);

        % Coalesce into Tall
        for r = 1:height(Tlegacy)
            hash = string(Tlegacy.hash(r));
            entry = Tlegacy.entry(r);
            ts    = Tlegacy.timestamp(r);
            Tall  = catalog_append(simDir, hash, entry, ts); % append does save; but we hold on to Tall ref
        end
        fprintf('[merge] Merged %4d rows from %s\n', height(Tlegacy), short_path(f, simDir));
    end

    % 3) Optional: prune entries with missing solution files
    Tall = catalog_load(simDir); % re-load final after multiple appends
    keep = true(height(Tall),1);
    for r = 1:height(Tall)
        fsol = fullfile(solDir, Tall.hash(r) + ".mat");
        if ~exist(fsol,'file')
            warning('[merge] Missing solution for hash=%s; keeping row but consider rebuilding or deleting.', Tall.hash(r));
            % If you prefer to drop: keep(r) = false;
        end
    end
    % Tall = Tall(keep,:); catalog_save(simDir, Tall);

    fprintf('[merge] Done. Unified catalog rows: %d (file: %s)\n', height(Tall), fullfile(simDir,'catalog.mat'));
end

% -------------- helpers --------------

function T = normalize_to_entry_table(Tin)
    if istable(Tin) && any(Tin.Properties.VariableNames == "entry")
        % Already normalized
        T = Tin;
        % Type fixes
        T.hash = string(T.hash);
        if ~isduration(T.timestamp) && ~isdatetime(T.timestamp)
            T.timestamp = datetime(T.timestamp);
        end
        return;
    end
    % Assume legacy wide schema; fold to entry
    T = fold_legacy_columns_local(Tin);
end

function T = fold_legacy_columns_local(Tin)
    names = Tin.Properties.VariableNames;
    assert(any(strcmp(names,'hash')) && any(strcmp(names,'timestamp')), 'Legacy catalog missing hash or timestamp');
    T = table; T.hash = string(Tin.hash); T.timestamp = Tin.timestamp;
    E = repmat(struct('params',struct(),'meta',struct()), height(Tin), 1);
    paramKeys = ["H0_1","H0_2","x1","v","k1","k2","kG"];
    for r=1:height(Tin)
        P = struct(); M = struct();
        for j=1:numel(names)
            nm = names{j};
            if any(strcmp(nm, {'hash','timestamp'})), continue; end
            val = Tin.(nm)(r);
            if any(paramKeys == string(nm))
                P.(nm) = val;
            else
                switch nm
                    case 'shape', M.shapeLabel = val;
                    case 'model_version', M.version = val;
                    otherwise, M.(nm) = val;
                end
            end
        end
        E(r).params = P;
        E(r).meta   = M;
    end
    T.entry = E;
end

function s = short_path(p, base)
    s = erase(string(p), string(base)+filesep);
end
