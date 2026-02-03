% File: tools/migrate_simresults.m
% Purpose:
%   - Consolidate scattered per-run result .mat files into SimResults/solutions/
%   - (Re)build SimResults/catalog.csv and SimResults/catalog.mat
%   - Leave sources in place (copy, don't move) so you can verify before pruning
%
% How to run:
%   >> cd(project_root)
%   >> addpath('tools');  % if needed
%   >> migrate_simresults

function migrate_simresults()
    root = pwd;
    simDir = fullfile(root,'SimResults');
    solDir = fullfile(simDir,'solutions');
    if ~exist(solDir,'dir'); mkdir(solDir); end

    % 1) discover candidate result files in legacy layout
    legacyGlobs = { ...
        fullfile(simDir,'*','results','*.mat'), ...
        fullfile(simDir,'*','MATS','*.mat') ...         % optional; skip if not results
    };
    files = {};
    for g = legacyGlobs
        files = [files; cellstr(string(safe_glob(g{1})))]; %#ok<AGROW>
    end
    files = unique(files);

    % Also include any results already in SimResults/solutions
    existingFlat = cellstr(string(safe_glob(fullfile(solDir,'*.mat'))));
    fprintf('[migrate] Found %d legacy candidates, %d in flat store already.\n', ...
        numel(files), numel(existingFlat));

    % 2) scan & index -> rows for catalog table
    rows = [];
    seen = containers.Map('KeyType','char','ValueType','logical');

    for i = 1:numel(files)
        src = files{i};
        try
            S = load(src);
        catch ME
            warning('[migrate] Skipping unreadable file: %s (%s)', src, ME.message);
            continue;
        end

        % Heuristics: find params/meta/result
        params = detect_params(S);
        meta   = detect_meta(S);
        result = detect_result(S);

        % Compute or reuse hash
        if isfield(meta,'hash') && ~isempty(meta.hash)
            hash = string(meta.hash);
        else
            % If old files lack model_version, mark as 'legacy'
            model_version = 'legacy';
            if isfield(meta,'version') && ~isempty(meta.version)
                model_version = meta.version;
            end
            hash = string(simpleDataHash(params, model_version));
            meta.hash = hash;
            meta.version = model_version;
        end

        if isKey(seen, char(hash))
            % Already indexed in this migration pass
            continue;
        end
        seen(char(hash)) = true;

        % 3) copy into flat store if missing (atomic write)
        dst = fullfile(solDir, hash + ".mat");
        if ~exist(dst,'file')
            tmp = dst + ".tmp";
            try
                save(tmp, 'result','meta','-v7.3');
                movefile(tmp, dst, 'f');
            catch ME
                if exist(tmp,'file'); delete(tmp); end
                warning('[migrate] Failed to write %s (%s)', dst, ME.message);
                continue;
            end
        end

        % 4) construct catalog row
        rows = [rows; catalog_row(hash, params, meta)]; %#ok<AGROW>

        if mod(i,50)==0
            fprintf('[migrate] Processed %d/%d files…\n', i, numel(files));
        end
    end

    % 5) Merge with any existing catalog.csv (dedupe by hash)
    catCsv = fullfile(simDir,'catalog.csv');
    Tnew   = rows2table(rows);
    if exist(catCsv,'file')
        Told = readtable(catCsv, 'TextType','string');
        Tall = outerjoin(Told, Tnew, 'Keys','hash','MergeKeys',true);
        % Prefer non-missing fields from Tnew when conflicts arise
        T = prefer_right_nonmissing(Tall);
    else
        T = Tnew;
    end

    % Sort by timestamp (oldest first) for stability
    if any(strcmp('timestamp', T.Properties.VariableNames))
        T = sortrows(T, 'timestamp');
    end

    % 6) Write catalog.csv and catalog.mat
    writetable(T, catCsv);
    save(fullfile(simDir,'catalog.mat'), 'T');

    fprintf('[migrate] Done. Catalog rows: %d. Flat solutions: %d\n', height(T), numel(safe_glob(fullfile(solDir,'*.mat'))));
end

% ---------- helpers ----------

function list = safe_glob(pattern)
    d = dir(pattern);
    list = fullfile({d.folder}, {d.name});
end

function params = detect_params(S)
    % Try common field names; fallback to empty struct
    candidates = {'params','p','parameters'};
    params = struct();
    for c = candidates
        if isfield(S, c{1})
            params = S.(c{1});
            return;
        end
    end
    % If params embedded inside meta
    if isfield(S,'meta') && isfield(S.meta,'params')
        params = S.meta.params;
    end
end

function meta = detect_meta(S)
    if isfield(S,'meta'); meta = S.meta; return; end
    meta = struct();
    % Try to salvage known scalars
    take = {'energy','pressure','residualMax','meshN','shapeLabel','version','hash'};
    for k = take
        if isfield(S, k{1}), meta.(k{1}) = S.(k{1}); end
    end
end

function result = detect_result(S)
    if isfield(S,'result'); result = S.result; return; end
    % fallback: pack common arrays if they exist
    result = struct();
    take = {'s','y','bcResidual','deResidual'};
    for k = take
        if isfield(S, k{1}), result.(k{1}) = S.(k{1}); end
    end
end

function h = simpleDataHash(params, model_version)
    % Keep only physics knobs; ignore solver counters/timestamps.
    key = struct('solverABI', string(model_version));
    % Copy selected fields if they exist
    keys = {'H0_1','H0_2','x1','v','k1','k2','kG'};
    for i=1:numel(keys)
        k = keys{i};
        if isfield(params,k), key.(k) = params.(k); end
    end
    bytes = getByteStreamFromArray(key);
    algorithm = 'MD5'; % Other options: SHA-256
    h = hash_bytes(bytes, algorithm);

end

function h = hash_bytes(bytes, algo)
    if nargin < 2, algo = 'MD5'; end
    md = java.security.MessageDigest.getInstance(algo);
    md.update(uint8(bytes));
    raw = typecast(md.digest, 'uint8');
    h = lower(reshape(dec2hex(raw)',1,[]));
end

function row = catalog_row(hash, params, meta)
    ts = datetime('now','Timezone','UTC','Format','yyyy-MM-dd HH:mm:ss');
    row = struct( ...
        'hash',         string(hash), ...
        'timestamp',    ts, ...
        'H0_1',         field_or_nan(params,'H0_1'), ...
        'H0_2',         field_or_nan(params,'H0_2'), ...
        'x1',           field_or_nan(params,'x1'), ...
        'v',            field_or_nan(params,'v'), ...
        'k1',           field_or_nan(params,'k1'), ...
        'k2',           field_or_nan(params,'k2'), ...
        'kG',           field_or_nan(params,'kG'), ...
        'energy',       field_or_nan(meta,'energy'), ...
        'pressure',     field_or_nan(meta,'pressure'), ...
        'residualMax',  field_or_nan(meta,'residualMax'), ...
        'meshN',        field_or_nan(meta,'meshN'), ...
        'shape',        field_or_string(meta,'shapeLabel'), ...
        'model_version',field_or_string(meta,'version') ...
    );
end

function v = field_or_nan(S, k)
    if isfield(S,k) && ~isempty(S.(k)) && isnumeric(S.(k)), v = double(S.(k));
    else, v = NaN;
    end
end

function v = field_or_string(S,k)
    if isfield(S,k) && ~isempty(S.(k)), v = string(S.(k));
    else, v = "";
    end
end

function T = rows2table(rows)
    if isempty(rows), T = cell2table(cell(0,15), 'VariableNames', ...
        ["hash","timestamp","H0_1","H0_2","x1","v","k1","k2","kG","energy","pressure","residualMax","meshN","shape","model_version"]); return; end
    T = struct2table(rows);
    % Ensure types
    T.hash = string(T.hash);
    if ~isa(T.timestamp,'datetime'), T.timestamp = datetime(T.timestamp); end
end

function T = prefer_right_nonmissing(Tall)
    % Tall has variables like hash, timestamp_Told, timestamp_Tnew
    % Keep 'hash' and coalesce *_Tnew over *_Told when not missing
    vars = Tall.Properties.VariableNames;
    base = unique(erase(vars(contains(vars,'_Told') | contains(vars,'_Tnew')), {'_Told','_Tnew'}), 'stable');
    out = table();
    out.hash = Tall.hash;
    for i=1:numel(base)
        b = base{i};
        vNew = Tall.(b + "_Tnew");
        vOld = Tall.(b + "_Told");
        if isempty(vOld); vOld = repmat(missing, height(Tall),1); end
        if isempty(vNew); vNew = repmat(missing, height(Tall),1); end
        v = vNew;
        m = ismissing(v);
        v(m) = vOld(m);
        out.(b) = v;
    end
    T = unique(out, 'rows', 'stable');
end
