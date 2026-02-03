% File: tools/rebuild_catalog_from_solutions.m
% Purpose:
%   Rebuild SimResults/catalog.mat (canonical entry-schema) by scanning
%   SimResults/hashed_results/*.mat and SimResults/solutions/*.mat.
%   Also write a wide CSV for analytics: SimResults/catalog.csv (and catalog_wide.csv).
%
% Canonical schema:
%   T = table(hash:string, timestamp:datetime(UTC), entry:cell)
%   where entry{i} = struct('params',struct,'meta',struct)
%
% How to run:
%   >> cd(project_root)
%   >> addpath('tools');
%   >> rebuild_catalog_from_solutions

function rebuild_catalog_from_solutions()

    projRoot = pwd;
    simDir   = fullfile(projRoot,'SimResults');

    if ~exist(simDir,'dir')
        error('Missing folder: %s', simDir);
    end

    % Candidate result folders
    dirA = fullfile(simDir,'hashed_results');
    dirB = fullfile(simDir,'solutions');

    files = [];
    if exist(dirA,'dir')
        files = [files; dir(fullfile(dirA,'*.mat'))]; %#ok<AGROW>
    end
    if exist(dirB,'dir')
        files = [files; dir(fullfile(dirB,'*.mat'))]; %#ok<AGROW>
    end

    if isempty(files)
        warning('No .mat result files found in %s or %s', dirA, dirB);
    end

    % Load existing canonical catalog if present (to reuse params when missing in files)
    Tprev = try_load_existing_catalog(simDir);
    prevMap = build_params_map(Tprev);  % containers.Map(hash -> params)

    % Collect canonical rows as struct array, then convert to table
    rows = repmat(struct('hash',"",'timestamp',datetime.empty,'entry',{{}}), 0, 1);

    for i = 1:numel(files)
        f = fullfile(files(i).folder, files(i).name);
        [~, base, ~] = fileparts(f);
        hash = string(base);

        % Skip if it looks like a catalog.mat or other non-hash file
        if hash == "catalog" || hash == "catalog_wide" || hash == "catalog_old"
            continue;
        end

        % Load file resiliently
        try
            S = load(f);
        catch ME
            warning('[rebuild] Skipping unreadable file: %s (%s)', f, ME.message);
            continue;
        end

        % Timestamp = file mtime
        ts = datetime(files(i).datenum, 'ConvertFrom','datenum', ...
                      'TimeZone','UTC', 'Format','yyyy-MM-dd HH:mm:ss');

        % Extract/derive params/meta/result
        [params, meta, result] = detect_payload(S);

        % If params missing, try reuse from previous catalog
        if isempty(fieldnames(params)) && isKey(prevMap, char(hash))
            params = prevMap(char(hash));
        end

        % Make sure meta has hash/version
        if ~isfield(meta,'hash') || strlength(string(meta.hash))==0
            meta.hash = hash;
        end
        if ~isfield(meta,'version') || strlength(string(meta.version))==0
            meta.version = "legacy";
        end

        % Derive common meta stats (mesh, residualMax, etc.)
        meta = derive_meta_stats(meta, result);

        entry = struct('params', params, 'meta', meta);

        rows(end+1,1) = struct('hash',hash,'timestamp',ts,'entry',{{entry}}); %#ok<AGROW>

        if mod(i,200)==0
            fprintf('[rebuild] Processed %d/%d files...\n', i, numel(files));
        end
    end

    % Canonical table
    T = canonical_table_from_rows(rows);

    % Drop duplicate hashes (keep newest timestamp)
    if ~isempty(T)
        [~, keepIdx] = unique(T.hash, 'last');
        T = T(sort(keepIdx), :);
        T = sortrows(T, 'timestamp');
    end

    % Save canonical catalog.mat atomically
    catalog_save(simDir, T);

    % Also write wide analytics CSV(s)
    Twide = make_wide_table(T);
    catCsvA = fullfile(simDir,'catalog.csv');        % keep expected name for analytics
    catCsvB = fullfile(simDir,'catalog_wide.csv');   % explicit name
    writetable(Twide, catCsvA);
    writetable(Twide, catCsvB);

    fprintf('[rebuild] Done.\n  canonical: %s\n  wide csv:   %s\n  wide csv:   %s\n  rows: %d\n', ...
        fullfile(simDir,'catalog.mat'), catCsvA, catCsvB, height(T));
end

% ===================== helpers =====================

function T = try_load_existing_catalog(simDir)
    catF = fullfile(simDir,'catalog.mat');
    if ~exist(catF,'file')
        T = table(string.empty(0,1), datetime.empty(0,1), cell.empty(0,1), ...
                  'VariableNames', ["hash","timestamp","entry"]);
        return;
    end
    try
        S = load(catF,'T');
        T = S.T;
    catch
        T = table(string.empty(0,1), datetime.empty(0,1), cell.empty(0,1), ...
                  'VariableNames', ["hash","timestamp","entry"]);
        return;
    end

    % If it's wide/legacy (no entry), ignore for params reuse
    if isempty(T) || ~any(T.Properties.VariableNames=="entry")
        T = table(string.empty(0,1), datetime.empty(0,1), cell.empty(0,1), ...
                  'VariableNames', ["hash","timestamp","entry"]);
        return;
    end

    % Enforce cell entry column
    if isstruct(T.entry)
        T.entry = num2cell(T.entry);
    end

    % Normalize timestamp timezone
    if ~isempty(T.timestamp)
        if isempty(T.timestamp.TimeZone), T.timestamp.TimeZone = 'UTC'; end
        T.timestamp.Format = 'yyyy-MM-dd HH:mm:ss';
    end
end

function mp = build_params_map(T)
    mp = containers.Map('KeyType','char','ValueType','any');
    if isempty(T) || ~any(T.Properties.VariableNames=="entry"), return; end
    for i=1:height(T)
        h = char(T.hash(i));
        e = T.entry{i};
        if isstruct(e) && isfield(e,'params') && isstruct(e.params)
            mp(h) = e.params;
        end
    end
end

function [params, meta, result] = detect_payload(S)
    params = struct();
    meta   = struct();
    result = struct();

    % Canonical-ish
    if isfield(S,'params') && isstruct(S.params), params = S.params; end
    if isfield(S,'meta')   && isstruct(S.meta),   meta   = S.meta;   end
    if isfield(S,'result') && isstruct(S.result), result = S.result; end

    % Sometimes an entry struct is saved
    if isempty(fieldnames(params)) && isfield(S,'entry') && isstruct(S.entry)
        if isfield(S.entry,'params') && isstruct(S.entry.params), params = S.entry.params; end
        if isfield(S.entry,'meta')   && isstruct(S.entry.meta),   meta   = S.entry.meta;   end
    end

    % Sometimes meta.params exists
    if isempty(fieldnames(params)) && isfield(meta,'params') && isstruct(meta.params)
        params = meta.params;
    end

    % Scavenge common flat fields into meta if meta empty
    if isempty(fieldnames(meta))
        take = {'label','E','P','BCmax','DEmax','mesh','meshN','residualMax','rMinAway','version','hash'};
        for k=1:numel(take)
            nm = take{k};
            if isfield(S,nm), meta.(nm) = S.(nm); end
        end
    end

    % Fallback: if result is missing but variables exist, pack minimal
    if isempty(fieldnames(result))
        if isfield(S,'sol'), result.sol = S.sol; end
        if isfield(S,'s'),   result.s   = S.s;   end
        if isfield(S,'y'),   result.y   = S.y;   end
        if isfield(S,'bcResidual'), result.bcResidual = S.bcResidual; end
        if isfield(S,'deResidual'), result.deResidual = S.deResidual; end
    end
end

function meta = derive_meta_stats(meta, result)
    % mesh / meshN
    if ~isfield(meta,'mesh') || isempty(meta.mesh) || ~isnumeric(meta.mesh)
        m = NaN;
        if isfield(result,'sol') && isstruct(result.sol) && isfield(result.sol,'x')
            m = numel(result.sol.x);
        elseif isfield(result,'s') && ~isempty(result.s)
            m = numel(result.s);
        elseif isfield(result,'y') && ~isempty(result.y)
            m = max(size(result.y));
        end
        meta.mesh = m;
    end

    if ~isfield(meta,'meshN') || isempty(meta.meshN) || ~isnumeric(meta.meshN)
        meta.meshN = meta.mesh;
    end

    % residualMax from stored residual vectors if present
    if ~isfield(meta,'residualMax') || isempty(meta.residualMax) || ~isnumeric(meta.residualMax)
        r = [];
        if isfield(result,'bcResidual') && ~isempty(result.bcResidual), r = [r; result.bcResidual(:)]; end %#ok<AGROW>
        if isfield(result,'deResidual') && ~isempty(result.deResidual), r = [r; result.deResidual(:)]; end %#ok<AGROW>
        if ~isempty(r)
            meta.residualMax = max(abs(r));
        else
            meta.residualMax = NaN;
        end
    end

    % Ensure canonical-ish names exist
    if ~isfield(meta,'label'), meta.label = ""; end
    if ~isfield(meta,'E'),     meta.E     = NaN; end
    if ~isfield(meta,'P'),     meta.P     = NaN; end
    if ~isfield(meta,'BCmax'), meta.BCmax = NaN; end
    if ~isfield(meta,'DEmax'), meta.DEmax = NaN; end
    if ~isfield(meta,'rMinAway'), meta.rMinAway = NaN; end
end

function T = canonical_table_from_rows(rows)
    % Create canonical table even if empty
    if isempty(rows)
        dt = datetime.empty(0,1);
        dt.TimeZone = 'UTC';
        dt.Format   = 'yyyy-MM-dd HH:mm:ss';
        T = table(string.empty(0,1), dt, cell.empty(0,1), ...
                  'VariableNames', ["hash","timestamp","entry"]);
        return;
    end

    hash = string({rows.hash})';
    ts   = vertcat(rows.timestamp);
    ent  = vertcat(rows.entry);  % already cell-of-entry
    T = table(hash, ts, ent, 'VariableNames', ["hash","timestamp","entry"]);

    % enforce timestamp UTC formatting
    if isempty(T.timestamp.TimeZone), T.timestamp.TimeZone = 'UTC'; end
    T.timestamp.Format = 'yyyy-MM-dd HH:mm:ss';

    % enforce entry cell type
    if isstruct(T.entry)
        T.entry = num2cell(T.entry);
    end
end

function Tw = make_wide_table(T)
    % Wide analytics view for CSV; safe even if params/meta are incomplete.
    n = height(T);

    % Preallocate columns
    hash      = T.hash;
    timestamp = T.timestamp;

    H0_1 = nan(n,1); H0_2 = nan(n,1);
    A    = nan(n,1); V    = nan(n,1);
    KA   = nan(n,1); KB   = nan(n,1); KG = nan(n,1);

    label = strings(n,1);
    E     = nan(n,1); P    = nan(n,1);
    BCmax = nan(n,1); DEmax = nan(n,1);
    mesh  = nan(n,1);
    rMinAway   = nan(n,1);
    residualMax = nan(n,1);
    version = strings(n,1);

    for i=1:n
        e = T.entry{i};
        if ~isstruct(e), continue; end

        p = struct(); m = struct();
        if isfield(e,'params') && isstruct(e.params), p = e.params; end
        if isfield(e,'meta')   && isstruct(e.meta),   m = e.meta;   end

        H0_1(i) = field_or_nan(p,'H0_1');
        H0_2(i) = field_or_nan(p,'H0_2');
        A(i)    = field_or_nan(p,'A');
        V(i)    = field_or_nan(p,'V');
        KA(i)   = field_or_nan(p,'KA');
        KB(i)   = field_or_nan(p,'KB');
        KG(i)   = field_or_nan(p,'KG');

        label(i) = string(field_or_string(m,'label'));
        E(i)     = field_or_nan(m,'E');
        P(i)     = field_or_nan(m,'P');
        BCmax(i) = field_or_nan(m,'BCmax');
        DEmax(i) = field_or_nan(m,'DEmax');
        mesh(i)  = field_or_nan(m,'mesh');
        rMinAway(i) = field_or_nan(m,'rMinAway');
        residualMax(i) = field_or_nan(m,'residualMax');
        version(i) = string(field_or_string(m,'version'));
    end

    Tw = table(hash, timestamp, H0_1, H0_2, A, V, KA, KB, KG, ...
               label, E, P, BCmax, DEmax, mesh, rMinAway, residualMax, version);
end

function v = field_or_nan(S, k)
    if isstruct(S) && isfield(S,k) && ~isempty(S.(k)) && isnumeric(S.(k))
        v = double(S.(k));
    else
        v = NaN;
    end
end

function v = field_or_string(S, k)
    if isstruct(S) && isfield(S,k) && ~isempty(S.(k))
        v = string(S.(k));
    else
        v = "";
    end
end

function catalog_save(simDir, T)
% Atomic write of canonical catalog.mat
    f   = fullfile(simDir, 'catalog.mat');
    tmp = [f '.tmp'];
    save(tmp, 'T', '-v7');   % use -v7.3 if the table will exceed 2 GB
    movefile(tmp, f, 'f');
end