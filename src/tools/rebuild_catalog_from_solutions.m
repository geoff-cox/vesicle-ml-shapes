% File: tools/rebuild_catalog_from_solutions.m
% Purpose:
%   Rebuild sim-results/catalog.csv and sim-results/catalog.mat by scanning
%   sim-results/solutions/*.mat, using the filename stem as the hash.
%
% How to run:
%   >> cd(project_root)
%   >> addpath('tools');
%   >> rebuild_catalog_from_solutions

function rebuild_catalog_from_solutions()
    simDir = fullfile(pwd,'sim-results');
    solDir = fullfile(simDir,'solutions');
    if ~exist(solDir,'dir')
        error('Missing folder: %s', solDir);
    end

    files = dir(fullfile(solDir,'*.mat'));
    if isempty(files)
        warning('No .mat files found in %s', solDir);
    end

    rows = repmat(empty_row(), 0, 1);

    for i = 1:numel(files)
        f = fullfile(files(i).folder, files(i).name);
        [~, base, ~] = fileparts(f);
        hash = string(base);  % source of truth

        % Load once, but be resilient to odd files
        try
            S = load(f);
        catch ME
            warning('[rebuild] Skipping unreadable file: %s (%s)', f, ME.message);
            continue;
        end

        % Recover params/meta/result if present
        params = detect_params(S);
        meta   = detect_meta(S);
        result = detect_result(S);

        % If meta.hash/version missing, set from context
        if ~isfield(meta,'hash') || strlength(string(meta.hash))==0
            meta.hash = hash;
        end
        if ~isfield(meta,'version') || strlength(string(meta.version))==0
            meta.version = "legacy";
        end

        % Derive quick stats if missing (residual, meshN, etc.)
        meta = derive_meta_stats(meta, result);

        % Timestamp = file mtime (stable, deterministic)
        ts = datetime(files(i).datenum, 'ConvertFrom','datenum', ...
                      'TimeZone','UTC', 'Format','yyyy-MM-dd HH:mm:ss');

        % Build catalog row
        row = struct( ...
            'hash',          hash, ...
            'timestamp',     ts, ...
            'H0_1',          field_or_nan(params,'H0_1'), ...
            'H0_2',          field_or_nan(params,'H0_2'), ...
            'x1',            field_or_nan(params,'x1'), ...
            'v',             field_or_nan(params,'v'), ...
            'k1',            field_or_nan(params,'k1'), ...
            'k2',            field_or_nan(params,'k2'), ...
            'kG',            field_or_nan(params,'kG'), ...
            'energy',        field_or_nan(meta,'energy'), ...
            'pressure',      field_or_nan(meta,'pressure'), ...
            'residualMax',   field_or_nan(meta,'residualMax'), ...
            'meshN',         field_or_nan(meta,'meshN'), ...
            'shape',         field_or_string(meta,'shapeLabel'), ...
            'model_version', field_or_string(meta,'version') ...
        );

        rows(end+1,1) = row; %#ok<AGROW>

        if mod(i,100)==0
            fprintf('[rebuild] Processed %d/%d files…\n', i, numel(files));
        end
    end

    % Make table, sort by timestamp, drop dup hashes (keep newest mtime)
    T = struct2table(rows);
    if ~isempty(T)
        [~, keepIdx] = unique(T.hash, 'last');  % if duplicates, keep last (most recent mtime)
        T = T(sort(keepIdx), :);
        T = sortrows(T, 'timestamp');
    else
        % Ensure correct variable names/types even if empty
        T = cell2table(cell(0,15), 'VariableNames', ...
             ["hash","timestamp","H0_1","H0_2","x1","v","k1","k2","kG","energy","pressure","residualMax","meshN","shape","model_version"]);
    end

    % Write outputs
    catCsv = fullfile(simDir,'catalog.csv');
    catMat = fullfile(simDir,'catalog.mat');

    writetable(T, catCsv);
    save(catMat, 'T');

    fprintf('[rebuild] Done. Wrote %d rows to:\n  %s\n  %s\n', height(T), catCsv, catMat);
end

% ---------------- helpers ----------------

function params = detect_params(S)
    candidates = {'params','p','parameters'};
    params = struct();
    for c = candidates
        if isfield(S, c{1})
            params = S.(c{1});
            return;
        end
    end
    if isfield(S,'meta') && isfield(S.meta,'params')
        params = S.meta.params;
    end
end

function meta = detect_meta(S)
    if isfield(S,'meta')
        meta = S.meta;
        return
    end
    meta = struct();
    % Scavenge common fields if they were saved flat
    take = {'energy','pressure','residualMax','meshN','shapeLabel','version','hash'};
    for k = take
        if isfield(S, k{1}), meta.(k{1}) = S.(k{1}); end
    end
end

function result = detect_result(S)
    if isfield(S,'result')
        result = S.result;
        return
    end
    result = struct();
    % Fallback: pack common arrays if they exist
    take = {'s','y','bcResidual','deResidual'};
    for k = take
        if isfield(S, k{1}), result.(k{1}) = S.(k{1}); end
    end
end

function meta = derive_meta_stats(meta, result)
    % meshN
    if ~isfield(meta,'meshN') || isempty(meta.meshN) || ~isnumeric(meta.meshN)
        if isfield(result,'s') && ~isempty(result.s)
            meta.meshN = numel(result.s);
        elseif isfield(result,'y') && ~isempty(result.y)
            % y can be [stateDim x N] or [N x stateDim]
            dims = size(result.y);
            meta.meshN = max(dims);
        else
            meta.meshN = NaN;
        end
    end

    % residualMax
    if ~isfield(meta,'residualMax') || isempty(meta.residualMax) || ~isnumeric(meta.residualMax)
        r = [];
        if isfield(result,'bcResidual') && ~isempty(result.bcResidual)
            r = [r; result.bcResidual(:)];
        end
        if isfield(result,'deResidual') && ~isempty(result.deResidual)
            r = [r; result.deResidual(:)];
        end
        if ~isempty(r)
            meta.residualMax = max(abs(r));
        else
            meta.residualMax = NaN;
        end
    end

    % Default placeholders if energy/pressure/shapeLabel missing
    if ~isfield(meta,'energy'),   meta.energy = NaN;     end
    if ~isfield(meta,'pressure'), meta.pressure = NaN;   end
    if ~isfield(meta,'shapeLabel') || isempty(meta.shapeLabel)
        meta.shapeLabel = "";
    end
end

function v = field_or_nan(S, k)
    if isfield(S,k) && ~isempty(S.(k)) && isnumeric(S.(k))
        v = double(S.(k));
    else
        v = NaN;
    end
end

function v = field_or_string(S, k)
    if isfield(S,k) && ~isempty(S.(k))
        v = string(S.(k));
    else
        v = "";
    end
end

function r = empty_row()
    r = struct( ...
        'hash',          "", ...
        'timestamp',     datetime.empty, ...
        'H0_1',          NaN, ...
        'H0_2',          NaN, ...
        'x1',            NaN, ...
        'v',             NaN, ...
        'k1',            NaN, ...
        'k2',            NaN, ...
        'kG',            NaN, ...
        'energy',        NaN, ...
        'pressure',      NaN, ...
        'residualMax',   NaN, ...
        'meshN',         NaN, ...
        'shape',         "", ...
        'model_version', "" ...
    );
end
