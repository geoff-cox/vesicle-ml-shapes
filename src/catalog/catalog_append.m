function T = catalog_append(simDir, hash, entry, ts)
% CATALOG_APPEND  Append/update SimResults/catalog.mat (MAT-only, UTC timestamps)
%
% T = catalog_append(simDir, hash, entry)       — load, append, save
% T = catalog_append(T,      hash, entry)       — in-memory only (no I/O)
%
% When the first argument is a table (returned from a previous call or
% catalog_load), the function works purely in memory.  The caller is
% responsible for calling catalog_save(simDir, T) when ready.
%
% When the first argument is a char/string directory path, the function
% loads from disk, appends, and writes back — identical to the original
% behaviour.

    inMemory = istable(simDir);

    % --- normalize incoming timestamp to UTC ---
    if nargin < 4 || isempty(ts)
        ts = datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss');
    else
        if ~isdatetime(ts), ts = datetime(ts); end
        if isempty(ts.TimeZone)
            ts.TimeZone = 'UTC';
        elseif ~strcmpi(ts.TimeZone,'UTC')
            try
                ts = datetime(ts,'TimeZone','UTC');
            catch
                ts.TimeZone = 'UTC';
            end
        end
        ts.Format = 'yyyy-MM-dd HH:mm:ss';
    end

    % --- load catalog (may be missing/empty) ---
    if inMemory
        T = simDir;           % first arg is already a table
    else
        T = catalog_load(simDir);
    end

    % --- ensure schema; when empty, make UTC-aware timestamp column ---
    if isempty(T) || ~any(T.Properties.VariableNames=="entry")
        dt = datetime.empty(0,1);
        dt.TimeZone = 'UTC';
        dt.Format   = 'yyyy-MM-dd HH:mm:ss';
        T = table(string.empty(0,1), dt, cell.empty(0,1), ...
                  'VariableNames', ["hash","timestamp","entry"]);
    end

    % --- normalize existing timestamp column to UTC (handles legacy) ---
    if ~isempty(T.timestamp)
        if isempty(T.timestamp.TimeZone)
            T.timestamp.TimeZone = 'UTC';
        elseif ~strcmpi(T.timestamp.TimeZone,'UTC')
            try
                T.timestamp = datetime(T.timestamp,'TimeZone','UTC');
            catch
                T.timestamp.TimeZone = 'UTC';
            end
        end
        T.timestamp.Format = 'yyyy-MM-dd HH:mm:ss';
    end

    % --- append or merge ---
    k = find(T.hash == string(hash), 1, 'first');
    if isempty(k)
        newRow = table(string(hash), ts, {entry}, ...
                       'VariableNames', ["hash","timestamp","entry"]);
        T = [T; newRow];   % both sides UTC
    else
        T.entry{k} = merge_entry(T.entry{k}, entry);
        if nargin >= 4 && ~isempty(ts)
            T.timestamp(k) = ts;
        end
    end

    % --- persist to disk only when operating in file mode ---
    if ~inMemory
        catalog_save(simDir, T);
    end
end

% ---------- helpers ----------
function e = merge_entry(a,b)
    e = a;
    e.params = merge_struct(getfield_safe(a,'params'), getfield_safe(b,'params'));
    e.meta   = merge_struct(getfield_safe(a,'meta'),   getfield_safe(b,'meta'));
end

function out = merge_struct(A,B)
    if ~isstruct(A), A = struct(); end
    if ~isstruct(B), B = struct(); end
    out = struct();
    keys = unique([fieldnames(A); fieldnames(B)]);
    for i=1:numel(keys)
        k = keys{i};
        va = getfield_safe(A,k);
        vb = getfield_safe(B,k);
        if isempty_value(vb), out.(k)=va; else, out.(k)=vb; end
    end
end

function v = getfield_safe(S,k)
    if isstruct(S) && isfield(S,k), v = S.(k); else, v = []; end
end

function tf = isempty_value(v)
    tf = isempty(v) ...
      || (isstring(v) && all(strlength(v)==0)) ...
      || (ischar(v)    && isempty(strtrim(v))) ...
      || (isnumeric(v) && (isempty(v) || (isscalar(v) && isnan(v))));
end
