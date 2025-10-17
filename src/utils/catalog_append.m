function T = catalog_append(simDir, hash, entry, ts)
    if nargin < 4 || isempty(ts)
        ts = datetime('now','Timezone','UTC','Format','yyyy-MM-dd HH:mm:ss');
    end
    T = catalog_load(simDir);
    k = find(T.hash == string(hash), 1, 'first');
    if isempty(k)
        newRow = table(string(hash), ts, entry, 'VariableNames', ["hash","timestamp","entry"]);
        T = [T; newRow];
    else
        T.entry(k) = merge_entry(T.entry(k), entry);
        % Keep earliest timestamp unless explicit ts is provided
        if nargin >= 4 && ~isempty(ts)
            T.timestamp(k) = ts;
        end
    end
    catalog_save(simDir, T);
end

% ---------- helpers ----------

function e = merge_entry(a, b)
    e = a;
    e.params = merge_struct(getfield(a,'params'), getfield(b,'params')); %#ok<GFLD>
    e.meta   = merge_struct(getfield(a,'meta'),   getfield(b,'meta'));
end

function out = merge_struct(A, B)
    if ~isstruct(A), A = struct(); end
    if ~isstruct(B), B = struct(); end
    keys = unique([fieldnames(A); fieldnames(B)]);
    for i=1:numel(keys)
        k = keys{i};
        vb = []; if isfield(B,k), vb = B.(k); end
        if isempty_value(vb), out.(k) = getfield(A,k); %#ok<GFLD>
        else,                  out.(k) = vb;
        end
    end
end

function tf = isempty_value(v)
    tf = isempty(v) || (isstring(v) && strlength(v)==0) || (ischar(v) && isempty(strtrim(v))) ...
        || (isnumeric(v) && (isempty(v) || (isscalar(v) && isnan(v))));
end