function hex = simpleDataHash(params, model_version)
    if nargin < 2 || isempty(model_version), model_version = "BVP-v3.1"; end
    key = struct();
    key.model_version = string(model_version);
    key.H0_1 = defaultArg(params,'H0_1',NaN);
    key.H0_2 = defaultArg(params,'H0_2',NaN);
    key.A    = defaultArg(params,'A',   NaN);
    key.V    = defaultArg(params,'V',   NaN);
    key.KA   = defaultArg(params,'KA',  NaN);
    key.KB   = defaultArg(params,'KB',  NaN);
    key.KG   = defaultArg(params,'KG',  NaN);
    % DO NOT include: delta, opts, timestamps, etc.

    % Stable, version-proof hashing
    hex = simpleDataHash_internal(key, 'SHA-256');  % your MessageDigest-backed helper
end


function hex = simpleDataHash_internal(x, method)
% SIMPLEDATAHASH  Compute a stable content hash for MATLAB data structures.
%   hex = simpleDataHash(x, method)
%   - x       : any MATLAB variable (struct, array, cell, etc.)
%   - method  : 'MD5', 'SHA-1', 'SHA-256', etc. (default = 'SHA-256')
%
% Notes:
%   • Uses getByteStreamFromArray for stable serialization.
%   • Does not depend on normalizeForHash (legacy removed).
%   • For structs, sorts fieldnames alphabetically for deterministic hashing.
%   • For nested structs/cells, applies recursively.

    if nargin < 2 || isempty(method)
        method = 'SHA-256';
    end

    % --- normalize for deterministic field order ---
    x = normalize_fields(x);

    % --- compute hash ---
    bytes = getByteStreamFromArray(x);
    md = java.security.MessageDigest.getInstance(method);
    md.update(uint8(bytes));  % safe for all MATLAB versions
    digest = typecast(md.digest(), 'uint8');
    hex = lower(reshape(dec2hex(digest)', 1, []));
end

% ---------- local helper ----------
function y = normalize_fields(x)
    if isstruct(x)
        % sort fields for consistent ordering
        f = sort(fieldnames(x));
        for i = 1:numel(f)
            [x.(f{i})] = normalize_fields(x.(f{i}));
        end
        y = orderfields(x, f);
    elseif iscell(x)
        y = cellfun(@normalize_fields, x, 'UniformOutput', false);
    else
        y = x;  % numeric, char, logical, etc.
    end
end