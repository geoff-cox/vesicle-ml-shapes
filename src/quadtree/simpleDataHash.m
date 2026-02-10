function hex = simpleDataHash(params, arg2, arg3)
    % SIMPLEDATAHASH  Stable content hash for physics-defining parameters.
    %
    % Supported call patterns:
    %   hex = simpleDataHash(params)
    %   hex = simpleDataHash(params, model_version)
    %   hex = simpleDataHash(params, model_version, method)
    %   hex = simpleDataHash(params, method)              % backward-compat for mistaken usage
    %
    % Rules:
    %   - Hash includes: model_version, H0_1, H0_2, A, V, KA, KB, KG
    %   - Does NOT include: delta, opts, timestamps, etc.
    %
    % If model_version is not provided, we try:
    %   params.model_version, params.version, else default "BVP-v3.1".

    if nargin < 2, arg2 = []; end
    if nargin < 3, arg3 = []; end

    allowed = upper(string(["MD5","SHA-1","SHA-256","SHA-384","SHA-512"]));

    % ---- interpret arguments ----
    method = "SHA-256";
    model_version = "";

    if nargin == 1
        % nothing
    elseif nargin == 2
        s2 = string(arg2);
        if any(upper(s2) == allowed)
            method = upper(s2);
        else
            model_version = s2;
        end
    else % nargin >= 3
        model_version = string(arg2);
        if ~isempty(arg3), method = upper(string(arg3)); end
    end

    if strlength(model_version) == 0
        model_version = infer_model_version(params);
    end

    % ---- build canonical key ----
    key = struct();
    key.model_version = string(model_version);
    key.H0_1 = defaultArg(params,'H0_1',NaN);
    key.H0_2 = defaultArg(params,'H0_2',NaN);
    key.A    = defaultArg(params,'A',   NaN);
    key.V    = defaultArg(params,'V',   NaN);
    key.KA   = defaultArg(params,'KA',  NaN);
    key.KB   = defaultArg(params,'KB',  NaN);
    key.KG   = defaultArg(params,'KG',  NaN);

    % ---- stable serialization & hash ----
    hex = simpleDataHash_internal(key, char(method));
end

function mv = infer_model_version(params)
    mv = "BVP-v3.1";
    if isstruct(params)
        if isfield(params,'model_version') && strlength(string(params.model_version))>0
            mv = string(params.model_version); return;
        end
        if isfield(params,'version') && strlength(string(params.version))>0
            mv = string(params.version); return;
        end
    end
end

function hex = simpleDataHash_internal(x, method)
    if nargin < 2 || isempty(method)
        method = 'SHA-256';
    end
    x = normalize_fields(x);
    bytes = getByteStreamFromArray(x);
    md = java.security.MessageDigest.getInstance(method);
    md.update(uint8(bytes));
    digest = typecast(md.digest(), 'uint8');
    hex = lower(reshape(dec2hex(digest)', 1, []));
end

function y = normalize_fields(x)
    if isstruct(x)
        f = sort(fieldnames(x));
        for i = 1:numel(f)
            [x.(f{i})] = normalize_fields(x.(f{i}));
        end
        y = orderfields(x, f);
    elseif iscell(x)
        y = cellfun(@normalize_fields, x, 'UniformOutput', false);
    else
        y = x;
    end
end
