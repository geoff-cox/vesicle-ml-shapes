function h = param_hash_v2(model_version, P)
%PARAM_HASH_V2 Deterministic hash for equation-defining parameters.
% Includes exactly:
%   model_version, H0_1, H0_2, A, V, KA, KB, KG
%
% Usage:
%   h = param_hash_v2(sim.SP.ModelVersion, paramsStruct);

    mv = string(model_version);

    key = join([
        "model_version=" + mv
        "H0_1=" + f(P.H0_1)
        "H0_2=" + f(P.H0_2)
        "A="    + f(P.A)
        "V="    + f(P.V)
        "KA="   + f(P.KA)
        "KB="   + f(P.KB)
        "KG="   + f(P.KG)
    ], ";");

    h = hash_sha256(key);
end

function s = f(x)
    if isempty(x) || ~isscalar(x) || ~isfinite(x)
        error('param_hash_v2:badparam','All hashed params must be finite scalars.');
    end
    s = string(sprintf('%.17g', double(x)));
end
