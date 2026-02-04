function [j, sha] = canon_json(x)
%CANON_JSON Canonical JSON string with stable field ordering.
% Returns:
%   j   - JSON (1x1 string)
%   sha - SHA-256 hex of UTF-8 JSON bytes (1x1 string)

    y = canon_normalize(x);
    j = string(jsonencode(y));  % compact; stable because fields ordered in y
    sha = hash_sha256(j);
end
