% -------------------------------------------------------------------------
% EXTRACTED HELPER for "simpleDataHash"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function hex = simpleDataHash(x, method)
    if nargin < 2, method = 'SHA-256'; end
    x = normalizeForHash(x);
    bytes  = getByteStreamFromArray(x);
    md     = java.security.MessageDigest.getInstance(method);
    md.update(typecast(uint8(bytes),'int8'));
    digest = typecast(md.digest(),'uint8');
    hex    = lower(reshape(dec2hex(digest)',1,[]));
end