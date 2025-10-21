% -------------------------------------------------------------------------
% EXTRACTED HELPER for "nearest_by_sign"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function idx = nearest_by_sign(cacheKeys, target)
    tgtSign = sign(target);
    same = all(sign(cacheKeys) == tgtSign, 2);
    pool = find(same);
    if isempty(pool), pool = 1:size(cacheKeys,1); end
    [~,j] = min(sum((cacheKeys(pool,:) - target).^2,2));
    idx = pool(j);
end