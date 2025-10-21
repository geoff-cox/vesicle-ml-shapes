% -------------------------------------------------------------------------
% EXTRACTED HELPER for "correctToBoundary"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function [pNext, nOut, success, cache] = correctToBoundary(S, pPred, nIn, pad, epsBis, cache)
    success = false;
    n = unit(nIn);

    pL = pPred - pad*n;  pR = pPred + pad*n;
    if ~inBox(pL, [S.grid.A_lim; S.grid.B_lim]) || ~inBox(pR, [S.grid.A_lim; S.grid.B_lim])
        pNext = pPred; nOut = n; return;
    end

    [lL, cache] = getLabel(S, pL, cache);
    [lR, cache] = getLabel(S, pR, cache);

    if lL == lR
        n = -n;
        pL = pPred - pad*n;  pR = pPred + pad*n;
        [lL, cache] = getLabel(S, pL, cache);
        [lR, cache] = getLabel(S, pR, cache);
        if lL == lR
            pad2 = 2*pad;
            pL = pPred - pad2*n; pR = pPred + pad2*n;
            [lL, cache] = getLabel(S, pL, cache);
            [lR, cache] = getLabel(S, pR, cache);
            if lL == lR, pNext = pPred; nOut = n; return; end
        end
    end

    a = pL; b = pR; la = lL; lb = lR;
    for it = 1:40
        m = 0.5*(a+b);
        [lm, cache] = getLabel(S, m, cache);
        if lm == la, a = m; else, b = m; end
        if norm(a-b) < epsBis, break; end
    end
    pNext = 0.5*(a+b);
    nOut  = unit(b - a);
    success = true;
end