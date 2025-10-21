% -------------------------------------------------------------------------
% EXTRACTED HELPER for "mergeCaches"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function C = mergeCaches(C, D)
    if isempty(D) || isempty(D.keys), return; end
    for i=1:size(D.keys,1)
        k = D.keys(i,:);
        if isempty(C.keys) || ~ismembertol(k, C.keys, 1e-12, 'ByRows',true)
            C.keys   = [C.keys; k];
            C.labels = [C.labels; D.labels(i)];
            C.E      = [C.E; D.E(i)];
            C.P      = [C.P; D.P(i)];
            C.sols{end+1} = D.sols{i};
        end
    end
end