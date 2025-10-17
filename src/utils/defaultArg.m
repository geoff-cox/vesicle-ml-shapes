% -------------------------------------------------------------------------
% EXTRACTED HELPER for "defaultArg"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function v = defaultArg(s, field, vDefault)
    if isfield(s, field), v = s.(field); else, v = vDefault; end
end