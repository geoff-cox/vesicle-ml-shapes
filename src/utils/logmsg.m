% -------------------------------------------------------------------------
% EXTRACTED HELPER for "logmsg"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function logmsg(S, fmt, varargin)
    msg = sprintf(fmt, varargin{:});
    if S.log.fid > 0
        fprintf(S.log.fid, '%s\n', msg);
    end
    if S.print2cmdwin
        fprintf('%s\n', msg);
    end
end