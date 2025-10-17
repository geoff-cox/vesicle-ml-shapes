% -------------------------------------------------------------------------
% EXTRACTED HELPER for "labelFromSolution"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function [label, E_total, P_osm] = labelFromSolution(sol)
    E_total = sol.y(9,end) - sol.y(18,end) + 0.5*sol.y(4,end);
    P_osm   = sol.parameters(1);
    rA = sol.y(4,:); rB = sol.y(13,:);
    rr = [rA rB];

    if exist('findpeaks','file')
        np = numel(findpeaks(rr));
    else
        % very light peak count without toolbox
        np = sum( rr(2:end-1) > rr(1:end-2) & rr(2:end-1) >= rr(3:end) );
    end

    if     np <= 2, label = 1;
    elseif np <= 4, label = 2;
    else            label = 3;
    end
end