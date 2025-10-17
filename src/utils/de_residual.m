% -------------------------------------------------------------------------
% EXTRACTED HELPER for "de_residual"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:12
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function [rMax, rComp, worstIdx] = de_residual(sol, odefun)
    % Second-order central difference on nonuniform grid, ignoring pole buffers.

    lam = sol.parameters(:);
    s   = sol.x;             % s in [0, pi]
    Y   = sol.y;
    n   = numel(s);
    if n < 3, error('de_residual: need at least 3 mesh points'); end

    % interior nodes (2..n-1)
    h   = diff(s);
    dY  = zeros(size(Y,1), n-2);
    for j = 2:n-1
        h0 = s(j)   - s(j-1);
        h1 = s(j+1) - s(j);
        dY(:,j-1) = ( -(h1/(h0*(h0+h1)))*Y(:,j-1) ...
                      + ((h1-h0)/(h0*h1))*Y(:,j) ...
                      + (h0/(h1*(h0+h1)))*Y(:,j+1) );
    end

    sC = s(2:end-1);
    YC = Y(:,2:end-1);
    
    % --- mask out pole buffers ---
    % estimate delta from first and last spacing if not known here
    % but better: pass it; we can infer from odefun only with extra plumbing.
    % We'll be conservative and drop 2*mean(h) near each pole as a buffer:
    pad = 2*mean(h);
    mask = (sC > pad) & (sC < (pi - pad));
    if ~any(mask)
        % fallback: evaluate everywhere
        mask = true(size(sC));
    end
    sC = sC(mask);
    YC = YC(:,mask);
    dY = dY(:,mask);

    F = zeros(size(YC));
    for j = 1:numel(sC)
        F(:,j) = odefun(sC(j), YC(:,j), lam);
    end

    R = dY - F;
    if isempty(R)
        rComp = zeros(size(Y,1),1);
    else
        rComp = max(abs(R),[],2);
    end
    [rMax, worstIdx] = max(rComp);
end