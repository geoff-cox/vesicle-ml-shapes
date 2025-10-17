% -------------------------------------------------------------------------
% EXTRACTED HELPER for "continuationSolve"
%   - Source: sim_driver_quad_tree_full.m
%   - Extracted: 2025-10-11 11:50:11
%   - Sub-helpers (nested functions) are retained in this file.
% -------------------------------------------------------------------------

function sol = continuationSolve(S, p0, sol0, pT, stepMax)
    % Two-leg path: (p0_1,p0_2) -> (pT_1,p0_2) -> (pT_1,pT_2)
    path = [linspace(p0(1), pT(1), max(2,ceil(abs(pT(1)-p0(1))/stepMax)))'  repmat(p0(2),[],1)];
    path = [path;  [repmat(pT(1), max(2,ceil(abs(pT(2)-p0(2))/stepMax)))'  linspace(p0(2), pT(2), max(2,ceil(abs(pT(2)-p0(2))/stepMax)))']];

    sol = sol0;
    Par = S.params;
    baseDelta = S.solver.delta;

    for k = 2:size(path,1)
        target = path(k,:);
        % adaptive attempt: try up to 5 backtracks with smaller step and larger delta
        ok = false; localStep = norm(target - path(k-1,:));
        for attempt = 1:5
            % build handles with a slightly larger delta near failures
            bump = min(3, 1 + 0.5*(attempt-1));     % 1, 1.5, 2, 2.5, 3 x
            Par.H0    = target;
            Par.delta = max(1e-3, min(0.03, baseDelta*bump));

            odefun = @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par);
            bcfun  = @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par);

            try
                sol = bvp6c(odefun, bcfun, sol, S.solver.opts);
                [BCmax,BCrep] = bc_diagnostics(sol,bcfun);
                if BCmax <= 1e-6
                    ok = true;
                    break;
                else
                    % reject and backtrack
                    target = path(k-1,:) + 0.5*(target - path(k-1,:)); % halve step
                end
            catch
                % backtrack on failure
                target = path(k-1,:) + 0.5*(target - path(k-1,:));
            end
        end
        if ~ok
            error('continuationSolve: could not reach waypoint (%.3f,%.3f)', target(1),target(2));
        end
        if S.plot.enable
            plot(S.plot.ax, Par.H0(1), Par.H0(2), '.', 'MarkerSize',6); drawnow limitrate;
        end
    end
end