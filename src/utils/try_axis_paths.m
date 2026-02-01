function [ok, initSol, fromH0] = try_axis_paths(warmSol, Hfrom, Htgt, sim, stepMax)
% TRY_AXIS_PATHS  Axis-aligned fallback when 2D continuation stalls.
% Try two one-dimensional homotopies: (H0_1 then H0_2) and (H0_2 then H0_1).
% Returns first that completes at least one accepted step (ok=true) along either axis.

    ok = false; initSol = []; fromH0 = Hfrom;

    % Path A: H0_1 -> H0_2
    [okA, solA] = oned_homotopy(warmSol, Hfrom, [Htgt(1) Hfrom(2)], sim, stepMax);
    if okA
        [okB, solB] = oned_homotopy(solA, [Htgt(1) Hfrom(2)], Htgt, sim, stepMax);
        if okB, ok = true; initSol = solB; fromH0 = Htgt; return; end
    end

    % Path B: H0_2 -> H0_1
    [okB1, solB1] = oned_homotopy(warmSol, Hfrom, [Hfrom(1) Htgt(2)], sim, stepMax);
    if okB1
        [okB2, solB2] = oned_homotopy(solB1, [Hfrom(1) Htgt(2)], Htgt, sim, stepMax);
        if okB2, ok = true; initSol = solB2; fromH0 = Htgt; return; end
    end
end

function [ok, solOut] = oned_homotopy(solIn, p0, p1, sim, stepMax)
% March in fixed-size steps along one axis; accept first step that passes gates.
    ok = false; solOut = solIn;
    if all(p0==p1), ok = true; return; end
    nSteps = max(1, ceil(norm(p1-p0)/stepMax));
    for k=1:nSteps
        alpha = k/nSteps;
        Hk = (1-alpha)*p0 + alpha*p1;
        [accept, solK] = try_one_shot(solOut, Hk, sim);
        if accept
            ok = true; solOut = solK; return;  % return as soon as we get one accepted step
        else
            % If Newton failed hard, coarsen the mesh and retry once at same Hk
            [accept2, solK2] = try_one_shot(coarsen_mesh(solOut, 0.5), Hk, sim);
            if accept2
                ok = true; solOut = solK2; return;
            end
        end
    end
end

function [accept, sol] = try_one_shot(solInit, H0, sim)
% Single attempt at a given H0 using the current ladder's deltas/tols.
    accept = false; sol = [];
    TH = sim.TH; MP = sim.MP;
    [aS,bS] = computePhaseScales(MP.A);
    Par0 = struct('H0',H0,'A',MP.A,'V',MP.V,'KA',MP.KA,'KB',MP.KB,'KG',MP.KG, ...
                  'aS',aS,'bS',bS,'delta',sim.TH.delta);

    deltas = sim.TH.delta_list;  % e.g., [0.01, 0.02, 0.005]
    opts   = {sim.TH.opts, bvpset(sim.TH.opts,'RelTol',1e-5,'AbsTol',1e-7)};

    for d = 1:numel(deltas)
      for o = 1:numel(opts)
          Par = Par0; Par.delta = deltas(d);
          odefun = @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par);
          bcfun  = @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par);
          try
              sol1 = bvp6c(odefun, bcfun, solInit, opts{o});
              [BCmax,~] = bc_diagnostics(sol1, bcfun);
              [DEmax,~,~] = de_residual(sol1, odefun);
              rmin = local_min_radius_interior(sol1);
              if (BCmax <= TH.BCmax) && (DEmax <= TH.DEmaxHard) && (rmin >= TH.rMin)
                  accept = true; sol = sol1; return;
              else
                  solInit = sol1; % warm next try
              end
          catch
              % continue
          end
      end
    end
end
