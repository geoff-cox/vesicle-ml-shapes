function [result, meta] = solveAtParams_v1(params, sim, warm)
% SOLVEATPARAMS  Solve the vesicle BVP at given H0 using multi-rung continuation.
% Supports verbose logging via sim.SP.Verbose (true/false)
% and optional homotopy recording via sim.SP.SaveHomotopy (true/false).

    arguments
        params (1,1) struct
        sim    (1,1) struct
        warm   (1,1) struct = struct()
    end

    % ---------- configuration ----------
    TH     = sim.TH;
    delta0 = TH.delta;
    opts0  = TH.opts;

    verbose      = isfield(sim.SP,'Verbose')      && sim.SP.Verbose;
    saveHomotopy = isfield(sim.SP,'SaveHomotopy') && sim.SP.SaveHomotopy;

    MP   = sim.MP;
    Htgt = [params.H0_1, params.H0_2];
    [aS,bS] = computePhaseScales(MP.A);

    basePar = struct('A',MP.A,'V',MP.V,'KA',MP.KA,'KB',MP.KB,'KG',MP.KG, ...
                     'aS',aS,'bS',bS,'delta',delta0);

    say = @(fmt,varargin) if_verbose(verbose,fmt,varargin{:});

    say('\n=== solveAtParams: starting solve at H0 = [%+.3f, %+.3f] ===', Htgt);

    % ---------- initial guess ----------
    initSol = [];
    if isfield(warm,'result') && isstruct(warm.result) && isfield(warm.result,'sol')
        initSol = warm.result.sol;
        say('Warm-start loaded from previous solved point.');
    end
    if isempty(initSol)
        seedParams = struct('A',MP.A,'V',MP.V,'KG',MP.KG,'KA',MP.KA,'KB',MP.KB);
        initSol = initialGuessFromFile(seedParams, Htgt);
        say('Using seed initial shape (InitialShapes folder).');
    end

    if isfield(warm,'fromParams') && all(isfield(warm.fromParams,{'H0_1','H0_2'}))
        H0_from = [warm.fromParams.H0_1, warm.fromParams.H0_2];
    else
        H0_from = Htgt;  % no known predecessor
    end
    Hfrom = H0_from;    % (alias used below)

    % Reseed if we cross an axis between Hfrom and Htgt
    crossesZero = ( (H0_from(1) <= 0 && Htgt(1) >= 0) || (H0_from(1) >= 0 && Htgt(1) <= 0) ) ...
               || ( (H0_from(2) <= 0 && Htgt(2) >= 0) || (H0_from(2) >= 0 && Htgt(2) <= 0) );
    
    if crossesZero
        say('    ↺ re-seed at [0,0] (crossed axis), trying again…');
        seedParams = struct('A',MP.A,'V',MP.V,'KG',MP.KG,'KA',MP.KA,'KB',MP.KB);
        seedSol = initialGuessFromFile(seedParams, [0,0]);  % your InitialShapes
        if ~isempty(seedSol), initSol = seedSol; end
    end

    % ---------- attempt configuration ----------
    stepCaps = [Inf, 0.50, 0.25, 0.12, 0.06, 0.03];  % H0 stepping ladders
    if isfield(TH,'delta_list') && ~isempty(TH.delta_list)
        deltas = TH.delta_list(:)';
    else
        % default: baseline, slightly larger, and smaller
        deltas = [delta0, min(2*delta0, 0.02), max(0.5*delta0, 5e-3)];
    end
    optSets  = {opts0, bvpset(opts0,'RelTol',1e-5,'AbsTol',1e-7)};

    homolog   = [];
    accepted  = false;
    sol       = [];
    usedPar   = basePar;
    
    % allow one axis-path rescue per solve (avoids loops)
    axisFallbackTried = false;

    % ---------- main attempt loop ----------
    rung = 0;
    for s = 1:numel(stepCaps)
        stepMax = stepCaps(s);
        rung = rung + 1;
        say('\n-- Rung %d/%d: stepMax = %.3g --', rung, numel(stepCaps), stepMax);

        if isfinite(stepMax)
            L = norm(Htgt - H0_from);
            nSteps = max(1, ceil(L/stepMax));
            Hpath = H0_from + (1:nSteps).'/nSteps .* (Htgt - H0_from);
            say('  continuation path: %d steps (total ∆H0 = %.3g)', nSteps, L);
        else
            Hpath = Htgt;  % direct attempt
            say('  direct attempt at target H0');
        end

        for dd = 1:numel(deltas)
            deltaNow = deltas(dd);
            for oo = 1:numel(optSets)
                optsNow = optSets{oo};
                say('    δ = %.4g | RelTol = %.1e', deltaNow, bvpget(optsNow,'RelTol'));
                curSol = initSol;

                % Build a segment Hprev -> Htarget with cap stepMax; bisect on failure.
                Hprev   = H0_from;
                Htarget = Htgt;
                if isfinite(stepMax)
                    maxSeg = stepMax;
                else
                    maxSeg = norm(Htarget - Hprev); % direct attempt if Inf
                end
                minSeg = 0.01;  % configurable: sim.TH.minH0Step
                
                while true
                    v = Htarget - Hprev;
                    L = norm(v);
                    if L <= maxSeg
                        next = Htarget;
                    else
                        next = Hprev + v * (maxSeg / L);
                    end
                
                    usedPar = basePar; 
                    usedPar.H0    = next; 
                    usedPar.delta = deltaNow;

                    odefun = @(s_,y_,lam_) BendV_Lag_EIGp_DE_impl(s_,y_,lam_,usedPar);
                    bcfun  = @(ya,yb,lam_) BendV_Lag_EIGp_BC_impl(ya,yb,lam_,usedPar);
                
                    try
                        curSol = bvp6c(odefun, bcfun, curSol, optsNow);
                        % success: advance

                        if saveHomotopy
                            [BCi,~]   = bc_diagnostics(curSol, bcfun);
                            [DEi,~,~] = de_residual(curSol, odefun);
                            homolog(end+1) = struct( ...
                                'H0'    , H0, ...
                                'delta' , usedPar.delta, ...
                                'mesh'  , numel(curSol.x), ...
                                'BCmax' , BCi, ...
                                'DEmax' , DEi); %#ok<AGROW>
                        end

                        if all(next == Htarget)
                            break; % reached this rung's target
                        end
                        Hprev = next;
                    catch
                        temp = maxSeg;
                        % failure: halve the segment length, retry unless we're below minSeg
                        maxSeg = max(maxSeg/2, minSeg);
                        say('    ✗ step failed: singular Jacobian | halving step %.4g → %.4g (δ=%.4g, RT=%.1e)', temp, maxSeg, deltaNow, bvpget(optsNow,'RelTol'));
                        if maxSeg <= minSeg + eps
                            curSol = [];    % give up on this δ/opts combo
                            break
                        end
                    end
                end

                % verify final target
                if ~isempty(curSol)
                    bcfunT  = @(ya,yb,lam_) BendV_Lag_EIGp_BC_impl(ya,yb,lam_,usedPar);
                    odefunT = @(s_,y_,lam_) BendV_Lag_EIGp_DE_impl(s_,y_,lam_,usedPar);
                    [BCmax,~]   = bc_diagnostics(curSol, bcfunT);
                    [DEmax,~,~] = de_residual(curSol, odefunT);
                    rmin        = local_min_radius_interior(curSol);
                    say('      Final gates: BC=%.2e | DE=%.2e | rmin=%.2e', BCmax,DEmax,rmin);

                    if (BCmax <= TH.BCmax) && (DEmax <= TH.DEmaxHard) && (rmin >= TH.rMin)
                        sol = curSol;
                        accepted = true;
                        say('      ✅ accepted at rung %d, δ=%.3g', rung, deltaNow);
                        break
                    else
                        say('      rejected by gates; BC=%.2e (max %.2e) | DE=%.2e (max %.2e) | rmin=%.2e (min %.2e).', ...
                            BCmax,TH.BCmax,DEmax,TH.DEmaxHard,rmin,TH.rMin);
                        initSol = curSol; % warm the next attempt
                    end
                end

                if accepted, break; end
            end
            if accepted, break; end
        end
        if accepted, break; end
    end

    if ~accepted
        error('solveAtParams: failed to converge after all rungs');
    end

    % ---------- label & metrics ----------
    [label, E_total, P_osm] = labelFromSolution(sol);
    bcfunT  = @(ya,yb,lam_) BendV_Lag_EIGp_BC_impl(ya,yb,lam_,usedPar);
    odefunT = @(s_,y_,lam_) BendV_Lag_EIGp_DE_impl(s_,y_,lam_,usedPar);
    [BCmax,~]   = bc_diagnostics(sol, bcfunT);
    [DEmax,~,~] = de_residual(sol, odefunT);

    result = struct('sol', sol, 'mesh', numel(sol.x));
    meta   = struct('label',label,'E',E_total,'P',P_osm, ...
                    'BCmax',BCmax,'DEmax',DEmax,'mesh',result.mesh);

    if saveHomotopy
        meta.homotopy = homolog;
    end

    say('=== solveAtParams: completed at H0=[%+.3f,%+.3f] ===\n', Htgt);
end

% ---------- local helpers ----------
function if_verbose(flag,fmt,varargin)
    if flag
        fprintf([fmt '\n'], varargin{:});
    end
end