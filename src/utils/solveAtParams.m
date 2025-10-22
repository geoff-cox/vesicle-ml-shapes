function [result, meta] = solveAtParams(params, sim, warm)
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

    verbose = isfield(sim.SP,'Verbose') && sim.SP.Verbose;
    saveHomotopy = isfield(sim.SP,'SaveHomotopy') && sim.SP.SaveHomotopy;

    MP   = sim.MP;
    Htgt = [params.H0_1, params.H0_2];
    [aS,bS] = computePhaseScales(MP.A);

    basePar = struct('A',MP.A,'V',MP.V,'KA',MP.KA,'KB',MP.KB,'KG',MP.KG, ...
                     'aS',aS,'bS',bS,'delta',delta0);

    say = @(fmt,varargin) if_verbose(verbose,fmt,varargin{:}); %#ok<NASGU>

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
        H0_from = Htgt;
    end

    % Before starting a segment Hprev -> Htarget:
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
    deltas = unique([ ...
        max(0.25*delta0, 5e-3), ...
        max(0.5*delta0 , 5e-3), ...
        delta0, ...
        min(2*delta0,   5e-2), ...
        min(3*delta0,   8e-2) ...
    ]);
    optSets  = {opts0, bvpset(opts0,'RelTol',1e-5,'AbsTol',1e-7)};

    homolog = [];
    accepted = false; sol = []; usedPar = basePar;

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

                % Build a segment from Hprev -> HtargetStep, try, and on failure bisect until min step reached.
                Hprev = H0_from;
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
                
                    usedPar = basePar; usedPar.H0 = next; usedPar.delta = deltaNow;
                    odefun = @(s_,y_,lam_) BendV_Lag_EIGp_DE_impl(s_,y_,lam_,usedPar);
                    bcfun  = @(ya,yb,lam_) BendV_Lag_EIGp_BC_impl(ya,yb,lam_,usedPar);
                
                    try
                        curSol = bvp6c(odefun, bcfun, curSol, optsNow);
                        % success: advance

                        if saveHomotopy
                            [BCi,~]   = bc_diagnostics(curSol, bcfun);
                            [DEi,~,~] = de_residual(curSol, odefun);
                            homItem = struct('H0',H0,'delta',usedPar.delta, ...
                                             'mesh',numel(curSol.x),'BCmax',BCi,'DEmax',DEi);
                            homolog = [homolog; homItem]; %#ok<AGROW>
                        end

                        if all(next == Htarget), break; end   % reached this rung's target
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






                % for i = 1:size(Hpath,1)
                %     H0 = Hpath(i,:);
                % 
                %     usedPar = basePar; usedPar.H0 = H0; usedPar.delta = deltaNow;
                %     odefun = @(s_,y_,lam_) BendV_Lag_EIGp_DE_impl(s_,y_,lam_,usedPar);
                %     bcfun  = @(ya,yb,lam_) BendV_Lag_EIGp_BC_impl(ya,yb,lam_,usedPar);
                % 
                %     try
                %         curSol = bvp6c(odefun, bcfun, curSol, optsNow);
                %         say('      ✓ step %d/%d converged @ H0=[%+.3f,%+.3f], mesh=%d', ...
                %             i, size(Hpath,1), H0(1), H0(2), numel(curSol.x));
                %     catch ME
                %         say('      ✗ step %d/%d failed: %s', i, size(Hpath,1), ME.message);
                %         curSol = []; break;
                %     end
                % 
                %     if saveHomotopy
                %         [BCi,~]   = bc_diagnostics(curSol, bcfun);
                %         [DEi,~,~] = de_residual(curSol, odefun);
                %         homItem = struct('H0',H0,'delta',usedPar.delta, ...
                %                          'mesh',numel(curSol.x),'BCmax',BCi,'DEmax',DEi);
                %         homolog = [homolog; homItem]; %#ok<AGROW>
                %     end
                % end

                % verify final target
                if ~isempty(curSol)
                    odefunT = @(s_,y_,lam_) BendV_Lag_EIGp_DE_impl(s_,y_,lam_,usedPar);
                    bcfunT  = @(ya,yb,lam_) BendV_Lag_EIGp_BC_impl(ya,yb,lam_,usedPar);
                    [BCmax,~]   = bc_diagnostics(curSol, bcfunT);
                    [DEmax,~,~] = de_residual(curSol, odefunT);
                    rmin        = local_min_radius_interior(curSol);
                    say('      Final gates: BC=%.2e | DE=%.2e | rmin=%.2e', BCmax,DEmax,rmin);

                    if (BCmax <= TH.BCmax) && (DEmax <= TH.DEmaxHard) && (rmin >= TH.rMin)
                        sol = curSol; accepted = true;
                        say('      ✅ accepted at rung %d, δ=%.3g', rung, deltaNow);
                        goto_done = true; %#ok<NASGU>
                    else
                        say('      rejected by gates; BC=%.2e (max %.2e) | DE=%.2e (max %.2e) | rmin=%.2e (min %.2e).', BCmax,TH.BCmax,DEmax,TH.DEmaxHard,rmin,TH.rMin);
                        initSol = curSol;
                    end
                end

                if exist('goto_done','var'), break; end
            end
            if exist('goto_done','var'), break; end
        end
        if exist('goto_done','var'), break; end
    end

    if ~accepted
        error('solveAtParams: failed to converge after all rungs');
    end

    % ---------- label & metrics ----------
    [label, E_total, P_osm] = labelFromSolution(sol);
    odefunT = @(s_,y_,lam_) BendV_Lag_EIGp_DE_impl(s_,y_,lam_,usedPar);
    bcfunT  = @(ya,yb,lam_) BendV_Lag_EIGp_BC_impl(ya,yb,lam_,usedPar);
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

% ---------- local helper ----------
function if_verbose(flag,fmt,varargin)
    if flag
        fprintf([fmt '\n'], varargin{:});
    end
end


% 
% 
% 
% function [result, meta] = solveAtParams(params, sim, warm)
% % SOLVEATPARAMS  Compute a vesicle solution at (H0_1,H0_2) with robust gates.
% % Contract (refactor):
% %   INPUT:
% %     params : struct with fields at least H0_1, H0_2.
% %              Optional physics knobs: A, V, KA, KB, KG, etc.
% %              Optional solver knobs:  delta, opts (bvpset struct)
% %     warm   : optional struct for warm-starts. Recognized fields:
% %              - warm.result.sol  or warm.sol   : a prior bvpinit/bvp6c solution struct
% %              - warm.meta        : prior meta (ignored here except maybe for seeds)
% %   OUTPUT:
% %     result : struct with fields {s,y,bcResidual,deResidual} OR the raw bvp6c 'sol'
% %              (below we return the raw bvp6c 'sol' as 'result.sol' and fill quick diags)
% %     meta   : struct with quick scalars for catalog:
% %              label (string), E (double), P (double),
% %              BCmax (double), DEmax (double), mesh (double)
% 
%     arguments
%         params (1,1) struct
%         sim    (1,1) struct
%         warm   (1,1) struct = struct()
%     end
% 
%     % --- thresholds & numerical knobs ---
%     TH     = sim.TH;
%     delta0 = TH.delta;
%     opts0  = TH.opts;
% 
%     % --- physics bundle (fixed for this run) ---
%     MP = sim.MP;
%     H0 = [params.H0_1, params.H0_2];
%     [aS,bS] = computePhaseScales(MP.A);
%     Par = struct('H0',H0, 'A',MP.A,'V',MP.V,'KA',MP.KA,'KB',MP.KB,'KG',MP.KG, ...
%                  'aS',aS,'bS',bS, 'delta',delta0);
% 
%     % --- initial guess / warm start ---
%     initSol = [];
%     if isfield(warm,'result') && isstruct(warm.result) && isfield(warm.result,'sol')
%         initSol = warm.result.sol;
%     elseif isfield(warm,'sol') && isstruct(warm.sol)
%         initSol = warm.sol;
%     end
%     if isempty(initSol)
%         % fallback seed from InitialShapes matching physics
%         seedParams = struct('A',MP.A,'V',MP.V,'KG',MP.KG,'KA',MP.KA,'KB',MP.KB);
%         initSol = initialGuessFromFile(seedParams, H0);
%     end
% 
%     % --- (optional) continuation from a solved H0 toward target H0 ---
%     canContinue = isfield(warm,'fromParams') && all(isfield(warm.fromParams,{'H0_1','H0_2'}));
%     if canContinue
%         Hfrom = [warm.fromParams.H0_1, warm.fromParams.H0_2];
%         Htgt  = [params.H0_1, params.H0_2];
% 
%         fprintf('solveAtParams: 0.15 step continuation from H0 = %s to H0 = %s\n', mat2str(Hfrom), mat2str(H0));
% 
%         % Only continue if we actually have a SOLVED state as initSol
%         hasSolvedInit = (isfield(warm,'result') && isstruct(warm.result) && isfield(warm.result,'sol')) ...
%                      || (isfield(warm,'sol')     && isstruct(warm.sol));
%         if hasSolvedInit && norm(Htgt - Hfrom) > 0.15
%             try
%                 if isfield(warm,'result') && isfield(warm.result,'sol')
%                     initSol = continuation_towards_H0(warm.result.sol, Hfrom, Htgt, sim, 0.15);
%                     fprintf('solveAtParams: continuation solution accepted\n');
%                 end
%             catch
%                 fprintf('solveAtParams: continuation errored - falling back\n');
%                 % best-effort; fall back to the provided warm start without continuation
%                 if isfield(warm,'result'), initSol = warm.result.sol; else, initSol = warm.sol; end
%             end
%         else
%             fprintf('continuation condition "hasSolvedInit && norm(Htgt - Hfrom) > 0.15" is false\n');
%             % distance small or no proven solved source -> skip continuation
%             if isempty(initSol)
%                 if isfield(warm,'result'), initSol = warm.result.sol; else, initSol = warm.sol; end
%             end
%         end
%     end
% 
%     % --- attempt ladder (delta / tolerances) ---
%     attempts = [ ...
%       struct('delta', delta0,                        'opts', opts0); ...
%       struct('delta', max(0.5*delta0, 1e-3),         'opts', opts0); ...
%       struct('delta', max(0.25*delta0, 1e-3),        'opts', opts0); ...
%       struct('delta', delta0, 'opts', bvpset(opts0,'RelTol',1e-5,'AbsTol',1e-7)) ...
%     ];
% 
%     % suppress bvp6c diagnostic warnings
%     warnState = warning('off','MATLAB:bvp6c:RelTolNotMet');
% 
%     solved = false; sol = [];
%     for a = 1:numel(attempts)
%         Par.delta = attempts(a).delta;          % numerical, not in hash
%         Par.H0    = H0;                         % ensure H0 is current target
%         odefun = @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par);
%         bcfun  = @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par);
%         fprintf('solveAtParams: attempting to solve at %s with δ = %g\n',mat2str(H0),Par.delta)
%         try
%             sol = bvp6c(odefun, bcfun, initSol, attempts(a).opts);
% 
%             % gates
%             [BCmax,~]       = bc_diagnostics(sol, bcfun);
%             [DEmax,~,~]     = de_residual(sol, odefun);
%             rmin            = local_min_radius_interior(sol);
% 
%             if (BCmax <= TH.BCmax) && (DEmax <= TH.DEmaxHard) && (rmin >= TH.rMin)
%                 fprintf('solveAtParams: solution accepted gates passed\n')
%                 solved = true; break;
%             else
%                 fprintf('solveAtParams: solution rejected by gates\n')
%                 initSol = sol;  % warm next rung from current iterate
%             end
%         catch
%             fprintf('solveAtParams: solution rejected due to error\n')
%             % try next rung
%         end
%     end
%     % restore original warning state after success
%     warning(warnState);
% 
%     if ~solved, error('solveAtParams: failed to converge'); end
% 
%     % --- label & metrics ---
%     [label, E_total, P_osm] = labelFromSolution(sol);
%     [BCmax,~] = bc_diagnostics(sol, @(ya,yb,lam) BendV_Lag_EIGp_BC_impl(ya,yb,lam,Par));
%     [DEmax,~,~] = de_residual(sol, @(s,y,lam) BendV_Lag_EIGp_DE_impl(s,y,lam,Par));
% 
%     result = struct('sol', sol, 'mesh', numel(sol.x));
%     meta   = struct('label',label,'E',E_total,'P',P_osm, ...
%                     'BCmax',BCmax,'DEmax',DEmax,'mesh',result.mesh);
% end
