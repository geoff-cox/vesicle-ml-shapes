function warm = pickWarmStart(params, sim, simDir, T)
% PICKWARMSTART  Choose a warm-start for (H0_1,H0_2) with fixed physics in sim.MP.
% Priority:
%   1) Nearest already-solved point (matching physics) from catalog -> hashed_results/<hash>.mat
%   2) Seed file from InitialShapes/ with matching physics (A,V,KG,KA,KB)
% Returns:
%   warm.result.sol  OR  warm.sol   (either is accepted by solveAtParams)
%
% Usage:
%   warm = pickWarmStart(task.params, sim, simDir, catalog_load(simDir))

    warm = struct();
    if nargin < 4 || isempty(T), T = catalog_load(simDir); end

    % ---------- physics to match (use small tolerances for floats) ----------
    MP = sim.MP;
    tol.A  = 1e-9;
    tol.V  = 1e-9;
    % KA/KB/KG are typically integers/small; equality is fine:
    eq = @(a,b,eps) (isnumeric(a) && isnumeric(b) && abs(a-b) <= eps) || isequal(a,b);

    % ---------- 1) nearest solved neighbor (same physics) ----------
    if any(T.Properties.VariableNames == "entry") && ~isempty(T)
        H1 = cellfun(@(e) getp(e,'H0_1'), T.entry);
        H2 = cellfun(@(e) getp(e,'H0_2'), T.entry);
        Av = cellfun(@(e) getp(e,'A'),    T.entry);
        Vv = cellfun(@(e) getp(e,'V'),    T.entry);
        KAv= cellfun(@(e) getp(e,'KA'),   T.entry);
        KBv= cellfun(@(e) getp(e,'KB'),   T.entry);
        KGv= cellfun(@(e) getp(e,'KG'),   T.entry);

        samePhys = arrayfun(@(i) ...
            eq(Av(i),MP.A,tol.A) && eq(Vv(i),MP.V,tol.V) && KAv(i)==MP.KA && KBv(i)==MP.KB && KGv(i)==MP.KG, ...
            (1:numel(Av))');

        % Keep rows that have meta from a SOLVE (seeds have meta.type="seed")
        isSolve = false(height(T),1);
        for i=1:height(T)
            m = T.entry{i}.meta;
            isSolve(i) = ~(isstruct(m) && isfield(m,'type') && strcmpi(string(m.type),'seed'));
        end

        mask = samePhys & isSolve & isfinite(H1) & isfinite(H2);
        if any(mask)
            H  = [params.H0_1 params.H0_2];
            d2 = (H1(mask)-H(1)).^2 + (H2(mask)-H(2)).^2;
            [~,ixRel] = min(d2);
            ix = find(mask);
            ix = ix(ixRel);

            % load the hashed result as warm-start
            meta = T.entry{ix}.meta;
            if isfield(meta,'hash')
                f = fullfile(simDir,'hashed_results', meta.hash + ".mat");
                if exist(f,'file')==2
                    S = load(f,'result','sol'); % support either field name
                    if isfield(S,'result') && isstruct(S.result) && isfield(S.result,'sol')
                        warm.result    = struct('sol', S.result.sol);
                        % record the neighbor's parameters for continuation
                        warm.fromParams = struct('H0_1', H1(ix), 'H0_2', H2(ix));
                        return
                    elseif isfield(S,'sol') && isstruct(S.sol)
                        warm.result    = struct('sol', S.sol);
                        warm.fromParams = struct('H0_1', H1(ix), 'H0_2', H2(ix));
                        return
                    end
                end
            end
        end

    end

    % ---------- 2) fall back to seed file from InitialShapes/ ----------
    % We stored seeds in catalog via bootstrap() as entry.meta.type="seed" and meta.source=filename
    if any(T.Properties.VariableNames == "entry") && ~isempty(T)
        % physics match only; seeds do not include H0
        isSeed = false(height(T),1);
        for i=1:height(T)
            m = T.entry{i}.meta;
            isSeed(i) = isstruct(m) && isfield(m,'type') && strcmpi(string(m.type),'seed');
        end

        Av = cellfun(@(e) getp(e,'A'),  T.entry);
        Vv = cellfun(@(e) getp(e,'V'),  T.entry);
        KAv= cellfun(@(e) getp(e,'KA'), T.entry);
        KBv= cellfun(@(e) getp(e,'KB'), T.entry);
        KGv= cellfun(@(e) getp(e,'KG'), T.entry);

        samePhys = arrayfun(@(i) ...
            eq(Av(i),MP.A,tol.A) && eq(Vv(i),MP.V,tol.V) && KAv(i)==MP.KA && KBv(i)==MP.KB && KGv(i)==MP.KG, ...
            (1:numel(Av))');

        mask = isSeed & samePhys;
        if any(mask)
            ix = find(mask,1,'first');  % any matching seed is fine; they're all H0=[0,0]
            src = "";
            if isfield(T.entry{ix}.meta,'source')
                src = string(T.entry{ix}.meta.source);
            end
            ishapesDir = fullfile(fileparts(simDir), 'InitialShapes'); % project_root/InitialShapes
            f = src;
            if ~isempty(src) && exist(f,'file')~=2
                f = fullfile(ishapesDir, src);
            end
            if exist(f,'file')==2
                tmp = load(f);
                if isfield(tmp,'Version') && numel(tmp.Version)>=1 && isfield(tmp.Version(1),'Solution')
                    warm.sol = tmp.Version(1).Solution;
                    return
                end
            end
        end
    end

    % ---------- no warm-start available ----------
    % leave 'warm' empty; solveAtParams will fall back to initialGuessFromFile if needed
end

% ---- local helpers ----
function v = getp(e, k)
    if isstruct(e) && isfield(e,'params') && isstruct(e.params) && isfield(e.params, k)
        v = double(e.params.(k));
    else
        v = NaN;
    end
end
