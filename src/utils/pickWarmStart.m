function warm = pickWarmStart(params, sim, simDir, T)
% PICKWARMSTART  Choose a warm-start for (H0_1,H0_2) with fixed physics in sim.MP.
% Priority:
%   1) Nearest already-solved point (matching physics) from catalog -> hashed_results/<hash>.mat
%   2) Seed file from InitialShapes/ with matching physics (A,V,KG,KA,KB)
% Modern-only: expects hashed_results/<hash>.mat with 'result' struct.

    warm = struct();
    if nargin < 4 || isempty(T), T = catalog_load(simDir); end
    MP = sim.MP; Htgt = [params.H0_1 params.H0_2];

    % 1) nearest solved neighbor with identical physics
    if ~isempty(T) && any(T.Properties.VariableNames=="entry")
        H1 = cellfun(@(e) fget(e,'params','H0_1'), T.entry);
        H2 = cellfun(@(e) fget(e,'params','H0_2'), T.entry);
        Av = cellfun(@(e) fget(e,'params','A'),   T.entry);
        Vv = cellfun(@(e) fget(e,'params','V'),   T.entry);
        KAv= cellfun(@(e) fget(e,'params','KA'),  T.entry);
        KBv= cellfun(@(e) fget(e,'params','KB'),  T.entry);
        KGv= cellfun(@(e) fget(e,'params','KG'),  T.entry);

        isSolve = cellfun(@(e) ~isstruct(e.meta) || ~isfield(e.meta,'type') ...
                                   || ~strcmpi(string(e.meta.type),'seed'), T.entry);

        mask = isSolve & Av==MP.A & Vv==MP.V & KAv==MP.KA & KBv==MP.KB & KGv==MP.KG ...
                     & isfinite(H1) & isfinite(H2);
        if any(mask)
            d2 = (H1(mask)-Htgt(1)).^2 + (H2(mask)-Htgt(2)).^2;
            [~,krel] = min(d2); idx = find(mask); ix = idx(krel);
            if isfield(T.entry{ix}.meta,'hash')
                f = fullfile(simDir,'hashed_results', T.entry{ix}.meta.hash + ".mat");
                if exist(f,'file')==2
                    S = load(f,'result');  % modern only
                    if isfield(S,'result') && isfield(S.result,'sol')
                        warm.result = S.result;
                        warm.fromParams = struct('H0_1', H1(ix), 'H0_2', H2(ix));
                        return
                    end
                end
            end
        end
    end

    % 2) seeds from InitialShapes (registered by bootstrap with meta.type="seed")
    if ~isempty(T) && any(T.Properties.VariableNames=="entry")
        isSeed = cellfun(@(e) isstruct(e.meta) && isfield(e.meta,'type') ...
                                 && strcmpi(string(e.meta.type),'seed'), T.entry);
        Av = cellfun(@(e) fget(e,'params','A'),   T.entry);
        Vv = cellfun(@(e) fget(e,'params','V'),   T.entry);
        KAv= cellfun(@(e) fget(e,'params','KA'),  T.entry);
        KBv= cellfun(@(e) fget(e,'params','KB'),  T.entry);
        KGv= cellfun(@(e) fget(e,'params','KG'),  T.entry);

        mask = isSeed & Av==MP.A & Vv==MP.V & KAv==MP.KA & KBv==MP.KB & KGv==MP.KG;
        if any(mask)
            ix = find(mask,1,'first');
            src = ""; if isfield(T.entry{ix}.meta,'source'), src = string(T.entry{ix}.meta.source); end
            ishapesDir = fullfile(fileparts(simDir),'InitialShapes');
            f = src; if ~isempty(src) && exist(f,'file')~=2, f = fullfile(ishapesDir, src); end
            if exist(f,'file')==2
                tmp = load(f);
                if isfield(tmp,'Version') && numel(tmp.Version)>=1 && isfield(tmp.Version(1),'Solution')
                    warm.sol = tmp.Version(1).Solution;
                    return
                end
            end
        end
    end
end

% ---- local helpers ----
function v = fget(e, a, b)
    if isstruct(e) && isfield(e,a) && isstruct(e.(a)) && isfield(e.(a),b)
        v = double(e.(a).(b)); else, v = NaN; end
end