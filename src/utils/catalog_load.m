% File: src/utils/catalog_io.m
% Canonical 3-col MAT-only catalog: hash, timestamp, entry(params,meta)

function T = catalog_load(simDir)
    catF = fullfile(simDir,'catalog.mat');
    if ~exist(catF,'file')
        T = table(string.empty(0,1), datetime.empty(0,1), repmat(struct('params',struct(),'meta',struct()),0,1), ...
                  'VariableNames', ["hash","timestamp","entry"]);
        save(catF,'T'); return;
    end
    S = load(catF,'T');
    T = S.T;
    % Normalize any legacy wide catalogs into the 3-col form
    if ~any(T.Properties.VariableNames == "entry")
        T = fold_legacy_wide_to_entry(T);
        save(catF,'T');
    end
end

function T = fold_legacy_wide_to_entry(Tw)
    % Map your 15-column legacy schema into entry.params/meta
    % Expected columns (case-sensitive): 
    % timestamp, hash, H0_1, H0_2, A, V, KA, KB, KG, label, E, P, BCmax, DEmax, mesh
    names = Tw.Properties.VariableNames;
    assert(all(ismember({'timestamp','hash'}, names)), 'Legacy catalog missing timestamp/hash.');
    T = table; T.hash = string(Tw.hash); T.timestamp = Tw.timestamp;
    E = repmat(struct('params',struct(),'meta',struct()), height(Tw), 1);
    for r = 1:height(Tw)
        P = struct(); M = struct();
        P.H0_1 = safe_get(Tw, 'H0_1', r);
        P.H0_2 = safe_get(Tw, 'H0_2', r);
        P.A    = safe_get(Tw, 'A',    r);
        P.V    = safe_get(Tw, 'V',    r);
        P.KA   = safe_get(Tw, 'KA',   r);
        P.KB   = safe_get(Tw, 'KB',   r);
        P.KG   = safe_get(Tw, 'KG',   r);

        M.label      = safe_get(Tw, 'label',  r, "string");
        M.E          = safe_get(Tw, 'E',      r);
        M.P          = safe_get(Tw, 'P',      r);
        M.BCmax      = safe_get(Tw, 'BCmax',  r);
        M.DEmax      = safe_get(Tw, 'DEmax',  r);
        M.mesh       = safe_get(Tw, 'mesh',   r);
        M.version    = "legacy";   % tag as legacy unless overwritten later
        M.hash       = string(T.hash(r));

        E(r).params = P;
        E(r).meta   = M;
    end
    T.entry = E;
end

function v = safe_get(T, name, r, kind)
    if nargin < 4, kind = "double"; end
    if any(strcmp(T.Properties.VariableNames, name))
        v = T.(name)(r);
        if kind == "string", v = string(v); end
    else
        if kind == "string"
            v = "";
        else
            v = NaN;
        end
    end
end
