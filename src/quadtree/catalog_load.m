function T = catalog_load(simDir)
    catF = fullfile(simDir,'catalog.mat');
    if ~exist(catF,'file')
        dt = datetime.empty(0,1);
        dt.TimeZone = 'UTC';
        dt.Format = 'yyyy-MM-dd HH:mm:ss';
        T = table(string.empty(0,1), dt, cell.empty(0,1), ...
                  'VariableNames', ["hash","timestamp","entry"]);
        save(catF,'T');
        return;
    end

    S = load(catF,'T');
    T = S.T;

    % Normalize any legacy wide catalogs into the 3-col form
    if ~any(T.Properties.VariableNames == "entry")
        T = fold_legacy_wide_to_entry(T);
        save(catF,'T');
    end

    % ENFORCE: entry is a cell column
    if isstruct(T.entry)
        T.entry = num2cell(T.entry);
    end

    % ENFORCE: UTC timestamps
    if ~isempty(T.timestamp)
        if isempty(T.timestamp.TimeZone), T.timestamp.TimeZone = 'UTC'; end
        T.timestamp.Format = 'yyyy-MM-dd HH:mm:ss';
    end
end

function T = fold_legacy_wide_to_entry(Tw)
% FOLD_LEGACY_WIDE_TO_ENTRY
% Converts a table with columns like H0_1,H0_2,A,V,KA,KB,KG,label,E,P,... into:
%   hash | timestamp | entry{r} = struct('params',P,'meta',M)

    if ~any(Tw.Properties.VariableNames=="hash")
        error('catalog_load:legacy', 'Legacy catalog missing "hash" column.');
    end
    if ~any(Tw.Properties.VariableNames=="timestamp")
        % Create placeholder timestamps
        dt = repmat(datetime('now','TimeZone','UTC'), height(Tw), 1);
        dt.Format = 'yyyy-MM-dd HH:mm:ss';
        Tw.timestamp = dt;
    end

    T = table;
    T.hash = string(Tw.hash);
    T.timestamp = Tw.timestamp;

    paramsFields = ["H0_1","H0_2","A","V","KA","KB","KG","aS","bS"];
    metaFields   = ["label","E","P","BCmax","DEmax","mesh","rMinAway","type","source","version","model_version"];

    E = repmat(struct('params',struct(),'meta',struct()), height(Tw), 1);

    varNames = string(Tw.Properties.VariableNames);

    for r = 1:height(Tw)
        P = struct(); M = struct();

        for c = 1:numel(varNames)
            nm = varNames(c);
            if nm=="hash" || nm=="timestamp", continue; end

            val = Tw{r, nm};
            if iscell(val), val = val{1}; end

            if any(nm==paramsFields)
                P.(char(nm)) = val;
            elseif any(nm==metaFields)
                M.(char(nm)) = val;
            else
                if ~isfield(M,'extra'), M.extra = struct(); end
                M.extra.(char(nm)) = val;
            end
        end

        % Ensure meta.hash exists for warm-start loading
        M.hash = string(T.hash(r));

        E(r).params = P;
        E(r).meta   = M;
    end

    T.entry = num2cell(E);
end
