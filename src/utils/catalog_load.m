function T = catalog_load(simDir)
    catF = fullfile(simDir,'catalog.mat');
    if ~exist(catF,'file')
        dt = datetime.empty(0,1);
        dt.TimeZone = 'UTC';
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
end

function T = fold_legacy_wide_to_entry(Tw)
    % ... same as your version up to creation of E(r).params/meta ...

    T = table; 
    T.hash = string(Tw.hash); 
    T.timestamp = Tw.timestamp;

    E = repmat(struct('params',struct(),'meta',struct()), height(Tw), 1);
    for r = 1:height(Tw)
        % build P, M ...
        E(r).params = P;
        E(r).meta   = M;
    end

    T.entry = num2cell(E);  % <-- critical
end
