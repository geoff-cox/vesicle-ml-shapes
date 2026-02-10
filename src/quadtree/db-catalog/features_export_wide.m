function T = features_export_wide(db, pack_id, varargin)
%FEATURES_EXPORT_WIDE Pivot feature_values (long) into a wide MATLAB table.
%
%   T = features_export_wide(db, pack_id)
%   T = features_export_wide(db, pack_id, 'OnlyOK', true, 'OnlyBest', true, 'Typed', true)
%
% Output table includes solution metadata columns plus one column per feature name.
%
% Options:
%   OnlyOK   (default true): include only solutions.status='ok'
%   OnlyBest (default true): include only solutions.is_best=1
%   Typed    (default true): decode value_real/value_int/value_text into a single MATLAB column
%                            by preferring real then int then text; if false, returns only value_real.
%
% Notes:
% - Columns are created dynamically based on feature_defs for the pack.
% - Joins solutions + param_points to include H0_1/H0_2 and physics params.
% - For features not present for a solution, values are missing/NaN/""
% This pivots feature_values into a wide MATLAB table.
%
% It supports:
%
% * filtering to only status='ok'
% * choosing whether to use best solutions only (is_best=1)
% * returning either only value_real or “typed” columns (real/int/text) (default is typed)
%
% Example usage:
%
%    T = features_export_wide(db, "geom_v1", 'OnlyOK', true, 'OnlyBest', true, 'Typed', true);
%
%    % Quick sanity plot: e.g., r_neck across parameter space
%    scatter(T.H0_1, T.H0_2, 30, T.r_neck, 'filled'); colorbar;
%    title('r\_neck over (H0\_1,H0\_2)');
%

    p = inputParser;
    addParameter(p,'OnlyOK',true,@islogical);
    addParameter(p,'OnlyBest',true,@islogical);
    addParameter(p,'Typed',true,@islogical);
    parse(p,varargin{:});

    onlyOK = p.Results.OnlyOK;
    onlyBest = p.Results.OnlyBest;
    typed = p.Results.Typed;

    pack_id = string(pack_id);

    % 1) Get feature names (and dtype) in stable order
    defSql = "SELECT name, dtype FROM feature_defs WHERE pack_id=" + sql_lit(pack_id) + ...
             " ORDER BY COALESCE(ord, 1e9), name;";
    defs = fetch(db, defSql);

    if isempty(defs)
        error('features_export_wide:nopack','No feature_defs found for pack_id=%s', pack_id);
    end

    featNames = string(defs(:,1));
    featDtype = string(defs(:,2));

    % 2) Pull base solution rows (best/ok filter)
    where = "WHERE 1=1 ";
    if onlyOK
        where = where + "AND s.status='ok' ";
    end
    if onlyBest
        where = where + "AND s.is_best=1 ";
    end

    baseSql = ...
        "SELECT " + ...
        "s.solution_id, s.run_id, s.param_hash, " + ...
        "p.model_version, p.H0_1, p.H0_2, p.A, p.V, p.KA, p.KB, p.KG, " + ...
        "s.label, s.energy_E, s.pressure_P, s.bc_max, s.de_max, s.mesh_n, s.wall_time_s, s.created_at " + ...
        "FROM solutions s JOIN param_points p ON p.param_hash = s.param_hash " + ...
        where + ...
        "ORDER BY p.H0_1, p.H0_2, s.created_at;";

    base = fetch(db, baseSql);
    if isempty(base)
        T = table();
        return;
    end

    % Convert base into MATLAB table
    baseCols = ["solution_id","run_id","param_hash","model_version","H0_1","H0_2","A","V","KA","KB","KG", ...
                "label","energy_E","pressure_P","bc_max","de_max","mesh_n","wall_time_s","created_at"];
    T = cell2table(base, 'VariableNames', cellstr(baseCols));

    % Normalize types for numeric columns
    numCols = ["H0_1","H0_2","A","V","KA","KB","KG","energy_E","pressure_P","bc_max","de_max","mesh_n","wall_time_s"];
    for c = numCols
        if iscell(T.(c))
            T.(c) = cellfun(@toDoubleOrNaN, T.(c));
        end
    end
    T.mesh_n = round(T.mesh_n);

    % 3) Fetch feature values for this pack for all selected solutions
    idList = string(T.solution_id);
    idIn = sql_in_list(idList);

    if typed
        fvSql = "SELECT solution_id, name, value_real, value_int, value_text " + ...
                "FROM feature_values WHERE pack_id=" + sql_lit(pack_id) + " AND solution_id IN " + idIn + ";";
    else
        fvSql = "SELECT solution_id, name, value_real " + ...
                "FROM feature_values WHERE pack_id=" + sql_lit(pack_id) + " AND solution_id IN " + idIn + ";";
    end

    fv = fetch(db, fvSql);
    if isempty(fv)
        % Create empty feature columns
        T = add_empty_feature_columns(T, featNames, featDtype, typed);
        return;
    end

    % 4) Initialize columns
    T = add_empty_feature_columns(T, featNames, featDtype, typed);

    % 5) Map solution_id -> row index
    idToRow = containers.Map('KeyType','char','ValueType','int32');
    for i=1:numel(idList)
        idToRow(char(idList(i))) = int32(i);
    end

    % 6) Fill values
    for r = 1:size(fv,1)
        sid = string(fv{r,1});
        nm  = string(fv{r,2});
        if ~isKey(idToRow, char(sid)), continue; end
        i = idToRow(char(sid));

        if typed
            v = prefer_typed_value(fv{r,3}, fv{r,4}, fv{r,5});
        else
            v = fv{r,3};
        end

        vn = matlab.lang.makeValidName(nm);
        T.(vn)(i) = cast_into_column(v, T.(vn));
    end
end

% ---------- helpers ----------
function x = toDoubleOrNaN(v)
    if isempty(v), x = NaN; return; end
    if isnumeric(v), x = double(v); return; end
    if ischar(v) || isstring(v)
        t = str2double(v);
        if isnan(t), x = NaN; else, x = t; end
        return;
    end
    x = NaN;
end

function in = sql_in_list(ids)
    % Build: ('a','b','c')
    lits = arrayfun(@(s) sql_lit(s), ids);
    in = "(" + strjoin(lits, ",") + ")";
end

function T = add_empty_feature_columns(T, featNames, featDtype, typed)
    for k=1:numel(featNames)
        nm = matlab.lang.makeValidName(featNames(k));
        dt = featDtype(k);

        if ~typed
            T.(nm) = NaN(height(T),1);
            continue;
        end

        switch dt
            case "real"
                T.(nm) = NaN(height(T),1);
            case "int"
                T.(nm) = NaN(height(T),1); % store as double for plotting ease
            case "text"
                T.(nm) = strings(height(T),1);
            otherwise
                % blob not supported here (keep placeholder)
                T.(nm) = strings(height(T),1);
        end
    end
end

function v = prefer_typed_value(vr, vi, vt)
    if ~isempty(vr) && ~(isnumeric(vr) && ~isfinite(double(vr)))
        v = vr; return;
    end
    if ~isempty(vi)
        v = vi; return;
    end
    if ~isempty(vt)
        v = vt; return;
    end
    v = [];
end

function out = cast_into_column(v, col)
    % Return a scalar assignable to col(i)
    if isempty(v)
        if isstring(col), out = ""; else, out = NaN; end
        return;
    end

    if isstring(col)
        out = string(v);
    else
        if isnumeric(v)
            out = double(v);
        else
            t = str2double(string(v));
            if isnan(t), out = NaN; else, out = t; end
        end
    end
end
