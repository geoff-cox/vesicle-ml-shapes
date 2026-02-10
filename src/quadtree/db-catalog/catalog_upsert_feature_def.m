function catalog_upsert_feature_def(db, pack_id, name, dtype, unit, description, ord)
%CATALOG_UPSERT_FEATURE_DEF Insert/update a feature definition.

    if nargin < 6, description = ""; end
    if nargin < 7, ord = []; end

    sql = "INSERT INTO feature_defs(pack_id,name,dtype,unit,description,ord) VALUES (" + ...
          sql_lit(pack_id) + "," + sql_lit(name) + "," + sql_lit(dtype) + "," + ...
          sql_lit_or_null(unit) + "," + sql_lit_or_null(description) + "," + sql_lit_or_null(ord) + ...
          ") ON CONFLICT(pack_id,name) DO UPDATE SET " + ...
          "dtype=excluded.dtype, unit=excluded.unit, description=excluded.description, ord=excluded.ord;";
    db_exec(db, sql);
end

function s = sql_lit_or_null(x)
    if isempty(x), s="NULL"; else, s=sql_lit(x); end
end
