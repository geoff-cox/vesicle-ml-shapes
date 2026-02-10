function catalog_upsert_feature_value(db, solution_id, pack_id, name, value, computed_at)
%CATALOG_UPSERT_FEATURE_VALUE Store one feature in long format.
% value may be numeric scalar, integer scalar, string/char, uint8 vector (blob).

    if nargin < 6 || isempty(computed_at)
        computed_at = iso_utc(datetime('now','TimeZone','UTC'));
    end

    vr = "NULL"; vi = "NULL"; vt = "NULL"; vb = "NULL";

    if isnumeric(value) && isscalar(value)
        if isfinite(value)
            vr = sql_lit(double(value));
        end
    elseif islogical(value) && isscalar(value)
        vi = sql_lit(int64(value));
    elseif (isstring(value) || ischar(value))
        vt = sql_lit(string(value));
    elseif isa(value,'uint8')
        % Store as hex string in SQL literal? Better: store NULL here and prefer artifacts or vectors.
        % If you truly want blobs in SQLite from MATLAB, use a binary parameter interface (not used here).
        error('catalog_upsert_feature_value:blob',...
              'BLOB insertion via pure SQL literals is not supported in this helper. Use feature_vectors or artifact files.');
    else
        error('catalog_upsert_feature_value:unsupported','Unsupported feature value type: %s', class(value));
    end

    sql = "INSERT INTO feature_values(solution_id,pack_id,name,value_real,value_int,value_text,value_blob,computed_at) VALUES (" + ...
          sql_lit(solution_id) + "," + sql_lit(pack_id) + "," + sql_lit(name) + "," + ...
          vr + "," + vi + "," + vt + "," + vb + "," + sql_lit(computed_at) + ...
          ") ON CONFLICT(solution_id,pack_id,name) DO UPDATE SET " + ...
          "value_real=excluded.value_real, value_int=excluded.value_int, value_text=excluded.value_text, " + ...
          "computed_at=excluded.computed_at;";

    db_exec(db, sql);
end
