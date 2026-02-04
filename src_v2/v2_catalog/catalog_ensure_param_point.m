function param_hash = catalog_ensure_param_point(db, model_version, P)
%CATALOG_ENSURE_PARAM_POINT Insert param_points row if missing; return param_hash.

    param_hash = param_hash_v2(model_version, P);

    % Fast existence check
    q = "SELECT param_hash FROM param_points WHERE param_hash=" + sql_lit(param_hash) + " LIMIT 1;";
    hit = db_fetch1(db, q);
    if ~isempty(hit)
        param_hash = string(hit);
        return;
    end

    created_at = iso_utc(datetime('now','TimeZone','UTC'));
    sql = "INSERT INTO param_points(param_hash,model_version,H0_1,H0_2,A,V,KA,KB,KG,created_at) VALUES (" + ...
          sql_lit(param_hash) + "," + sql_lit(string(model_version)) + "," + ...
          sql_lit(P.H0_1) + "," + sql_lit(P.H0_2) + "," + sql_lit(P.A) + "," + sql_lit(P.V) + "," + ...
          sql_lit(P.KA) + "," + sql_lit(P.KB) + "," + sql_lit(P.KG) + "," + sql_lit(created_at) + ");";

    % Insert may race; tolerate UNIQUE failure
    try
        db_exec(db, sql);
    catch ME
        if ~contains(ME.message, "UNIQUE", "IgnoreCase", true)
            rethrow(ME);
        end
    end
end
