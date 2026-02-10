function artifact_id = catalog_add_artifact(db, solution_id, kind, format, relpath, sha256, bytes)
%CATALOG_ADD_ARTIFACT Register an artifact (e.g., solution HDF5 file).

    artifact_id = new_uuid();
    created_at = iso_utc(datetime('now','TimeZone','UTC'));

    sql = "INSERT INTO artifacts(artifact_id,solution_id,kind,format,path,bytes,sha256,created_at) VALUES (" + ...
          sql_lit(artifact_id) + "," + sql_lit(solution_id) + "," + sql_lit(kind) + "," + sql_lit(format) + "," + ...
          sql_lit(relpath) + "," + sql_lit_or_null(bytes) + "," + sql_lit_or_null(sha256) + "," + sql_lit(created_at) + ");";

    db_exec(db, sql);
end

function s = sql_lit_or_null(x)
    if isempty(x), s="NULL"; else, s=sql_lit(x); end
end
