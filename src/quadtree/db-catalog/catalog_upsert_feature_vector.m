function catalog_upsert_feature_vector(db, solution_id, pack_id, vector_name, vec, dtype, order_struct)
%CATALOG_UPSERT_FEATURE_VECTOR Insert/update vector into feature_vectors.
%
% vec: numeric vector
% dtype: 'float32'|'float64'|'int32'|'int64'
% order_struct: optional struct encoded into JSON (e.g., names, meaning)

    if nargin < 7, order_struct = struct(); end

    computed_at = iso_utc(datetime('now','TimeZone','UTC'));

    [u8, dtype, n] = pack_vector_blob(vec, dtype);
    blobLit = sql_blob_lit(u8);

    [order_json, ~] = canon_json(order_struct);

    sql = "INSERT INTO feature_vectors(solution_id,pack_id,vector_name,dtype,length,blob,order_json,computed_at) VALUES (" + ...
          sql_lit(solution_id) + "," + sql_lit(pack_id) + "," + sql_lit(vector_name) + "," + ...
          sql_lit(dtype) + "," + sql_lit(n) + "," + blobLit + "," + sql_lit(order_json) + "," + sql_lit(computed_at) + ...
          ") ON CONFLICT(solution_id,pack_id,vector_name) DO UPDATE SET " + ...
          "dtype=excluded.dtype, length=excluded.length, blob=excluded.blob, order_json=excluded.order_json, computed_at=excluded.computed_at;";

    exec(db, sql);
end
