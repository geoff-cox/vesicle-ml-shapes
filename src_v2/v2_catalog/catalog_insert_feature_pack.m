function catalog_insert_feature_pack(db, pack_id, description, registry_yaml)
%CATALOG_INSERT_FEATURE_PACK Insert a feature pack (idempotent).
    created_at = iso_utc(datetime('now','TimeZone','UTC'));
    reg_sha = hash_sha256(registry_yaml);

    sql = "INSERT OR IGNORE INTO feature_packs(pack_id,created_at,description,registry_yaml,registry_sha256) VALUES (" + ...
          sql_lit(pack_id) + "," + sql_lit(created_at) + "," + sql_lit(description) + "," + ...
          sql_lit(registry_yaml) + "," + sql_lit(reg_sha) + ");";
    db_exec(db, sql);
end
