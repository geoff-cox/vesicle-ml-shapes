% Assume db already open and solution_id exists.
pack_id = "embed_v1";
vector_name = "umap_d2";
vec = single([0.123; -0.456]);   % 2D embedding

% Register pack (one-time)
catalog_insert_feature_pack(db, pack_id, "Embedding vectors", "pack: embed_v1");

% Upsert vector (BLOB)
catalog_upsert_feature_vector(db, solution_id, pack_id, vector_name, vec, "float32", struct('meaning',"UMAP 2D"));

% Fetch it back
v2 = catalog_fetch_feature_vector(db, solution_id, pack_id, vector_name);
