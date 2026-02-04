-- Vesicle Catalog Schema v2 (SQLite)
-- Authoritative scalar index + provenance + features + artifacts pointers.
-- All timestamps are ISO-8601 UTC strings (e.g., 2026-02-03T18:42:10Z).

PRAGMA foreign_keys = ON;

-- ----------------------------
-- RUNS (provenance / config)
-- ----------------------------
CREATE TABLE IF NOT EXISTS runs (
  run_id        TEXT PRIMARY KEY,
  started_at    TEXT NOT NULL,
  finished_at   TEXT,
  status        TEXT NOT NULL CHECK(status IN ('running','finished','aborted')),
  model_version TEXT NOT NULL,

  -- exploration config
  h0_1_min      REAL NOT NULL,
  h0_1_max      REAL NOT NULL,
  h0_2_min      REAL NOT NULL,
  h0_2_max      REAL NOT NULL,
  qt_max_depth  INTEGER NOT NULL,
  qt_max_cells  INTEGER NOT NULL,

  -- optional global physics (if constant across run)
  A            REAL,
  V            REAL,
  KA           REAL,
  KB           REAL,
  KG           REAL,

  -- canonical JSON snapshot of SP/TH/MP + any metadata
  sim_json      TEXT NOT NULL,
  sim_json_sha256 TEXT NOT NULL,

  notes         TEXT
);

CREATE INDEX IF NOT EXISTS runs_started_idx ON runs(started_at);
CREATE INDEX IF NOT EXISTS runs_status_idx  ON runs(status);

-- ----------------------------
-- SEEDS (optional but useful)
-- ----------------------------
CREATE TABLE IF NOT EXISTS seeds (
  seed_id       TEXT PRIMARY KEY,
  name          TEXT NOT NULL,
  source        TEXT,                 -- e.g. file name
  seed_hash     TEXT NOT NULL UNIQUE, -- deterministic hash of seed identity/content if desired
  created_at    TEXT NOT NULL,
  seed_json     TEXT,                 -- optional: seed metadata/config
  notes         TEXT
);

-- ----------------------------------------
-- PARAMETER POINTS (dedup anchor)
-- One row per unique equation-defining params
-- ----------------------------------------
CREATE TABLE IF NOT EXISTS param_points (
  param_hash    TEXT PRIMARY KEY,   -- SHA-256 hex
  model_version TEXT NOT NULL,

  H0_1          REAL NOT NULL,
  H0_2          REAL NOT NULL,
  A             REAL NOT NULL,
  V             REAL NOT NULL,
  KA            REAL NOT NULL,
  KB            REAL NOT NULL,
  KG            REAL NOT NULL,

  created_at    TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS param_points_h0_idx ON param_points(H0_1, H0_2);

-- ----------------------------------------
-- SOLUTIONS (multiple per param_hash)
-- Each row is one attempt/branch/solution artifact
-- ----------------------------------------
CREATE TABLE IF NOT EXISTS solutions (
  solution_id   TEXT PRIMARY KEY,
  run_id        TEXT NOT NULL REFERENCES runs(run_id) ON DELETE CASCADE,
  param_hash    TEXT NOT NULL REFERENCES param_points(param_hash) ON DELETE CASCADE,

  -- branch provenance (optional but strongly recommended)
  seed_hash       TEXT,      -- references seeds.seed_hash if used
  delta_path_hash TEXT,      -- hash of homotopy schedule/options (optional)
  branch_tag      TEXT,      -- free-form: "branchA", "bud", etc.
  parent_solution_id TEXT,   -- continuation parent (optional)

  -- status
  status        TEXT NOT NULL CHECK(status IN ('ok','fail','reject')),
  is_best       INTEGER NOT NULL DEFAULT 0 CHECK(is_best IN (0,1)),

  fail_class    TEXT,      -- e.g. 'bvp_fail'|'nan'|'singular'|'timeout'
  fail_msg      TEXT,

  -- scalar science outputs (optional, but common)
  label         TEXT,
  energy_E      REAL,
  pressure_P    REAL,

  -- solver diagnostics
  bc_max        REAL,
  de_max        REAL,
  mesh_n        INTEGER,
  n_refines     INTEGER,
  wall_time_s   REAL,

  created_at    TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS solutions_run_idx    ON solutions(run_id);
CREATE INDEX IF NOT EXISTS solutions_param_idx  ON solutions(param_hash);
CREATE INDEX IF NOT EXISTS solutions_status_idx ON solutions(status);
CREATE INDEX IF NOT EXISTS solutions_label_idx  ON solutions(label);

-- At most one "best" solution per param point:
CREATE UNIQUE INDEX IF NOT EXISTS solutions_best_per_param
  ON solutions(param_hash)
  WHERE is_best = 1;

-- ----------------------------------------
-- ARTIFACTS (HDF5, plots, embeddings, etc.)
-- ----------------------------------------
CREATE TABLE IF NOT EXISTS artifacts (
  artifact_id   TEXT PRIMARY KEY,
  solution_id   TEXT NOT NULL REFERENCES solutions(solution_id) ON DELETE CASCADE,

  kind          TEXT NOT NULL,  -- 'solution_h5'|'embedding'|'plot'|'other'
  format        TEXT NOT NULL,  -- 'h5'|'png'|'bin'|'json'|...
  path          TEXT NOT NULL,  -- relative path under results root
  bytes         INTEGER,
  sha256        TEXT,

  created_at    TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS artifacts_solution_idx ON artifacts(solution_id);
CREATE INDEX IF NOT EXISTS artifacts_kind_idx     ON artifacts(kind);

-- ----------------------------------------
-- FEATURE PACKS (versioned)
-- ----------------------------------------
CREATE TABLE IF NOT EXISTS feature_packs (
  pack_id         TEXT PRIMARY KEY,   -- e.g. "geom_v1"
  created_at      TEXT NOT NULL,
  description     TEXT NOT NULL,
  registry_yaml   TEXT NOT NULL,
  registry_sha256 TEXT NOT NULL
);

-- Feature definitions (the contract)
CREATE TABLE IF NOT EXISTS feature_defs (
  pack_id      TEXT NOT NULL REFERENCES feature_packs(pack_id) ON DELETE CASCADE,
  name         TEXT NOT NULL,
  dtype        TEXT NOT NULL CHECK(dtype IN ('real','int','text','blob')),
  unit         TEXT,
  description  TEXT,
  ord          INTEGER,     -- optional ordering for vectors / exports
  PRIMARY KEY(pack_id, name)
);

-- Feature values (long format; flexible & stable)
CREATE TABLE IF NOT EXISTS feature_values (
  solution_id  TEXT NOT NULL REFERENCES solutions(solution_id) ON DELETE CASCADE,
  pack_id      TEXT NOT NULL REFERENCES feature_packs(pack_id) ON DELETE CASCADE,
  name         TEXT NOT NULL,

  value_real   REAL,
  value_int    INTEGER,
  value_text   TEXT,
  value_blob   BLOB,

  computed_at  TEXT NOT NULL,

  PRIMARY KEY(solution_id, pack_id, name),
  FOREIGN KEY(pack_id, name) REFERENCES feature_defs(pack_id, name) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS feature_values_pack_idx ON feature_values(pack_id);
CREATE INDEX IF NOT EXISTS feature_values_name_idx ON feature_values(name);

-- Optional: dense vectors (embeddings, descriptors) stored directly in DB
CREATE TABLE IF NOT EXISTS feature_vectors (
  solution_id  TEXT NOT NULL REFERENCES solutions(solution_id) ON DELETE CASCADE,
  pack_id      TEXT NOT NULL REFERENCES feature_packs(pack_id) ON DELETE CASCADE,
  vector_name  TEXT NOT NULL,              -- e.g. "umap_d2", "pca_d8"
  dtype        TEXT NOT NULL CHECK(dtype IN ('float32','float64','int32','int64')),
  length       INTEGER NOT NULL,
  blob         BLOB NOT NULL,
  order_json   TEXT,                       -- optional: names/ordering
  computed_at  TEXT NOT NULL,
  PRIMARY KEY(solution_id, pack_id, vector_name)
);

-- ----------------------------------------
-- Scheduler state (restartability)
-- ----------------------------------------
CREATE TABLE IF NOT EXISTS scheduler_state (
  key          TEXT PRIMARY KEY,
  value_json   TEXT NOT NULL,
  updated_at   TEXT NOT NULL
);

-- Convenience view: “best” solutions joined with param_points
CREATE VIEW IF NOT EXISTS v_best_solutions AS
SELECT
  s.solution_id, s.run_id, s.param_hash,
  p.model_version, p.H0_1, p.H0_2, p.A, p.V, p.KA, p.KB, p.KG,
  s.seed_hash, s.delta_path_hash, s.branch_tag, s.parent_solution_id,
  s.status, s.label, s.energy_E, s.pressure_P,
  s.bc_max, s.de_max, s.mesh_n, s.n_refines, s.wall_time_s,
  s.created_at
FROM solutions s
JOIN param_points p ON p.param_hash = s.param_hash
WHERE s.is_best = 1;
