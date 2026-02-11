# Quadtree Parameter-Space Exploration — Pipeline Walkthrough

This document describes the simulation pipeline that explores equilibrium shapes
of two-phase vesicles over the spontaneous-curvature parameter plane
**(H₀⁽¹⁾, H₀⁽²⁾)** using an adaptive quadtree refinement strategy.

---

## Directory Structure

```
src/
├── script_driver_slim.m        ← ENTRY POINT: bootstrap → configure → run → cleanup
│
├── quadtree/                   ← Adaptive parameter-space exploration
│   ├── sim_explore_H0_quad_tree.m   Main driver loop (task dispatch, I/O, failure handling)
│   └── processQuadtree.m           Quadtree scheduler (cell subdivision, corner selection)
│
├── solver/                     ← BVP continuation solver & physics model
│   ├── solveAtParams.m              Multi-strategy continuation solver with acceptance gates
│   └── computePhaseScales.m         Phase arc-length scaling from area fraction A
│
├── catalog/                    ← Data persistence (MAT-based catalog)
│   ├── catalog_load.m               Load catalog table from disk (with legacy migration)
│   ├── catalog_save.m               Atomic write of catalog to disk
│   └── catalog_append.m             Append/merge entries (supports in-memory batching)
│
├── utils/                      ← Shared utilities
│   ├── simpleDataHash.m             SHA-256 hash for physics parameter keys
│   └── defaultArg.m                 Safe struct field access with defaults
│
├── bvp6c/                      ← 6th-order BVP solver (9 files — DO NOT MODIFY)
│   ├── bvp6c.m                      Core solver (adapted from MATLAB)
│   ├── bvpset.m / bvpget.m          Option handling
│   ├── bvpinit.m / deval.m          Initialization and interpolation
│   └── ntrp*.m / odenumjac.m        Interpolation & Jacobian helpers
│
├── initial-shapes/             ← Seed .mat files for continuation bootstrap
│   └── SIM_Node_*.mat               Pre-computed solutions at H₀ = (0, 0) for various (A, V, KG, KA, KB)
│
└── scripts/                    ← Diagnostics & visualization (run manually)
    ├── script_sanity_checks.m       Fast pipeline validation with mock solver
    ├── run_sanity_checks.m          Real solver smoke test (5 parameter points)
    ├── check_interface_curvature.m  Curvature jump diagnostic
    ├── plot_phase_map.m             Morphology scatter over (H₀⁽¹⁾, H₀⁽²⁾)
    ├── plot_residual_cloud.m        BC vs DE residual acceptance cloud
    └── preview_one_solution.m       2D vesicle profile visualization
```

---

## Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────────┐
│                     script_driver_slim.m                           │
│  ┌──────────┐    ┌──────────────┐    ┌────────────────────────┐   │
│  │bootstrap()│───▶│ sim_config() │───▶│sim_explore_H0_quad_tree│   │
│  └──────────┘    └──────────────┘    └───────────┬────────────┘   │
│   • addpath       • SP (settings)                │                │
│   • import seeds  • TH (thresholds)     ┌────────▼────────┐      │
│   • suppress warn • MP (physics)        │   MAIN LOOP     │      │
│                                         │   (up to 1e5    │      │
│                                         │    iterations)  │      │
│                                         └────────┬────────┘      │
│                                                  │               │
│  ┌──────────┐                                    │               │
│  │ cleanup()│◀───────────────────────────────────┘               │
│  └──────────┘                                                    │
└─────────────────────────────────────────────────────────────────────┘
```

### Main Loop Detail

Each iteration of `sim_explore_H0_quad_tree` follows this sequence:

```
1. processQuadtree(cache, T, MP)
   │
   ├─▶ Returns task = {H0_1, H0_2} for the next unsolved corner
   │   or ∅ if all cells are resolved
   │
   ▼
2. Compute SHA-256 hash of (model_version, H0_1, H0_2, A, V, KA, KB, KG)
   │
   ▼
3. Check: result file already on disk?
   ├─ YES → update catalog in memory → continue
   │
   ▼
4. Check: failure registry (cooldown / permanent skip)?
   ├─ Blocked → rotate queue → continue
   │
   ▼
5. pickWarmStart(params, sim, simDir, T)
   │
   ├─ Search catalog for nearest solved neighbor with matching physics
   ├─ Fall back to seed file from initial-shapes/
   │
   ▼
6. solveAtParams(params, sim, warm)
   │
   ├─ Multi-rung continuation: step caps × delta list × tolerance sets
   ├─ Acceptance gates: BCmax ≤ 1e-6, DEmax ≤ 0.2, rMin ≥ 1e-3
   │
   ├─ SUCCESS → save .mat → append to catalog → clear failure record
   └─ FAILURE → update failure registry with exponential backoff
   │
   ▼
7. Periodic checkpoint: save catalog + cache to disk every SAVE_EVERY iterations
```

### Quadtree Scheduler (`processQuadtree`)

```
┌───────────────────────────────────────────────────┐
│  Start: root cell covers full (H0_1, H0_2) bounds │
│  Default: [-1, 1] × [-1, 1]                       │
└───────────────┬───────────────────────────────────┘
                │
        ┌───────▼───────┐
        │ Pop next cell │◀──────────────────────────────────┐
        │ from queue    │                                   │
        └───────┬───────┘                                   │
                │                                           │
        ┌───────▼────────────┐                              │
        │ Any unsolved       │                              │
        │ corners?           │                              │
        ├─ YES → return task │                              │
        │   (H0_1, H0_2)    │                              │
        │   for that corner  │                              │
        │                    │                              │
        ├─ NO ↓              │                              │
        └────────────────────┘                              │
                │                                           │
        ┌───────▼───────────────┐                           │
        │ Uniform test:         │                           │
        │ • All labels equal?   │                           │
        │ • Energy spread ≤ eTol│                           │
        │ • Pressure spread     │                           │
        │   ≤ pTol              │                           │
        └───────┬───────────────┘                           │
                │                                           │
         ┌──────┴──────┐                                    │
    UNIFORM        NOT UNIFORM                              │
    or max depth   and depth < max                          │
         │              │                                   │
    ┌────▼────┐   ┌─────▼──────┐                            │
    │  Store  │   │ Subdivide  │                            │
    │ in cells│   │ into 4     ├────────────────────────────┘
    │ (done)  │   │ children   │  (enqueue SW, SE, NE, NW)
    └─────────┘   └────────────┘
```

### Solver Strategy (`solveAtParams`)

The solver uses a **multi-rung continuation** approach to reach the target
(H₀⁽¹⁾, H₀⁽²⁾) from a warm-start solution:

```
For each step-cap rung (∞, 0.50, 0.25, 0.12, 0.06, 0.03):
  For each delta (0.01, 0.015, 0.02, 0.008, 0.005):
    For each tolerance set (strict, relaxed):
      March from H0_from → H0_target in fixed-size steps
      At each step: bvp6c solve → check acceptance gates
      If step fails → halve step size (down to minH0Step)

If all rungs fail → error (caught by driver, recorded in failure registry)
```

**Acceptance gates** (all must pass):
| Gate | Threshold | Description |
|------|-----------|-------------|
| `BCmax` | ≤ 1e-6 | Maximum boundary condition residual |
| `DEmax` | ≤ 2e-1 | Maximum ODE residual (interior) |
| `rMin` | ≥ 1e-3 | Minimum radius away from poles |
| `rNeck` | ≥ 0 | Neck radius (disabled by default) |

---

## Data Cataloging

### Design Choice: MAT-Based Catalog

The catalog uses MATLAB's native `.mat` format stored as a 3-column table:

| Column | Type | Description |
|--------|------|-------------|
| `hash` | string | SHA-256 of physics parameters (unique key) |
| `timestamp` | datetime (UTC) | When the entry was created/updated |
| `entry` | cell of struct | `{params: {...}, meta: {...}}` |

**Why MAT-based (not SQLite)?**

1. **Zero dependencies** — no SQLite MEX files or Java bridges needed; works
   on any MATLAB R2020b+ installation out of the box.
2. **Native MATLAB types** — tables, structs, and datetime objects serialize
   without conversion; preserves full numeric precision.
3. **Atomic writes** — uses temp-file + `movefile` pattern to prevent
   corruption on crash or power loss.
4. **In-memory batching** — the driver keeps the catalog in memory and only
   writes to disk every `SAVE_EVERY` iterations (default 50), avoiding
   per-iteration I/O overhead.
5. **Easy export** — `catalog_load` returns a standard MATLAB table, ready
   for analysis, Parquet export, or Python ingestion via `scipy.io.loadmat`.

**Catalog entry structure:**
```matlab
entry.params   % H0_1, H0_2, A, V, KA, KB, KG, aS, bS
entry.meta     % label, E, P, BCmax, DEmax, mesh, rMinAway, rNeck, hash, version
```

### Individual Solution Files

Each accepted solution is saved to `SimResults/hashed_results/<hash>.mat`
containing:
- `result.sol` — full BVP solution struct (mesh, state vector, parameters)
- `meta` — diagnostics, label, energy, pressure, residuals

### Cache File

`SimResults/cache.mat` stores the quadtree state between runs:
- `cache.QT.queue` — pending quadtree cells
- `cache.QT.cells` — resolved (uniform or max-depth) cells
- `cache.failures` — failure registry with exponential backoff
- `cache.config` — bounds, maxDepth, maxCells, tolerances

---

## Running the Simulation

### Quick Start

```matlab
cd src
script_driver_slim
```

### Configuration

Edit the `sim_config()` function in `script_driver_slim.m`:

```matlab
SP.MaxIters = 1e5;          % max iterations (set lower for testing)
MP.A  = 0.50;               % area fraction (α-phase / total)
MP.V  = 0.72;               % reduced volume
MP.KA = 1.0;                % bending modulus, α-phase
MP.KB = 1.0;                % bending modulus, β-phase
MP.KG = 0.0;                % Gaussian bending modulus
```

### Sanity Checks

```matlab
% Fast mock test (no real BVP solves):
cd src/scripts
script_sanity_checks

% Real solver smoke test (5 points near origin):
run_sanity_checks
```

---

## Component Dependency Graph

```
script_driver_slim.m
 ├── computePhaseScales        (solver/)
 ├── bvpset                    (bvp6c/)
 ├── catalog_load              (catalog/)
 ├── catalog_append            (catalog/)
 ├── simpleDataHash            (utils/)
 │
 └── sim_explore_H0_quad_tree  (quadtree/)
      ├── catalog_load / catalog_save / catalog_append   (catalog/)
      ├── simpleDataHash                                  (utils/)
      ├── computePhaseScales                              (solver/)
      ├── defaultArg                                      (utils/)
      │
      ├── processQuadtree                                 (quadtree/)
      │    └── defaultArg                                 (utils/)
      │
      └── solveAtParams                                   (solver/)
           ├── computePhaseScales                         (solver/)
           ├── bvp6c / bvpset / bvpinit / deval           (bvp6c/)
           ├── defaultArg                                 (utils/)
           └── initial-shapes/*.mat                       (initial-shapes/)
```

---

## TODO: Next Steps for Overnight Simulations

### Must Do Before First Overnight Run

- [ ] **Verify seed files** — Run `run_sanity_checks` to confirm the solver
  converges from each seed shape at `H₀ = (0, 0)`.
- [ ] **Set appropriate bounds** — Choose `SP.H0Bounds` based on the physics
  regime of interest (default `[-1, 1] × [-1, 1]` may be too wide or narrow).
- [ ] **Tune `MaxIters`** — For an 8-hour overnight run, estimate iterations:
  `MaxIters ≈ (8 × 3600) / avg_solve_time`. A typical solve takes 1–60 s.
- [ ] **Configure `SAVE_EVERY`** — Set to ~50–100 for overnight runs. Lower
  values increase crash resilience but add I/O overhead.
- [ ] **Test crash recovery** — Kill a short run mid-flight and restart;
  verify the driver resumes correctly from `cache.mat` and skips existing
  `.mat` files.
- [ ] **Disk space check** — Each solution `.mat` file is 50–200 KB.
  Budget ~200 MB per 1,000 solutions.

### Recommended Improvements

- [ ] **Parallel warm-start lookup** — Pre-index the catalog by `(H0_1, H0_2)`
  to avoid `cellfun` scans on every warm-start selection.
- [ ] **Parquet export script** — Write a MATLAB script that converts
  `catalog.mat` to a Parquet file for analysis in Python/pandas.
- [ ] **Live progress dashboard** — Create a script that reads `catalog.mat`
  and plots the quadtree coverage in real time.
- [ ] **Multi-physics sweeps** — Parameterize the driver over `(A, V)` pairs
  to build a comprehensive shape atlas.
- [ ] **Implement ML feature extraction** — Follow `notebooks/feature-list-data-schema.yaml`
  to extract spectral descriptors, curvature statistics, and topological features
  from each solved profile.
- [ ] **Logging** — Implement `logmsg()` to write structured logs to `OPfile.txt`
  for post-mortem analysis.
