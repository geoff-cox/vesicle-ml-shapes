# GitHub Copilot Instructions for vesicle-ml-shapes

## Repository Overview
**Purpose**: Computational framework for exploring equilibrium shapes of two-phase vesicles using MATLAB numerical continuation and machine learning.  
**Size**: ~13 MB, 45 MATLAB files, 1 Python template  
**Languages**: MATLAB (primary, R2020b+), Python 3.x (ML pipeline - templates only, not yet functional)  
**No special toolboxes required** - uses only base MATLAB installation.

## Critical: Environment & Build

### MATLAB Environment Setup
**ALWAYS run these commands at the start of any MATLAB session:**
```matlab
restoredefaultpath; rehash toolboxcache;
addpath(genpath(fullfile(pwd,'src')));
addpath(fullfile(pwd,'bvp6c-solver'));
addpath(fullfile(pwd,'initial-shapes'));
if ~exist('sim-results','dir'), mkdir('sim-results'); end
```
**Why**: MATLAB won't find functions without proper path setup. The above is from `script_driver_slim.m` bootstrap section.

### Running Simulations
**To run the main simulation (if MATLAB available):**
```matlab
script_driver_slim
```
This bootstraps paths, initializes catalog from seed shapes in `initial-shapes/`, and runs parameter exploration. **Note**: MATLAB may not be available in CI environments. This is a research codebase primarily for local execution.

### Python Environment (ML Pipeline)
**Status**: Python code in `notebooks/` is **TEMPLATE ONLY - NOT FUNCTIONAL**. See `notebooks/README.md`.  
**If implementing**: Would require `pip install numpy pandas scipy scikit-learn umap-learn matplotlib pyarrow`

### No CI/CD Pipeline
**There are NO GitHub Actions workflows or automated tests.** The `.github/` directory only contains this instructions file. Validation must be manual.

## Project Structure & Key Files

### Directory Layout
```
├── script_driver_slim.m          # MAIN ENTRY POINT - run simulations
├── script_driver.mlx             # Live Script version (interactive)
├── bvp6c-solver/                 # 6th-order BVP solver (9 .m files) - DO NOT MODIFY
│   └── bvp6c.m                   # Core solver adapted from MATLAB
├── initial-shapes/               # 4 seed .mat files for continuation bootstrap
│   └── SIM_Node_*.mat           # Named by parameters (A, V, KG, KA, KB, H0_1, H0_2)
├── sim-results/                  # Generated outputs (git-ignored solutions/)
│   ├── catalog.mat              # Master index (hash → params/meta) - 3-column table
│   ├── cache.mat                # Quad-tree state and failure registry
│   └── hashed_results/          # Individual .mat solution files (SHA-256 names)
├── src/
│   ├── sim_explore_H0_quad_tree.m    # Main parameter space explorer (263 lines)
│   ├── utils/                   # 29 utility functions
│   │   ├── solveAtParams_v2.m   # Multi-stage BVP continuation solver (240 lines)
│   │   ├── catalog_*.m          # Catalog I/O functions (load/save/append)
│   │   ├── processQuadtree.m    # Adaptive refinement scheduler
│   │   ├── BendV_Lag_EIGp_*.m   # Physics model (BC & DE implementations)
│   │   └── simpleDataHash.m     # SHA-256 hashing for unique IDs
│   └── tools/                   # 5 maintenance scripts (catalog rebuild, migration)
├── notebooks/                    # ML pipeline design (NOT IMPLEMENTED - see README.md)
│   ├── README.md                # MUST READ - explains template-only status
│   ├── feature-list-data-schema.yaml  # Feature spec for future ML work
│   └── pipeline-notebook.py     # Analysis template (requires data export first)
├── docs/                         # Research PDFs and technical guides
│   └── explore_H0_quad_tree_README.mlx
├── PROJECT_AUDIT.md              # Comprehensive architecture documentation (512 lines)
└── README.md                     # Quick start guide (101 lines)
```

### Files in Repository Root
`.gitignore`, `AUDIT_SUMMARY.md`, `PROJECT_AUDIT.md`, `README.md`, `script_driver.mlx`, `script_driver_slim.m`, `tool_driver.mlx`

## Validation & Testing

### No Automated Tests
**Critical**: There is NO test suite. Manual validation only:
1. **For solver changes**: Run `script_driver_slim` with small `MaxIters` (e.g., 5-10) and verify:
   - No errors/crashes
   - Results appear in `sim-results/hashed_results/`
   - Catalog updates correctly (`catalog.mat` row count increases)
2. **For catalog/utility changes**: Check catalog integrity after changes:
   ```matlab
   T = catalog_load('sim-results');
   disp(T);  % Should show hash, timestamp, entry columns
   ```
3. **Solution quality checks** (manual): Use `bc_diagnostics.m` to verify boundary conditions met.

### Common Validation Pitfalls
- **Path issues**: If MATLAB can't find functions, re-run bootstrap commands above.
- **Catalog corruption**: Always use `catalog_load`, `catalog_append`, `catalog_save` - never directly edit `.mat` files.
- **Hash collisions**: Extremely rare (SHA-256), but if suspected, check for duplicate hashes: `length(unique(T.hash)) == height(T)`

## Build/Run Commands Reference

### MATLAB Commands (Tested Approach)
**None of these can be tested in current environment (MATLAB not installed), but these are verified workflows from codebase:**
1. **Setup session**: See "MATLAB Environment Setup" above
2. **Run simulation**: `script_driver_slim` (modifies `sim_config()` to set `MaxIters`, parameters)
3. **Interactive tools**: `tool_driver.mlx` (MATLAB Live Script - requires GUI)
4. **Maintenance**: Run scripts from `src/tools/` if catalog needs rebuilding

### Python Commands
**Not applicable** - ML code is template-only. Future implementation would need data export from MATLAB first.

## Common Pitfalls & Workarounds

### MATLAB-Specific Issues
1. **"Undefined function" errors**: Missing path setup. Re-run bootstrap section from `script_driver_slim.m`.
2. **BVP solver warnings**: `warning('off','MATLAB:bvp6c:RelTolNotMet')` is set during bootstrap - these are expected for difficult parameter regions.
3. **Solver stalls**: The code has built-in mesh coarsening (`coarsen_mesh.m`) and multiple delta attempts - let it complete.
4. **File locking on Windows**: MATLAB may lock `.mat` files. Close all figures/workspaces before re-running.

### Data Management
1. **Never edit** `catalog.mat` or `cache.mat` manually - use provided functions.
2. **Hash consistency**: Changing physics parameters (A, V, KA, KB, KG) or hash logic in `simpleDataHash.m` invalidates entire catalog. Don't modify unless necessary.
3. **Seed shapes**: The 4 `.mat` files in `initial-shapes/` are imported to catalog on first run via `import_initial_shapes_into_catalog()` in `script_driver_slim.m`.

### Git Ignore
Repository-level `.gitignore` is currently configured to exclude only `sim-results/solutions` and `sim-results/delete_these` (legacy result paths).

## Key Architecture Facts

### Simulation Flow
1. **Bootstrap** → 2. **Config** (`sim.MP`, `sim.TH`, `sim.SP`) → 3. **Quad-tree loop** → 4. **Hash key generation** → 5. **Catalog lookup** → 6. **Warm-start selection** → 7. **BVP solve** (`solveAtParams_v2`) → 8. **Persist result** → 9. **Update catalog/cache**

### Critical Structures
- **`sim.MP`**: Physical params (A, V, KA, KB, KG) - fixed per run. H0_1, H0_2 vary per iteration.
- **`sim.TH`**: Solver thresholds (BCmax=1e-6, delta=0.01, RelTol=1e-6, AbsTol=1e-8, etc.)
- **`sim.SP`**: Policy (MaxIters, ModelVersion="BVP-v3.1", Verbose, LogToFile, etc.)
- **`catalog`**: MATLAB table with columns `hash` (string), `timestamp` (datetime), `entry` (struct with `.params` and `.meta`)
- **`cache`**: Struct with `.frontier` (quad-tree cells), `.failures` (backoff registry), `.config`

### Do Not Modify (Without Deep Understanding)
- **`bvp6c-solver/`**: High-order numerics - changes break convergence.
- **`BendV_Lag_EIGp_*_impl.m`**: Physics model - requires expertise in vesicle mechanics.
- **Hash generation logic**: Changing `simpleDataHash.m` invalidates catalog uniqueness.

## Expected Command Execution Times
- **Single BVP solve**: 0.1-60 seconds (highly variable, depends on parameter region)
- **Full parameter sweep** (1000 points): Hours to days (serial execution)
- **Catalog load/append**: Milliseconds (linear in catalog size, acceptable up to ~10K entries)

## Documentation Priority Order
1. **Start here**: This file
2. **Architecture deep-dive**: `PROJECT_AUDIT.md` (comprehensive, 512 lines)
3. **Quick reference**: `README.md` (basic usage)
4. **ML status**: `notebooks/README.md` (explains template-only state)
5. **Research context**: `docs/explore_H0_quad_tree_README.mlx` (quad-tree strategy)
6. **Code comments**: Inline docs in `src/utils/*.m` (most functions have headers)

## Trusting These Instructions
**These instructions are based on thorough codebase analysis.** Only search for additional information if you find errors, inconsistencies, or gaps for your specific task. The codebase structure is stable and well-documented via `PROJECT_AUDIT.md`.
