# Project Audit Report: vesicle-ml-shapes

**Date:** January 29, 2026  
**Repository:** geoff-cox/vesicle-ml-shapes

---

## Executive Summary

This project aims to systematically explore the equilibrium shape landscape of two-phase vesicles using numerical continuation methods, and to develop machine learning approaches for classifying and understanding these complex morphologies. The codebase consists of a mature MATLAB simulation engine that solves boundary value problems (BVPs) for vesicle shapes, coupled with preliminary Python infrastructure for ML-based shape analysis.

**Current Status:** The simulation framework is well-developed and functional, with initial parameter sweeps completed. The ML pipeline remains in the planning/design phase with schema and notebook scaffolds prepared but not yet implemented.

---

## 1. Project Goals

### Primary Objectives

1. **Systematic Parameter Space Exploration**
   - Map the (H₀₁, H₀₂) spontaneous curvature parameter space for two-phase vesicles
   - Use adaptive quad-tree refinement to efficiently explore regions with high morphological variation
   - Build a comprehensive catalog of equilibrium vesicle shapes under various physical constraints

2. **Shape Classification & Understanding**
   - Extract geometric and energetic features from solved vesicle profiles
   - Apply unsupervised learning (PCA, UMAP, clustering) to discover natural shape families
   - Identify morphological phase boundaries and transitions in parameter space
   - Enable rapid classification of new vesicle shapes

3. **Scientific Insight**
   - Understand how spontaneous curvature differences drive budding, tubulation, and other complex morphologies
   - Characterize hysteresis and multiple solution branches
   - Provide a computational framework for vesicle biophysics research

### Target Applications

- Cell membrane modeling (domain formation, budding, endocytosis)
- Artificial vesicle design
- Understanding lipid bilayer phase separation phenomena

---

## 2. Project Source Code

### 2.1 Overall Structure

```
vesicle-ml-shapes/
├── bvp6c-solver/          # High-order BVP solver (modified from MATLAB)
├── initial-shapes/         # Seed shapes for continuation (4 .mat files)
├── SimResults/            # Simulation outputs and catalog
│   ├── hashed_results/    # Individual solution .mat files (hash-named)
│   ├── catalog.mat        # Master index of all solutions
│   └── cache.mat          # Quad-tree state and work queue
├── src/
│   ├── utils/             # Core simulation utilities (30 .m files)
│   ├── tools/             # Maintenance/migration scripts (5 .m files)
│   └── sim_explore_H0_quad_tree.m  # Main driver
├── docs/                  # Research plan PDF and quad-tree README
├── notebooks/             # ML pipeline design (Python - templates only)
│   ├── README.md                      # Implementation status
│   ├── feature-list-data-schema.yaml  # ML feature specifications
│   └── pipeline-notebook.py           # Starter analysis template
├── script_driver_slim.m   # Runnable simulation entry point
├── script_driver.mlx      # Live script version
└── tool_driver.mlx        # Interactive tools
```

### 2.2 MATLAB Simulation Engine

#### Core Components

**A. Main Driver (`script_driver_slim.m`)**
- Entry point for running simulations
- Configures simulation parameters (physics, solver settings, thresholds)
- Bootstrap procedure: path setup, catalog initialization from seed shapes
- Delegates to `sim_explore_H0_quad_tree` for parameter exploration

**B. Parameter Space Explorer (`sim_explore_H0_quad_tree.m`)**
- Implements adaptive quad-tree scheduler for (H₀₁, H₀₂) space
- Main iteration loop:
  1. Query quad-tree for next parameter point to solve
  2. Generate physics-aware hash key (includes A, V, κₐ, κᵦ, κ_G, H₀₁, H₀₂)
  3. Check if solution already exists
  4. Find warm-start from nearby solved points or seed shapes
  5. Call BVP solver
  6. Persist results with atomic write (`.mat` files)
  7. Update catalog and cache
- Includes failure tracking with exponential backoff to avoid repeated failed solves

**C. BVP Solver (`solveAtParams.m`)**
- Multi-rung continuation solver for vesicle equilibrium shapes
- Handles homotopy from initial guess to target (H₀₁, H₀₂)
- Implements adaptive path strategies:
  - Straight continuation when path is smooth
  - Axis-aligned paths when crossing phase boundaries
  - Mesh coarsening on solver stalls
  - Multiple δ attempts with fallback logic
- Uses 6th-order BVP solver (`bvp6c`) for high accuracy
- Returns solution structure with geometry (r, z, angles) and metadata

**D. Physical Model (`BendV_Lag_EIGp_*_impl.m`)**
- Implements governing equations:
  - Bending energy: ∫ κ(H - H₀)² dA (per phase)
  - Optional Gaussian curvature term: ∫ κ_G K dA
  - Line tension at phase boundary
  - Area and volume constraints (Lagrange multipliers)
- Boundary conditions:
  - Pole regularity (south α-phase, north β-phase)
  - Neck interface conditions (continuity, force balance)
- Two versions: BC (boundary conditions) and DE (differential equations)

**E. Quad-Tree Scheduler (`processQuadtree.m`)**
- Adaptive refinement based on solution uniformity
- Corner-based cell logic: only subdivides cells where neighboring solutions differ significantly
- Uniformity tests:
  - Energy variation < ε_tol
  - Pressure variation < p_tol
  - Shape label consistency (shape classification metric)
- Configuration parameters: max depth (6), max cells (2000), tolerance thresholds

**F. Catalog System (`catalog_*.m`)**
- Persistent index of all solved points
- Structure: hash → {params, meta}
- Supports both seed shapes and computed solutions
- Enables quick lookup for warm-start selection

**G. Utilities (30+ helper functions)**
- `pickWarmStart.m`: Find best continuation starting point
- `initialGuessFromFile.m`: Load seed shapes from initial-shapes/
- `coarsen_mesh.m`: Adaptive mesh reduction for difficult solves
- `simpleDataHash.m`: Generate SHA-256 hashes for result identification
- `computePhaseScales.m`: Derive geometric scales from area fraction
- `uniformTest.m`: Assess if neighboring solutions are similar
- Various geometry helpers, diagnostics, and file I/O utilities

#### Key Parameters

**Physical Constants (sim.MP):**
- A: Area fraction of phase α (e.g., 0.50 = equal phases)
- V: Reduced volume (sphere = 1.0, typical vesicle ≈ 0.72)
- KA, KB: Bending rigidities of phases α and β (typically 1.0)
- KG: Gaussian bending modulus (typically 0.0)
- H₀₁, H₀₂: Spontaneous curvatures (sweep parameters)

**Solver Settings (sim.TH):**
- RelTol: 1e-6, AbsTol: 1e-8
- Delta: Continuation step size (0.01 default, adaptive)
- BCmax: Boundary condition residual threshold (1e-6)
- DEmax: Differential equation residual threshold (2e-1)

**Simulation Policy (sim.SP):**
- MaxIters: Maximum parameter points to attempt (15 in example)
- ModelVersion: "BVP-v3.1" (for reproducibility tracking)
- Verbose, LogToFile, SaveHomotopy: Diagnostic flags

### 2.3 Python ML Pipeline (Planned)

#### Architecture (from design documents)

**A. Data Schema (`feature-list-data-schema.yaml`)**
Defines comprehensive feature extraction approach:

1. **Spectral Descriptors** (20 Fourier modes)
   - Rotation/translation/scale-invariant shape representation
   - Complex coefficients from (r, z) profile: c_j = r_j + iz_j
   - Normalized via k=1 mode alignment

2. **Differential Geometry Statistics** (13 features)
   - Mean curvature H: mean, std, skew, kurtosis, L1, L2, pos/neg area
   - Gaussian curvature K: mean, std, L1, L2
   - Correlation between H and K along profile

3. **Topological & Landmark Features** (11 features)
   - Neck radius and tangent angle
   - Curvature extrema and locations
   - Number of lobes (from peak detection)
   - Invagination, tubulation, asymmetry scores

4. **Energetics & Constraints** (8 features)
   - Total energy and components (bending, line, Gaussian, surface)
   - Osmotic pressure
   - Total area and volume

5. **Control Parameters** (6 features)
   - Phase fraction (x₁), spontaneous curvatures (H₀₁, H₀₂)
   - Bending rigidity ratio, Gaussian flag, line tension

6. **Provenance Metadata**
   - Solution ID, run ID, branch ID (forward/backward/hysteresis)
   - Timestamps, optional human labels

**B. Analysis Pipeline (`pipeline-notebook.py`)**
Template notebook structure (not yet executed):

1. Load features from Parquet table (one row per solution)
2. Build feature matrix from geometry-only columns
3. Standardize features (Z-score normalization)
4. Dimensionality reduction:
   - PCA (keep 10 components, ~90% variance expected)
   - UMAP (2D for visualization)
5. Clustering:
   - K-means with silhouette-based k selection (sweep k=3..12)
   - Alternative: Gaussian Mixture Models or HDBSCAN
6. Visualizations:
   - PCA and UMAP scatter plots colored by cluster
   - Cluster medoids: representative shapes per family
   - Parameter plane heatmaps showing dominant cluster
7. Hysteresis analysis (if branch metadata available)
8. Supervised model training (optional):
   - Random Forest classifier for rapid shape labeling
   - Feature importance analysis

**C. Directory Layout (Planned)**
```
vesicle_ml/
├── metadata/
│   ├── schema.yaml
│   └── normalization.json
├── profiles/
│   └── <solution_id>.npz    # Resampled (r,z,s,H,K) arrays
├── features/
│   ├── features.parquet     # Tabular features
│   └── splits.json          # Train/val/test splits
├── figures/
└── params/                   # Optional raw parameter files
```

**Status:** No Python code has been executed yet. The `.yaml` schema and `.py` notebook are design documents only. The actual feature extraction, data conversion from MATLAB `.mat` to `.npz`/Parquet, and ML workflows are **not implemented**.

### 2.4 Code Quality Assessment

**Strengths:**
- Well-organized modular structure with clear separation of concerns
- Robust error handling with continuation strategies and fallback mechanisms
- Comprehensive hashing and caching system for reproducibility
- Atomic file writes prevent data corruption
- Extensive inline documentation in key solver functions
- Adaptive algorithms (quad-tree, mesh coarsening, continuation paths)

**Areas for Improvement:**
- No automated test suite (manual testing only)
- Limited code comments in some utility functions
- No continuous integration / automated validation
- Python ML code is template-only, not production-ready (see `notebooks/README.md` for status)
- Minimal version control history visible (limited commits in clone)

---

## 3. Next Steps to Make Progress

### Phase 1: Complete Simulation Sweep (Immediate - Weeks 1-2)

**Priority: High** - The simulation engine is ready; complete the parameter exploration.

1. **Run Extended Parameter Sweep**
   - Increase `MaxIters` in `script_driver_slim.m` from 15 to 1000-10000
   - Monitor simulation progress via `SimResults/OPfile.txt` logs
   - Estimated time: days to weeks depending on parameter range and hardware
   - Target: 500-2000 distinct solutions across (H₀₁, H₀₂) space

2. **Validate Solution Quality**
   - Implement diagnostic script to check:
     - Boundary condition residuals (< BCmax threshold)
     - Energy convergence
     - No NaN/Inf values in solutions
   - Create visualization tool for spot-checking random solution profiles
   - Document any systematic failure modes

3. **Expand Parameter Coverage (Optional)**
   - Vary area fraction A (currently fixed at 0.50)
   - Explore different reduced volumes V
   - Test non-zero Gaussian modulus κ_G
   - Each new global parameter set requires separate sweep + catalog

### Phase 2: Data Export & Feature Extraction (Weeks 3-4)

**Priority: High** - Bridge MATLAB simulation to Python ML pipeline.

4. **Implement MATLAB→Python Data Converter**
   - Create script `export_to_ml_format.m`:
     - Read all `.mat` files from `SimResults/hashed_results/`
     - Load catalog to get parameter metadata
     - For each solution:
       - Resample profile to N=512 points uniformly in arc-length
       - Apply Procrustes alignment (center, standardize orientation)
       - Compute curvatures H, K from geometry
       - Write to `.npz` file: {r, z, s, H, K, psi}
   - Write metadata JSON files (normalization stats, column descriptions)

5. **Feature Computation Pipeline**
   - Implement feature extraction in Python (`src/descriptors.py`):
     - Fourier descriptor calculation with invariance normalization
     - Differential geometry statistics (moments, norms)
     - Topology detection (peak finding for lobes, invagination heuristics)
     - Neck landmark identification
   - Process all profiles → single `features.parquet` table
   - Validate: check for NaNs, verify physical constraints (area, volume)

6. **Create Train/Val/Test Splits**
   - Random 70/15/15 split stratified by parameter region if possible
   - Generate 5-fold CV splits for robust evaluation
   - Document split methodology in `metadata/README.md`

### Phase 3: ML Analysis & Insights (Weeks 5-6)

**Priority: Medium-High** - Extract scientific value from simulation data.

7. **Execute Unsupervised Learning Pipeline**
   - Run `pipeline-notebook.py` (modify to actual data paths)
   - Generate PCA/UMAP embeddings
   - Perform k-means sweep, select optimal k (likely 4-8 clusters)
   - Visualize:
     - 2D scatter plots of shape space
     - Representative profiles (medoids) per cluster
     - Parameter plane heatmaps of cluster distribution

8. **Shape Family Characterization**
   - Manually inspect and label each cluster with descriptive names:
     - Prolate, oblate, budded, tubulated, multi-lobed, biconcave, etc.
   - Compute cluster statistics:
     - Average energies, curvatures, topological counts
     - Parameter ranges that produce each family
   - Identify sharp boundaries (first-order transitions) vs smooth crossovers

9. **Hysteresis & Multi-Stability Analysis**
   - If continuation branch data available:
     - Map forward vs backward sweep paths
     - Identify regions with multiple stable solutions
     - Quantify energy barriers between branches
   - Create phase diagrams showing regions of bistability

10. **Scientific Validation**
    - Compare findings to literature on vesicle shapes (e.g., Seifert 1997)
    - Check if known morphologies (stomatocytes, discocytes) are reproduced
    - Identify novel shapes or unexpected parameter regimes

### Phase 4: Automation & Tooling (Weeks 7-8)

**Priority: Medium** - Make the workflow reproducible and extensible.

11. **Dockerize Environment**
    - Create `Dockerfile` with MATLAB Runtime (or Octave if license-free needed)
    - Include Python 3.11+ with dependencies: numpy, scipy, scikit-learn, umap-learn, pandas, pyarrow, matplotlib
    - Document setup instructions in `README.md`

12. **Create End-to-End Workflow Scripts**
    - `run_simulation.sh`: Launch MATLAB driver, monitor completion
    - `export_data.sh`: Convert MATLAB results to ML format
    - `train_models.sh`: Execute full ML pipeline, save artifacts
    - Integrate with Make or similar build tool for dependency tracking

13. **Implement Automated Tests**
    - MATLAB unit tests:
      - Test uniformTest.m with synthetic data
      - Validate hash collision resistance
      - Test catalog append/load idempotency
    - Python unit tests:
      - Feature extraction on known shapes (sphere, cylinder)
      - Normalization correctness
      - Splitting reproducibility
    - Continuous integration via GitHub Actions (if feasible given MATLAB licensing)

14. **Performance Monitoring**
    - Add timing instrumentation to BVP solver
    - Log solver statistics: mesh size, iteration count, continuation steps
    - Create performance dashboard: average solve time vs parameter location

### Phase 5: Documentation & Dissemination (Weeks 9-10)

**Priority: Medium** - Share results with community and enable future work.

15. **Technical Documentation**
    - Write detailed `docs/SOLVER_GUIDE.md`:
      - Governing equations and discretization
      - Continuation strategies and when they're applied
      - Solver parameter tuning advice
    - Write `docs/ML_GUIDE.md`:
      - Feature design rationale
      - Interpretation of PCA components
      - How to extend with new features or models

16. **Research Artifacts**
    - Create Jupyter notebooks with polished visualizations:
      - Interactive shape browser (plotly or similar)
      - Parameter sweep animations
      - Cluster comparison galleries
    - Write methods section draft for paper/thesis
    - Prepare supplementary dataset with solution snapshots

17. **Code Publication**
    - Tag stable release (v1.0.0)
    - Archive on Zenodo with DOI
    - Write citation file (CITATION.cff)
    - Publish preprint (bioRxiv/arXiv) if results warrant

### Phase 6: Extensions & Future Directions (Months 3+)

**Priority: Low/Future** - Blue-sky ideas for continued research.

18. **Advanced ML Techniques**
    - Graph neural networks for vesicle topology
    - Variational autoencoders for generative shape design
    - Active learning to guide simulation to interesting regions
    - Transfer learning from synthetic to experimental microscopy data

19. **Physics Extensions**
    - Dynamics: Implement time-dependent Stokes flow solver
    - Non-equilibrium: Add membrane permeability, active pumping
    - Inhomogeneous fields: Spatially varying spontaneous curvature
    - Three-phase vesicles (ternary lipid mixtures)

20. **Experimental Integration**
    - Shape matching pipeline: fit model to microscopy images
    - Parameter inference: inverse problem (shape → physical constants)
    - Validation: compare predictions to GUV (giant unilamellar vesicle) experiments

---

## Appendix A: Key File Reference

### Essential MATLAB Files

| File | Lines | Purpose |
|------|-------|---------|
| `script_driver_slim.m` | 104 | Main entry point, configuration |
| `sim_explore_H0_quad_tree.m` | 263 | Parameter sweep driver |
| `solveAtParams_v2.m` | 227 | BVP solver with continuation |
| `processQuadtree.m` | ~70 | Adaptive scheduler |
| `BendV_Lag_EIGp_BC_impl.m` | ~90 | Boundary conditions |
| `BendV_Lag_EIGp_DE_impl.m` | ~120 | Differential equations |
| `pickWarmStart.m` | ~40 | Warm-start selection |
| `catalog_*.m` | ~100 | Persistent storage |

**Total MATLAB codebase:** ~90 files (including bvp6c-solver), ~3500 lines (excluding solver library)

### Configuration Files

- `.gitignore`: Excludes MATLAB autosave (`.asv`), OS files, build artifacts
- No `requirements.txt` yet (Python dependencies not pinned)
- No Dockerfile, Makefile, or CI config files

### Data Assets

- **initial-shapes/**: 4 seed shapes (A=0.50/0.60/0.75/0.90, V=0.72)
- **SimResults/**: Catalog with 8 initial seed entries + any solved points
- **docs/**: Research plan PDF outlines multi-phase vesicle theory and project motivation

---

## Appendix B: Estimated Resource Requirements

### Computational

- **Single BVP solve:** 0.1 - 60 seconds (highly variable)
  - Fast: Smooth continuation from good warm-start
  - Slow: Multiple fallback attempts, mesh coarsening, axis paths
- **Full parameter sweep:** 1-7 days for 1000 points (assuming 8-core workstation)
  - Parallelizable in principle (independent parameter points)
  - Current code is serial; MATLAB Parallel Computing Toolbox could parallelize
- **Feature extraction:** Minutes to hours for 1000 solutions (Python, single-threaded)
- **ML training:** Seconds to minutes (clustering/PCA/UMAP on ~50-dim x 1000-sample data)

### Storage

- **Per solution:** ~200-300 KB (MATLAB .mat with full trajectory)
- **1000 solutions:** ~200-300 MB
- **Exported ML data:**
  - Profiles (512 pts × 5 vars × 8 bytes): ~20 KB/solution → 20 MB total
  - Features table (Parquet): ~1 MB (highly compressible)
- **Total project footprint:** < 1 GB including code, docs, and initial results

### Human Effort

- **Phase 1-2 (Simulation & Export):** 2-3 weeks (mostly waiting for compute)
- **Phase 3 (ML Analysis):** 1-2 weeks of active work
- **Phase 4-5 (Infrastructure & Docs):** 2-3 weeks
- **Total to "first publishable results":** ~2 months at 50% FTE

---

## Appendix C: Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| BVP solver fails for large regions of parameter space | Medium | High | Implement more robust fallbacks; reduce refinement depth in problematic zones |
| Insufficient shape diversity for meaningful clustering | Low | Medium | Expand parameter ranges (vary A, V, κ_G); verified seed diversity exists |
| ML features lack discriminative power | Low | Medium | Schema includes 50+ diverse features; worst-case: use raw profile pixels with CNN |
| MATLAB licensing constraints for deployment | High | Low | Code runs fine on licensed workstation; consider Octave port for open-source future |
| Data management: catalog corruption or loss | Low | High | Implement backup script; catalog is append-only with atomic writes (already robust) |
| Python dependencies break over time | Medium | Low | Pin versions in requirements.txt; use virtual environment or container |

---

## Conclusions

The **vesicle-ml-shapes** project is scientifically well-motivated and technically sound. The MATLAB simulation engine is production-ready with sophisticated adaptive algorithms and robust data management. The next critical step is to **execute an extended parameter sweep** to generate sufficient data, followed by **feature extraction and ML analysis** to realize the project's scientific potential.

Key strengths include the mature BVP solver, quad-tree adaptive refinement, and comprehensive ML schema design. The main gap is the lack of implemented Python ML code—this is entirely scaffolding at present.

Recommended immediate actions:
1. Run large-scale simulation (increase MaxIters to 1000+)
2. Implement MATLAB→Python data export script
3. Execute feature extraction and initial clustering analysis

With focused effort over the next 2 months, this project can deliver novel insights into vesicle morphology and establish a powerful computational pipeline for membrane biophysics research.

---

**End of Audit Report**
