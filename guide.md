# Vesicle Parameter Space Exploration — Beginner’s Guide

This project explores equilibrium shapes of **two-phase vesicles** across a 2D parameter space:
- **H0_1**: spontaneous curvature in phase α
- **H0_2**: spontaneous curvature in phase β

The goal is to build a **catalog** of solved shapes and refine sampling adaptively where shapes change rapidly.

---

## What the simulation does

At a high level:

1. Define a rectangular box in (H0_1, H0_2)
2. Use a **quadtree** to subdivide regions of parameter space
3. For each cell, solve the shape equations at its **corners**
4. Decide whether the cell is "uniform" or "mixed"
   - If uniform → stop refining that region
   - If mixed → subdivide into 4 smaller cells
5. Save every solution to disk, and maintain a searchable catalog

### Why quadtree refinement?
A uniform grid wastes time in regions where the shape barely changes.
The quadtree concentrates computation near transitions (budding, tubulation, etc.).

Text visualization of refinement:

```
Start: one big cell
+------------------+
| |
| |
| |
+------------------+

If corners disagree, subdivide:
+---------+--------+
| | |
| C1 | C2 |
+---------+--------+
| | |
| C4 | C3 |
+---------+--------+
```


---

## Project outputs (what gets saved)

All results live in `SimResults/`:

- `SimResults/hashed_results/<hash>.mat`
  - Contains `result` (the BVP solution) and `meta` (diagnostics + label)
- `SimResults/catalog.mat`
  - A table with rows: `hash, timestamp, entry`
  - `entry.params` stores physics + H0 values
  - `entry.meta` stores label/energy/pressure/diagnostics
- `SimResults/cache.mat`
  - Quadtree state and work queue (lets you resume runs)

---

## Quick start

### 1) Run sanity checks (recommended before any long run)

Create/run:
```matlab
script_sanity_checks
```

Optional physics smoke test:
```matlab
script_sanity_checks('doPhysicsSmoke',true)
```

This validates the pipeline without waiting days.

### 2) Running a real exploration

The main entry point is the driver script (example):
```matlab
warnState = bootstrap();
sim = sim_config;
sim_explore_H0_quad_tree(sim);
cleanup(warnState);
```

Key parameters you will usually tune:

Parameter space box

Set:

```matlab
sim.SP.H0Bounds = [-1 1; -1 1];   % [H0_1min H0_1max; H0_2min H0_2max]
```

Quadtree refinement limits
```matlab
sim.SP.QTmaxDepth = 7;     % deeper = finer sampling
sim.SP.QTmaxCells = 6000;  % memory/time limiter
```
Saving frequency
```matlab
sim.SP.CacheSaveEvery = 50;
```
### 3) How to resume a run

Just run the driver again.
The scheduler loads `SimResults/cache.mat` and `SimResults/catalog.mat` and continues.

### 4) Troubleshooting checklist
The run keeps retrying the same point

Check `SimResults/cache.mat` exists and updates

Increase backoff and/or retries:
```matlab
sim.SP.MaxRetries   = 3;
sim.SP.BackoffSteps = [0 1 2 4 8 16];
```
Solutions look suspicious (self-intersections, pinching)

Increase validity gates:
```matlab
sim.TH.rMin     = 1e-3;
sim.TH.rNeckMin = 5e-4;   % optional
```
Solver struggles / mesh explodes

Reduce H0 step size:
```matlab
sim.TH.minH0Step = 0.005;
```

Increase NMax:
```matlab
sim.TH.opts = bvpset(sim.TH.opts,'NMax',2500);
```

### 5) Interpreting labels (current simple classifier)

The current `labelFromSolution` is a lightweight heuristic (peak count in radius profile).
This is a placeholder for richer features later (PCA/UMAP/clustering).

### 6) Suggested workflow

 1. Smoke test small bounds and low depth

 2. Run exploration mode to get a coarse map

 3. Tighten gates and rerun selective regions for a high-trust catalog

 4. Export catalog entries for ML feature extraction

---

## Next Steps to get the catalog ML-ready
1) Extend `meta` to store a small set of cheap, stable features per solve:
- `rNeck`, `min(r)` interior, max curvature proxy, number of extrema, etc.
2) Store a *downsampled profile* (e.g., 200 points of `(r(s), z(s))`) for quick plotting/embedding.
3) Add an “export” tool that writes a single flat table/CSV of meta + features for Python.

---
