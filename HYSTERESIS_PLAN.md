# Hysteresis Support — Multi-Solution Design

## Problem

The quadtree simulation explores equilibrium shapes of two-phase vesicles over the
`(H0_1, H0_2)` spontaneous-curvature parameter plane.  For a given physics
configuration `(A, V, KA, KB, KG)`, the system can admit **multiple equilibrium
shapes** at the same `(H0_1, H0_2)` — a phenomenon known as **hysteresis** or
**multistability**.

### Audit Finding

Before this change, the system **could not** capture hysteresis because:

1. **Hash = physics identity, not solution identity.**  The SHA-256 hash was computed
   from `(model_version, H0_1, H0_2, A, V, KA, KB, KG)`.  This means each
   `(H0_1, H0_2)` at fixed physics maps to exactly **one** hash.
2. **File-existence skip.**  On each iteration the driver checks
   `exist(fmat, 'file')` and skips re-solving if a result already exists on disk.
3. **Catalog merge.**  `catalog_append` merges entries with the same hash,
   overwriting metadata rather than adding a new row.
4. **Single-entry lookup.**  The quadtree scheduler's `lookup_in_catalog` returns
   only one entry per `(H0_1, H0_2)` (the last match).

Together, these guarantees mean that **only one solution per `(H0_1, H0_2)` could
ever be stored**, regardless of which continuation path was used to reach it.

### What Makes Hysteresis Possible Physically

Different warm-start solutions (i.e., different continuation paths through the
`(H0_1, H0_2)` plane) can converge to different local energy minima at the same
target point.  The BVP solver performs Newton-type continuation; which basin of
attraction it lands in depends on the initial guess.

---

## Solution: Branch Tags

### Concept

A **branch tag** is a user-defined string label (e.g., `"upper"`, `"lower"`,
`"path-A"`) that distinguishes independent solution families.  When set, the
branch tag becomes part of the hash key, so the same `(H0_1, H0_2)` produces
**different hashes** for different branches.

### Backward Compatibility

When `BranchTag` is absent or empty (the default), the hash is identical to the
original — no existing catalogs or result files are invalidated.

### How to Use

To explore a second solution branch, run the driver with a different branch tag:

```matlab
sim = sim_config();
sim.SP.BranchTag = "upper";   % or "lower", "path-A", etc.
sim_explore_H0_quad_tree(sim);
```

Each branch-tag run shares the same `cache.mat` file, but the quadtree state
(`cache.QT`) and failure registry are automatically reset when the branch tag
changes from the previous run.  This ensures each branch explores the
`(H0_1, H0_2)` plane independently.  The catalog will contain multiple
entries for the same `(H0_1, H0_2)`, distinguished by their hash (which includes
the branch tag) and the `meta.branch_tag` field.

### What Changed

| File | Change |
|------|--------|
| `src/utils/simpleDataHash.m` | Include `branch_tag` in hash key when non-empty |
| `src/quadtree/sim_explore_H0_quad_tree.m` | Read `SP.BranchTag`, add to hash key and meta; branch-aware warm-start |
| `src/quadtree/processQuadtree.m` | Branch-aware `corner_hash`, `lookup_in_catalog`, `refresh_corners_from_catalog` |
| `src/script_driver_slim.m` | Document `SP.BranchTag` option in `sim_config()` |

### Architecture

```
                    BranchTag = ""               BranchTag = "upper"
                    ──────────────               ───────────────────
Hash key:           {mv, H0_1, H0_2, A, V,      {mv, H0_1, H0_2, A, V,
                     KA, KB, KG}                  KA, KB, KG, branch_tag}

Result file:        <hash-A>.mat                 <hash-B>.mat
Catalog entry:      hash-A → entry               hash-B → entry
                    (meta.branch_tag absent)      (meta.branch_tag = "upper")
Warm-start:         any solved neighbor           prefer same-branch neighbor
Quadtree lookup:    matches all entries           matches only "upper" entries
```

### Multi-Branch Workflow

```
% Run 1: default branch (no tag)
sim = sim_config();
sim_explore_H0_quad_tree(sim);

% Run 2: explore "upper" branch with different warm-start strategy
sim2 = sim_config();
sim2.SP.BranchTag = "upper";
% Optionally customize warm-start or solver settings to bias toward upper branch
sim_explore_H0_quad_tree(sim2);

% Run 3: explore "lower" branch
sim3 = sim_config();
sim3.SP.BranchTag = "lower";
sim_explore_H0_quad_tree(sim3);
```

All three runs share the same `catalog.mat` and `hashed_results/` directory.
Solutions are distinguishable by their `meta.branch_tag` field.

---

## Future Work

- **Automatic bifurcation detection**: Identify regions where the quadtree's
  uniform test fails due to label disagreements and automatically launch
  secondary branch explorations.
- **Branch-aware visualization**: Extend `plot_phase_map.m` to overlay
  multiple branches with different colors/markers.
- **Energy comparison**: Post-processing script to compare energies across
  branches at the same `(H0_1, H0_2)` to identify ground states vs. metastable
  states.
- **Warm-start seeding across branches**: Allow one branch to seed another's
  warm-start by loading solutions from a different branch and perturbing them.
