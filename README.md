# vesicle-ml-shapes

A computational framework for exploring equilibrium shapes of two-phase vesicles using numerical continuation methods and machine learning.

## Overview

This project systematically explores the parameter space of two-phase vesicles by solving boundary value problems (BVPs) for equilibrium shapes under bending energy minimization. It combines:

- **MATLAB Simulation Engine**: Mature BVP solver with adaptive quad-tree refinement for parameter space exploration
- **Python ML Pipeline**: (Planned) Machine learning tools for shape classification and analysis

## Quick Start

### Running Simulations

```matlab
% In MATLAB, run the main driver:
script_driver_slim
```

This will:
1. Bootstrap the simulation environment
2. Configure physical parameters and solver settings
3. Explore the (H₀₁, H₀₂) parameter space using adaptive quad-tree refinement
4. Save results to `SimResults/hashed_results/`

### Key Entry Points

- **`script_driver_slim.m`**: Main entry point for running simulations
- **`script_driver.mlx`**: MATLAB Live Script version with interactive controls
- **`tool_driver.mlx`**: Interactive tools for analysis

## Repository Structure

```
vesicle-ml-shapes/
├── bvp6c-solver/          # High-order (6th-order) BVP solver
│   └── bvp6c.m            # Modified MATLAB BVP solver
├── initial-shapes/         # Seed shapes for continuation
│   └── SIM_Node_*.mat     # Bootstrap shapes at different area fractions
├── SimResults/            # Simulation outputs
│   ├── hashed_results/    # Individual solution .mat files (SHA-256 named)
│   ├── catalog.mat        # Master index of all solutions
│   └── cache.mat          # Quad-tree state and work queue
├── src/
│   ├── sim_explore_H0_quad_tree.m  # Main parameter space explorer
│   ├── utils/             # Core utilities (29 functions)
│   │   ├── solveAtParams_v2.m      # Multi-stage continuation solver
│   │   ├── processQuadtree.m       # Adaptive refinement
│   │   ├── catalog_*.m             # Catalog management
│   │   └── ...                     # Other helpers
│   └── tools/             # Maintenance scripts
├── docs/                  # Research documentation
│   ├── Vesicle Sim Research Plan.pdf
│   └── explore_H0_quad_tree_README.mlx
├── notebooks/             # ML pipeline (templates - see notebooks/README.md)
│   ├── README.md                    # Implementation status
│   ├── feature-list-data-schema.yaml
│   └── pipeline-notebook.py
├── script_driver_slim.m   # Runnable simulation entry point
├── script_driver.mlx      # Live script version
└── PROJECT_AUDIT.md       # Detailed project documentation
```

## Requirements

- **MATLAB** R2020b or later (no special toolboxes required)
- **Python** 3.x (for future ML components - currently templates only)

## Documentation

- **`PROJECT_AUDIT.md`**: Comprehensive project documentation including architecture, next steps, and technical details
- **`docs/`**: Research background and technical guides
- **`notebooks/README.md`**: Status of ML pipeline (currently planning/template stage)

## Features

- **Adaptive Parameter Exploration**: Quad-tree refinement focuses computational effort on regions with high morphological variation
- **Robust Continuation**: Multi-stage solver with fallback strategies for difficult parameter regions
- **Persistent Catalog**: SHA-256 hashing ensures unique identification and reproducibility
- **Warm-Start Selection**: Intelligent selection of initial guesses from nearby solved points

## Scientific Background

Vesicles are closed lipid bilayer membranes that exhibit complex morphologies driven by:
- Bending energy minimization with spontaneous curvature
- Phase separation between lipid domains
- Area and volume constraints

This framework enables systematic exploration of vesicle shapes across parameter space to understand budding, tubulation, and other morphological transitions.

## Citation

If you use this code in your research, please cite:
```
[Citation information to be added]
```

## License

[License information to be added]