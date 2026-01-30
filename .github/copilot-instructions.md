# GitHub Copilot Instructions for vesicle-ml-shapes

## Project Overview

This repository contains a computational framework for exploring equilibrium shapes of two-phase vesicles using numerical continuation methods and machine learning. The project combines:

1. **MATLAB Simulation Engine**: A mature BVP (Boundary Value Problem) solver that systematically explores the parameter space of spontaneous curvatures (H₀₁, H₀₂) using adaptive quad-tree refinement
2. **Python ML Pipeline**: (In development) Machine learning tools for shape classification and analysis

The simulation framework solves for vesicle equilibrium shapes under bending energy minimization with area and volume constraints, generating a comprehensive catalog of morphologies across parameter space.

## Project Structure

- **`bvp6c-solver/`**: High-order (6th-order) BVP solver, adapted from MATLAB's standard solver
- **`InitialShapes/`**: Seed shapes (`.mat` files) used as starting points for continuation
- **`SimResults/`**: Simulation outputs including:
  - `solutions/`: Individual solution `.mat` files (hash-named for uniqueness)
  - `catalog.mat`: Master index of all solved parameter points
  - `cache.mat`: Quad-tree state and work queue for parameter exploration
- **`src/`**: Core simulation code
  - `utils/`: 30+ helper functions including solvers, catalog management, geometry utilities
  - `tools/`: Maintenance and migration scripts
  - Main drivers for parameter exploration
- **`notebooks/`**: Python ML pipeline design (feature schemas, analysis templates)
- **`docs/`**: Research documentation and technical guides
- **`retired/`**: Deprecated code versions (do not modify)

## Language and Environment

- **Primary Language**: MATLAB (R2020b or later recommended)
- **Secondary Language**: Python 3.x (for ML components in `notebooks/`)
- **Key Dependencies**:
  - MATLAB base installation (no special toolboxes required for core functionality)
  - Python: numpy, pandas, scikit-learn (for ML pipeline when implemented)

## Key Entry Points

1. **`script_driver_slim.m`**: Main entry point for running simulations
   - Configures physical parameters and solver settings
   - Initializes catalog and bootstraps from seed shapes
   - Delegates to parameter space explorer

2. **`src/sim_explore_H0_quad_tree.m`**: Adaptive parameter exploration
   - Implements quad-tree scheduler for (H₀₁, H₀₂) space
   - Handles solution hashing, catalog lookups, and warm-start selection
   - Main iteration loop with automatic failure tracking

3. **`src/utils/solveAtParams_v2.m`**: Multi-stage continuation solver
   - Homotopy-based continuation from initial guess to target parameters
   - Adaptive path strategies for crossing phase boundaries
   - Mesh coarsening and fallback logic for difficult solves

## Coding Conventions

### MATLAB Code

- **Function Names**: Use camelCase (e.g., `pickWarmStart`, `coarsen_mesh`)
- **Variable Names**: Use descriptive names; common conventions:
  - `sim`: Structure containing simulation configuration (`.MP` for physical params, `.TH` for solver thresholds)
  - `sol`: Solution structure from BVP solver
  - `catalog`: Master index of solved points
  - `cache`: Quad-tree state and work queue
- **Comments**: Use `%` for single-line comments; include header blocks for functions with purpose, inputs, outputs
- **File I/O**: All simulation results use `.mat` format with hash-based naming for uniqueness
- **Error Handling**: Use try-catch blocks for solver calls; track failures in cache to avoid repeated attempts

### Python Code

- **Style**: Follow PEP 8 conventions
- **Type Hints**: Use where appropriate for clarity
- **Imports**: Group standard library, third-party, and local imports separately

## Important Constraints and Guidelines

### Do Not Modify

- **`bvp6c-solver/`**: This is a specialized 6th-order BVP solver adapted from MATLAB. Changes here could break numerical accuracy.
- **`retired/`**: Contains deprecated code for historical reference only
- **Physical Model Files** (`BendV_Lag_EIGp_*_impl.m`): These encode the governing equations for vesicle mechanics. Modifications require deep understanding of the underlying physics.

### When Making Changes

- **Preserve Hash Consistency**: The hash generation in `simpleDataHash.m` ensures unique identification of parameter combinations. Changes to hash logic will invalidate the entire catalog.
- **Catalog Integrity**: Always use the catalog utility functions (`catalog_append`, `catalog_load`, etc.) rather than directly modifying `catalog.mat`
- **Warm-Start Logic**: The continuation solver relies heavily on finding good initial guesses. Changes to `pickWarmStart.m` can significantly impact solve success rates.
- **Quad-Tree Parameters**: The refinement thresholds in `processQuadtree.m` control exploration density. Changes affect computational cost and result completeness.

### Testing and Validation

- **No Automated Test Suite**: Currently no formal test infrastructure exists
- **Manual Testing**: 
  - For solver changes: Run small parameter sweeps and verify energy convergence
  - For catalog changes: Check that results can be loaded and warm-starts selected correctly
  - For quad-tree changes: Monitor subdivision patterns and uniformity metrics
- **Validation Approaches**:
  - Compare energy values and pressure checks from `bc_diagnostics.m`
  - Verify solution continuity across parameter space
  - Check that boundary conditions are satisfied within tolerance

## Common Development Tasks

### Adding New Utility Functions

Place in `src/utils/` with clear function header:
```matlab
function output = myNewFunction(input1, input2)
% MYNEWFUNCTION Brief description of purpose
%
% Inputs:
%   input1 - Description of first input
%   input2 - Description of second input
%
% Outputs:
%   output - Description of output
%
% Example:
%   result = myNewFunction(x, y);
```

### Modifying Solver Parameters

Edit the configuration section in `script_driver_slim.m`:
- Physical parameters: `sim.MP` structure (A, V, KA, KB, KG, H01, H02)
- Solver settings: `sim.TH` structure (RelTol, AbsTol, Delta, BCmax)

### Extending the ML Pipeline

Work in `notebooks/`:
- Follow the schema defined in `feature-list-data-schema.yaml`
- Extract features from `.mat` files in `SimResults/solutions/`
- Use catalog.mat for metadata (parameters, convergence info)

## Performance Considerations

- **BVP Solves**: Most computationally expensive operation (seconds to minutes per solve)
- **Quad-Tree Overhead**: Negligible compared to BVP solves
- **Catalog Lookups**: Linear scan through hash keys; consider optimization if catalog exceeds 10,000+ entries
- **Warm-Start Selection**: Distance calculations scale with catalog size; current implementation acceptable for expected dataset sizes

## Known Issues and Limitations

- Solver can fail in regions with complex bifurcations or multi-stability
- Exponential backoff for failed parameters helps but doesn't guarantee coverage
- Phase boundary crossings require careful continuation paths
- No automatic mesh adaptation in current solver (manual coarsening only)

## Scientific Background

Vesicles are closed lipid bilayer membranes that can exhibit complex morphologies driven by:
- **Bending Energy**: Resistance to curvature deviations from spontaneous curvature H₀
- **Phase Separation**: Two lipid phases with different H₀ values create interfaces
- **Constraints**: Fixed area and enclosed volume

The equilibrium shape minimizes total energy subject to constraints, leading to rich behavior including budding, tubulation, and multi-lobed structures.

## Getting Help

- **Documentation**: See `PROJECT_AUDIT.md` for detailed architecture overview
- **Research Context**: Check `docs/` for technical background and quad-tree exploration strategy
- **Code Comments**: Most utility functions have inline documentation

## Suggested Workflow for Code Changes

1. **Understand the Context**: Read relevant sections of `PROJECT_AUDIT.md` and function headers
2. **Test Locally**: Run small-scale tests (e.g., single parameter point) before full sweeps
3. **Verify Results**: Check solution quality using `bc_diagnostics.m` and energy convergence
4. **Preserve Compatibility**: Ensure changes don't break catalog loading or warm-start selection
5. **Document Changes**: Update comments and documentation as needed

## Future Development Priorities

1. Complete Python ML pipeline implementation
2. Add automated tests for critical utility functions
3. Implement parallel solver execution for parameter sweeps
4. Enhance visualization tools for shape analysis
5. Develop web interface for exploring the shape catalog
