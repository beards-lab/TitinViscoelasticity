# TitinViscoelasticity - AI Agent Guidelines

This workspace is a MATLAB-based simulation and optimization project for Titin-mediated viscoelastic passive muscle mechanics.

## AI Agent Guidance

- Focus on MATLAB scripts and workspace-variable configuration, not building a new CLI or web app.
- Use `Model/RunCombinedModel.m` and `Model/OptimizeCOmbined.m` as the primary entry points for model execution and fitting.
- Place new analysis or plotting scripts in `Model/`, `DataProcessing/`, or `Figures/` and add a short header comment explaining the entry point.
- Preserve the existing workflow: scripts set `pCa`, `rampSet`, and parameter vectors in the MATLAB workspace.
- When modifying model code, avoid changing ODE tolerances or parameter indexing without validating with small test runs.

## Quick Start

- **Main model runner**: [Model/RunCombinedModel.m](Model/RunCombinedModel.m) – Set `pCa` and `rampSet` in workspace, then run
- **Parameter optimization**: [Model/OptimizeCOmbined.m](Model/OptimizeCOmbined.m) – Fits model to experimental data
- **Core model equations**: [Model/dXdT.m](Model/dXdT.m) – ODE state evolution (Titin unfolding dynamics)
- **Experimental data**: [Data/](Data/) – CSV files with force/time measurements at various pCa levels

## Key Concepts

### Configuration (Script-Based, Not Functions)
Scripts use **workspace variables** for configuration—don't pass parameters as function arguments:
```matlab
pCa = 4.4;  % Calcium concentration (or 11 for relaxed)
rampSet = [0.1, 1, 10, 100];  % Ramp velocities (seconds)
```

### Model State Space
- **State variables**: 2 × 15 spatial points × 10 folding states = 300 ODE equations
- **Solver**: ODE15s (stiff) with 1e-2 relative/absolute tolerances
- **Integration phases**: pre-ramp → ramp → decay
- **Core parameters**: 12 fitted values (Fss, n_ss, kp, np, kd, nd, alphaU, nU, mu, delU, kA, kD)

### Data Format
- **Experimental**: `AvgMava_pCa{X.X}_{TIME}s.csv` (force vs. time for one condition)
- **Model outputs**: `.mat` files storing optimized parameters and results
- **Conditions**: 6 pCa levels × 4 ramp rates = 24 scenarios

## Common Tasks

| Task | Entry Point | Notes |
|------|-------------|-------|
| Run model at single condition | `RunCombinedModel.m` | Set `pCa`, `rampSet` first |
| Fit parameters to data | `OptimizeCOmbined.m` | Uses fminsearch/surrogateopt; may take 10-30 min |
| Evaluate fitness | `evalCombined.m` | Compares model vs. experimental data |
| Parameter sensitivity | `removeParamsRecursive.m` | Tests which parameters matter most |
| Test dynamic scenarios | `SimRefolding.m` | Custom force/time protocols |
| Plot results | [Figures/](Figures/) folder | Various plotting scripts |
| Data processing | [DataProcessing/](DataProcessing/) folder | Averaging and calibration routines |

## Important Conventions

⚠️ **Before modifying model code:**
- Verify ODE tolerance changes don't affect output (they can strongly affect results)
- Check parameter index mappings via `ParamConversion.m` when working with vectors
- Test on a small dataset first (e.g., single pCa/ramp condition)

🔧 **When adding features:**
- Use workspace variables for config (don't add function parameters)
- Place new scripts in `Model/`, `DataProcessing/`, or `Figures/` as appropriate
- Document which `.m` files are entry points vs. utility functions

📊 **When running optimization:**
- Expect 10-30 minutes for parameter fitting
- Monitor cost function convergence
- Save `.mat` files with meaningful names (e.g., `fminsearch_pCa4.4_v2.mat`)

## Documentation

For deeper understanding of the project, see:
- [README.md](README.md) – Project overview
- [CODEBASE_ANALYSIS.md](CODEBASE_ANALYSIS.md) – Detailed architecture, workflows, and troubleshooting
- Individual script headers for implementation details

## Project License

See [LICENSE](LICENSE) for usage terms.
