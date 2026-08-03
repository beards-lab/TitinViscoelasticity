# TitinViscoelasticity MATLAB Project - Codebase Analysis

## Executive Summary
This is a **viscoelastic mechanical model of the Titin protein** in cardiac muscle, implementing a state-evolution model with parameter optimization against experimental force-relaxation data. The project combines:
- **ODE-based passive model** (differential equations solved with ODE15s)
- **Parameter optimization** (using genetic algorithms, fminsearch, surrogateopt)
- **Data fitting** (experimental CSV data vs model outputs)
- **Multi-condition experiments** (6 calcium levels, 4 ramp rates)

---

## 1. KEY ENTRY POINT SCRIPTS

### Primary Execution Flow
| Script | Purpose | How to Run |
|--------|---------|-----------|
| [RunCombinedModel.m](Model/RunCombinedModel.m) | **Main model runner** - executes passive model simulation | Set `pCa`, `params`, `rampSet`, `drawPlots` variables, then run |
| [OptimizeCOmbined.m](Model/OptimizeCOmbined.m) | **Parameter optimization** - fits parameters to experimental data | Set `runOptim = true`, configure `modSel` and `pCas` |
| [SimRefolding.m](Model/SimRefolding.m) | **Refolding experiments** - tests dynamic relaxation protocols | Test script with predefined scenarios |
| [removeParamsRecursive.m](Model/removeParamsRecursive.m) | **Sensitivity analysis** - identifies minimal parameter set needed | Recursive search to remove parameters while maintaining fit |

### Data Processing Entry Points
| Script | Purpose |
|--------|---------|
| [FigFitDecayOverlay.m](DataProcessing/FigFitDecayOverlay.m) | Fit power-law decay to force-time data; generates publication figures |
| [AverageRamps.m](DataProcessing/AverageRamps.m) | Average experimental ramps across trials |
| [evalPowerFit.m](DataProcessing/evalPowerFit.m) | Cost function for power-law decay fitting |
| [PlotFigPeaks.m](Figures/PlotFigPeaks.m) | Compare peak forces across conditions |

### Core Model Functions
| Function | Role |
|----------|------|
| [dXdT.m](Model/dXdT.m) | **State evolution ODE** - heart of the model |
| [fitHill.m](Model/fitHill.m) | Hill function fitting for parameter trends across pCa |
| [ParamConversion.m](Model/ParamConversion.m) | Converts between old "mod" format and new "param" format |

---

## 2. TYPICAL WORKFLOW FOR RUNNING THE MODEL

### Basic Model Execution
```matlab
% Step 1: Set parameters
params = [5.19, 12.8, 4345, 2.37, 4e+04, 2.74, 8.658e+05, 5.797, 0.678, 0.165, 0.005381, 0.383];
pCa = 4.51;  % Calcium concentration (log scale, 4.51 = high Ca, 11 = relaxed)
rampSet = [1 2 3 4];  % Which ramp rates to run (1=100s, 4=0.1s)

% Step 2: Configure visualization options
drawPlots = true;
plotDetailedPlots = false;
plotInSeparateFigure = true;

% Step 3: Run model
RunCombinedModel;  % Loads data, simulates, plots results

% Outputs: Force{j}, Time{j}, Length{j}, states{j} for each ramp j
```

### Parameter Optimization Workflow
```matlab
% 1. Define subset of parameters to optimize
modSel = [3 4 5 6 7 8 9];  % Indices of parameters to fit

% 2. Set optimization bounds
init = params(modSel);
lb = 0.01 * init;   % Lower bounds
ub = 20 * init;     % Upper bounds

% 3. Define cost function
evalLin = @(optMods) evalCombined(optMods, params, modSel, [4.4]);

% 4. Optimize
options = optimset('Display','iter', 'TolFun', 1e-3, 'MaxIter', 100);
x = fminsearch(evalLin, init, options);
params(modSel) = x;  % Update parameters

% For more expensive optimization:
[x, fval] = surrogateopt(evalLin, lb, ub, optimoptions('surrogateopt', 'MaxTime', 6*3600));
```

### Data Processing & Fitting Workflow
```matlab
% 1. Load experimental data
load ../data/pca11data.mat;  % Loads Farr (forces), Tarr (times)

% 2. Fit power-law decay model
x = [3.7242, 0.2039, 4.8357];  % [amplitude, exponent, offset]
[cost, rampShift] = evalPowerFit(x, Farr, Tarr, true);

% 3. Generate publication figures
figure; plot results
exportgraphics(gcf, '../Figures/output.png', 'Resolution', 150);
```

### Multi-pCa Parameter Identification
```matlab
% Models for different calcium levels
modSet = [mod4; mod5_5; mod5_8; mod6; mod6_2; mod11];
pcax = [4.4, 5.5, 5.75, 6, 6.2, 11];

% Optimize specific parameters across pCa conditions
modSel = [3 4 5 6 7 8 9];
isolateRunCombinedModelAllpCas(optMods, modSel, modSet, pcax, pcaSel);
```

---

## 3. DATA FORMAT AND NAMING CONVENTIONS

### Experimental Data Files (CSV)
**Location**: `Data/` folder

**Naming convention**: `AvgpCa{pCa}_{ramp}.csv` or `AvgRelaxed_{ramp}.csv`

**File structure**:
```
Time,F,SD
0.001,0.5,0.1
0.002,0.7,0.15
...
```
- **Time**: Time (seconds, relative to ramp start)
- **F**: Force (kPa)
- **SD**: Standard deviation across trials

**Example files**:
- `AvgMava_pCa4.4_0.1s.csv` - pCa 4.4, 0.1s ramp
- `AvgRelaxedMavaSet_10s.csv` - Relaxed (pCa 11), 10s ramp
- `AvgpCa6.0_0.1s.csv` - pCa 6.0, 0.1s ramp

### Model Data Files (MAT)
**Example files**:
- `pca4.4modeldata.mat` - Pre-computed model outputs for pCa 4.4
- `pca11modeldata.mat` - Pre-computed model outputs for pCa 11

**Contents**: Structures containing `Farr`, `Tarr`, `states`, `cost`

### Ramp Rate Naming
| Array Index | Duration | Velocity |
|-------------|----------|----------|
| j=1 | 100s | Lmax/100 (slowest) |
| j=2 | 10s | Lmax/10 |
| j=3 | 1s | Lmax/1 |
| j=4 | 0.1s | Lmax/0.1 (fastest) |

### pCa Levels
- **pCa 4.4 (high Ca)**: Strongly activated muscle
- **pCa 5.5-6.2**: Intermediate calcium levels
- **pCa 11 (no Ca)**: Relaxed muscle (baseline)

---

## 4. CORE MODEL STRUCTURE AND EQUATIONS

### Model Architecture
```
dXdT.m (State evolution function)
  ├─ Input: current state x, parameters (kA, kD, kd, Fp, RU, RF, mu)
  ├─ Output: state derivatives dx/dt
  └─ Solves with ODE15s in RunCombinedModel.m

RunCombinedModel.m (Model orchestrator)
  ├─ Loads experimental data
  ├─ Sets up spatial/strain grid (Nx=15, Ng=10 states)
  ├─ Computes force basis functions Fp(s,n)
  ├─ Computes unfolding rates RU(s,n)
  └─ Integrates ODE in multiple phases (pre-ramp, ramp, decay)
```

### State Variables in dXdT.m
```matlab
pu(Nx, Ng+1)  % Probability of unfolded states (unattached)
pa(Nx, Ng+1)  % Probability of unfolded states (attached, Ca-dependent)
L              % Sarcomere length
```

- **Nx = 15**: Spatial discretization points along half-sarcomere
- **Ng = 10**: Number of folding states in Titin molecule
- Each state represents degree of unfolding (n=0: fully folded, n=Ng: fully unfolded)

### Key Equations
1. **Proximal chain force**: `Fp = kp * max(0, s - slack)^np`
   - Force increases with strain and depends on folding state
   
2. **Unfolding rate**: `RU = alphaU * (s/L_0)^nU * (Ng - n)`
   - State-dependent: faster unfolding at higher strains
   
3. **Distal chain force**: `Fd = kd * max(0, L - s)^nd`
   - Exponential dependence on stretch
   
4. **Total force**: `Force = sum(Fd .* pu) + attachment force`

### Numerical Integration Approach
```matlab
% Three-phase protocol
times = [-100, 0; 0, Tend_ramp; Tend_ramp, Tend_ramp + 200];
velocities = {0, V, 0};  % Pre-ramp, ramp, decay

% ODE solver settings
opts = odeset('RelTol', 1e-2, 'AbsTol', 1e-2);
[t, x] = ode15s(@dXdT, times, x0, opts, Nx, Ng, ds, kA, kD, kd, ...);
```

---

## 5. CONFIGURATION PATTERNS

### Parameter Definition
**Core 12 parameters** (from RunCombinedModel.m):
```matlab
paramNames = {'\F_{ss}', 'n_ss', 'k_p', 'n_p', 'k_d', 'n_d', ...
              '\alpha_U', 'n_U', '\mu', '\delta_U', 'k_{A}', 'k_{D}'};

Fss  = params(1);   % Steady-state force level
n_ss = params(2);   % Steady-state exponent
kp   = params(3);   % Proximal chain stiffness
np   = params(4);   % Proximal chain exponent
kd   = params(5);   % Distal chain stiffness
nd   = params(6);   % Distal chain exponent
alphaU = params(7); % Unfolding rate constant
nU   = params(8);   % Unfolding rate exponent
mu   = params(9);   % Viscous damping
delU = params(10);  % Unfolding distance
kA   = params(11);  % PEVK attachment rate
kD   = params(12);  % PEVK detachment rate
```

### Hard-Coded Spatial/Temporal Parameters
```matlab
% Numerical grid
Lmax = 1.175 - 0.95;  % Half-sarcomere stretch (225 nm)
Nx = 15;              % Spatial discretization
Ng = 10;              % Number of Titin folding states
ds = Lmax / (Nx - 1); % Spatial step

% Ramp rates
rds = fliplr([0.1, 1, 10, 100]);  % Ramp durations (seconds)
Vlist = Lmax ./ rds;              % Corresponding velocities
```

### Configuration Variable Patterns
Scripts use **optional variable checks** for flexibility:
```matlab
% Example pattern in RunCombinedModel.m
if ~exist('drawPlots', 'var')
    drawPlots = true;  % Default value
end
if ~exist('pCa', 'var')
    pCa = 11;  % Default to relaxed
end
if ~exist('params', 'var')
    % Load default parameter set
    params = [5.19, 12.8, ...];
end
```

### MAT File Storage
- **Optimization results**: `fmisrch_Ca.mat`, `surropt_All.mat`
- **Color schemes**: `SoHot.mat`, `SoCool.mat` (for figures)
- **Pre-fit data**: `pca11modeldata.mat`, `pca4.4modeldata.mat`

---

## 6. DOCUMENTATION AND INLINE COMMENTS

### Code Documentation Style
- **Minimal formal documentation** (no function headers or comments)
- **Inline comments** explain algorithm sections but are sparse
- **Variable naming** mostly descriptive but uses physics conventions
  - `pu` = unfolded probability, unattached
  - `pa` = unfolded probability, attached
  - `RU` = unfolding rate
  - `RF` = refolding rate

### Learning Resources Within Codebase
1. **OptimizeCOmbined.m**: Shows parameter optimization patterns
2. **SimRefolding.m**: Examples of running different simulation protocols
3. **ParamConversion.m**: Maps old/new parameter naming conventions
4. **fitHill.m**: Hill function fitting methodology
5. **dXdT.m**: Core mathematics (PDE discretization, upwind differencing)

### Missing Documentation
- ❌ No README explaining model theory
- ❌ No parameter fitting guide or tutorial
- ❌ No docstrings for functions
- ❌ No workflow diagram

---

## 7. COMMON DEVELOPMENT TASKS

### Task: Run Base Model for Specific Condition
```matlab
clear; close all;
params = [5.19, 12.8, 4345, 2.37, 4e+04, 2.74, 8.658e+05, 5.797, 0.678, 0.165, 0.005381, 0.383];
pCa = 4.51;
rampSet = [4];  % Only fastest ramp
drawPlots = true;
RunCombinedModel;
```

### Task: Optimize Parameters Against Single pCa
```matlab
modSel = [3 4 5 6 7 8 9];  % Parameters to optimize
init = params(modSel);
evalLin = @(optMods) evalCombined(optMods, params, modSel, pCa);
options = optimset('Display','iter', 'TolFun', 1e-3, 'MaxIter', 100);
x = fminsearch(evalLin, init, options);
params(modSel) = x;
```

### Task: Compare Against Experimental Data
```matlab
% Run model
pCa = 4.51;
RunCombinedModel;  % Produces Force{j}, Time{j}

% Compare with experimental data
load('../Data/AvgMava_pCa4.4_0.1s.csv');
expt_data = readtable('../Data/AvgMava_pCa4.4_0.1s.csv');

figure; 
plot(Time{4}, Force{4}, 'b-', 'LineWidth', 2); hold on;
errorbar(expt_data.Time, expt_data.F, expt_data.SD, 'ko');
```

### Task: Sensitivity Analysis
See **OptimizeCOmbined.m** (~line 520):
```matlab
perturbation = 0.01;  % 1% perturbation
for i = modSel
    perturbedParams = params;
    perturbedParams(i) = params(i) * (1 + perturbation);
    
    baseOutput = evalCombined([], params, [], [4.4]);
    perturbedOutput = evalCombined([], perturbedParams, [], [4.4]);
    
    sensitivity(i) = (perturbedOutput - baseOutput) / (params(i) * perturbation);
end
```

### Task: Generate Publication Figure
```matlab
f = figure(1); clf;
f.Position = [300 200 7.2*96 7.2*96/1.5];  % Set size for 2-col journal

% Plot model vs data
for i_pca = 1:length(pCas)
    params = paramSet(i_pca, :);
    pCa = pcas(i_pca);
    RunCombinedModel;
    % ... plotting code
end

exportgraphics(f, '../Figures/ModelResults.png', 'Resolution', 150);
exportgraphics(f, '../Figures/ModelResults.eps');
```

---

## 8. KEY PATTERNS AND CONVENTIONS FOR AI AGENTS

### Pattern 1: Script-Based Configuration
- **No function arguments** for main entry points
- **Uses workspace variables** for configuration
- **Check-and-default pattern** for optional variables
- ⚠️ **AI Challenge**: Must carefully manage workspace state

### Pattern 2: Parameter Vectorization
- **Long parameter vectors** (12-23 elements)
- **Index-based access** rather than named struct fields
- **Requires tracking index meanings** (see ParamConversion.m)
- **Multiple "formats"**: old "mod" format vs. new "param" format

### Pattern 3: Multi-Condition Nested Loops
```matlab
for i_pca = 1:length(pcas)
    for j_ramp = rampSet
        % Compute for this condition
    end
end
```
- pCa levels and ramp rates are typical loop dimensions
- Data, models, and results indexed by both

### Pattern 4: Callback-Style Cost Functions
```matlab
evalLin = @(optMods) evalCombined(optMods, params, modSel, pCas);
x = fminsearch(evalLin, init, options);
```
- Cost functions defined as anonymous function handles
- Standard MATLAB optimization interface

### Pattern 5: Sparse Inline Plotting
- Heavy use of `nexttile()` and tiledlayout
- Figure positioning: `figure.Position = [x y width height]`
- Colormap customization stored in MAT files

---

## 9. RECOMMENDATIONS FOR AI AGENT APPROACH

### ✅ DO
1. **Always call RunCombinedModel as a script** (don't convert to function without testing)
2. **Document parameter index meanings** when modifying params
3. **Save parameter sets to MAT files** for reproducibility
4. **Use cost functions (evalCombined) for optimization** rather than direct model calls
5. **Respect the three-phase ODE integration** (pre-ramp, ramp, decay)
6. **Check for data existence** before calling RunCombinedModel
7. **Test single-ramp scenarios** before multi-ramp
8. **Maintain workspace state carefully** - scripts depend on prior assignments

### ❌ AVOID
1. Converting scripts to functions without understanding workspace dependencies
2. Modifying hard-coded spatial parameters (Nx, Ng, Lmax) without validation
3. Mixing old "mod" format with new "param" format without conversion
4. Running optimization without defining `runOptim = true`
5. Changing ODE tolerances (AbsTol, RelTol) without benchmarking
6. Assuming parameter sets are universal across pCa conditions
7. Direct array indexing without consulting ParamConversion.m

### 🔧 TROUBLESHOOTING CHECKLIST
| Issue | Check |
|-------|-------|
| Model doesn't run | Verify data files exist in `../Data/` |
| Wrong output | Confirm `pCa` and `params` are set before RunCombinedModel |
| Slow optimization | Check `runOptim=false` to disable GA, use fminsearch |
| Cost = Inf | Verify `params > 0` and all data loaded |
| Memory error | Reduce `Nx` or `rds` length; use fewer ramp rates |

---

## 10. PROJECT STATISTICS

| Metric | Value |
|--------|-------|
| **Total MATLAB files** | 15+ scripts + functions |
| **Total parameters** | 12 core, up to 23 in optimization |
| **Experimental conditions** | 6 pCa levels × 4 ramp rates = 24 combinations |
| **Typical run time** | 5-10s per single ramp (ODE solve) |
| **Total code lines** | ~3000 (heavily commented in places, sparse in others) |
| **Main solver** | ODE15s (stiff ODE solver) |
| **Optimization methods** | fminsearch, surrogateopt, genetic algorithm |

---

## 11. ENTRY POINT DECISION TREE

```
START
 ↓
Want to run model?
 ├─YES→ Set params, pCa, rampSet → RunCombinedModel
 └─NO→ Want to optimize?
      ├─YES→ Define modSel → OptimizeCOmbined (set runOptim=true)
      └─NO→ Want to analyze data?
           ├─YES→ Load data → evalPowerFit or FigFitDecayOverlay
           └─NO→ Want sensitivity test?
                ├─YES→ Use removeParamsRecursive
                └─NO→ Check SimRefolding for examples
```

---

## SUMMARY TABLE: Quick Reference

| Aspect | Details |
|--------|---------|
| **Main Model** | Passive Titin viscoelasticity (ODE-based) |
| **Core Function** | `dXdT.m` (state evolution) |
| **Solver** | ODE15s (2nd order stiff ODE) |
| **Entry Point** | `RunCombinedModel.m` (script, not function) |
| **Workflow** | Load data → set params → run ODE → compare results |
| **Parameters** | 12 core (extensible to 23) |
| **Conditions** | 6 pCa × 4 ramp rates |
| **Optimization** | fminsearch, surrogateopt, genetic algorithm |
| **Key Files** | OptimizeCOmbined, SimRefolding, ParamConversion |
| **Data Format** | CSV (experimental), MAT (pre-computed) |
| **Visualization** | MATLAB figures exported as PNG/EPS |
