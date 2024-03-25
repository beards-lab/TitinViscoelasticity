# TitinViscoelasticity
Codes accompanying paper titled Mesoscopic Scale Analysis of Titin-Mediated Viscoelastic Passive Muscle Mechanics

## Data
This folder includes source data (*/Raw*) as well as processed, averaged data for model identification (*Avg...\*.csv*). 

## DataProcessing
This folder includes data averaging and correcting routines,

## Figures
Contains plotting routines, if not plotted already as a part of other scripts elsewhere.

## Model
Contains scripts to optimized and run titin viscoelastic model under various conditions.
* OptimizeCOmbined - runs and optimizes the model
* RunCombinedModel - evaluates a single pCa level of the model for given ramp velocities
* dXdT - ODE equations of the model, solved by the ode solver
* evalPowerFit - plots power-law fit of the data and model output.
