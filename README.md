# Multicarving for Correct Inference in High-Dimensional Logistic Regression Models
R code used for the Master Thesis of Jack McManus Heller and Viktoria Kirichenko.

The code is adapted from the work of Schultheiss et al. presented in https://arxiv.org/abs/2006.04613.

The code consists of the following:
 - In /inference/
   - Code used to perform multicarving as described in Schultheiss et al. This has been adapted to capture coefficients.

 - In /simulation_setups/
   - 6 files for simulations:
      - carving_binomial_simulation_betabinomial.R for multicarving with Y data generated from a beta-binomial distribution
      - carving_binomial_simulation_lambdaMin.R for multicarving with a cross-validated Lasso calculated with Lambda_min
      - carving_binomial_simulation.R for multicarving with simulated data
      - carving_cloglog.R for multicarving with the misspecified link of Cloglog
      - carving_loanapp.R for multicarving with real data
      - loanapp_dataCreation.R preprocess data to fit into carving_loanapp.R

 - In /plotting/: Various scripts for creating different plots for the thesis.
 - In /tabulation/: Various scripts for tabulating results and calculating performance metrics.

