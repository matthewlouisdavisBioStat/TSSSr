# TSSSr

Thompson Sampling with Smoothing Splines for Hyperparameter Optimization for R. This package provides the supplementary materials for the dissertation of Matthew Louis Davis University of Iowa (2023), as well as an implementation of the proposed procedure. 

The package is centered around the ``set_optim_thompson" function, which can be applied to ModelSpecification objects in MachineShop package for equipping them with the proposed procedure. 

# Supplementary Materials

## RCode
- <data set>_best.R are files for running the stopping rule proposed for each data set. 
- run_<data set>.R are files for running the empirical evaluations vs. RGS, PSO, and BO.
- run_<data set>.R are files for running the empirical evaluations vs. RGS, Greedy Wahba, and the Adaptive Knot scheme.
- FitLogisticWithThompson< p >Param.R provides code for reproducing the computational complexity analyses for P = 2, 4, and 8 parameter logistic regression models.
- EmpiricalPlots.R provides code for plotting results.
- ProcessBest.R provides code for processing the stopping rule performances.
- Cookbook.R provides code for preparing recipes (pre-processed data sets used for modelling).
- NBARecipe.R provides code for preparing recipes for NBA data, which is more involved.

## Results

## TuneCNN

