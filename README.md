# TSSSr

Thompson Sampling with Smoothing Splines for Hyperparameter Optimization for R. This package provides the supplementary materials for the dissertation of Matthew Louis Davis University of Iowa (2023), as well as an implementation of the proposed procedure. 

The package is centered around the ``set_optim_thompson" function, which can be applied to ModelSpecification objects in MachineShop package for equipping them with the proposed procedure. 

Supplementary Materials and Data are provided as .zip files, and are referenced in the manuscript.

# Supplementary Materials

## RCode
- Boilerplate.Rmd is an rmarkdown containing example code for using TSSSr in practice.
- BasicImplementation.Rmd is an rmarkdown containing code used for the walkthrough presented.
- Cookbook.R provides code for preparing recipes (pre-processed data sets used for modelling).
- create_modspec.R is code for constructing model specification objects for the empirical evaluations.
- dataset_best.R are files for running the stopping rule proposed for each data set. 
- Demonstration.R provides code for reproducing the demonstrations in the overview of Chapter 2. 
- EmpiricalPlots.R provides code for plotting the empirical evaluations.
- FitLogisticWithThompsonXParam.R provides code for reproducing the computational complexity analyses for X = 2, 4, and 8 parameter logistic regression models.
- NBARecipe.R provides code for preparing recipes for NBA data.
- ProcessBest.R provides code for processing the stopping rule performances.
- ProcessCompComplexity.R provides code for reproducing the computational complexity analyses
- run_dataset.R are files for running the empirical evaluations vs. RGS, PSO, and BO.
- run_dataset_wahba.R are files for running the empirical evaluations vs. RGS, Greedy Wahba, and the adaptive knot modification.
- set_optim_thompson.R is the code used for equipping TSSS in the dissertation (included for reproducibility, although TSSSr is better documented).
- set_optim_thompson_customplots.R is the modified code used for the demonstration, in which plots are overlaid with previous grid search performance and points evaluated.
- set_optim_thompson_timed.R is the modified code used when conducting computational complexity analyses.
- set_optim_wahba.R is the code for equipping a Greedy Wahba optimizer.


## Results
- BestAverageAt100.xlsx contains the summary statistics on time and performance after 100 grid points were evaluated, for all empirical evaluations.
- ComputationalComplexity.xlsx provides the data for reproducing computational complexity analyses.
- Other .RData files are eponymous; the suffix "best" refers to stopping rule scenarios, the suffix "wahba" refers to scenario for which Greedy Wahba was evaluated. 
- StoppingRulePerformances.xlsx provides the eponymous performances when evaluating the stopping rule.

## TuneCNN
### CNN RCode:
- CNN_X.R for X = 1-8 provide the 8 implementations for running asynchronous parallel.
- CNN_Final.R for running the final CNN trained on the entire MNIST data set.

### CNN Results:
- FinalFitLog.txt provides information from the HPC system of Iowa outputted while running CNN_Final.R
- mnist_test_performance.RData contains data for the final performance metrics of the CNN.
- PredictedVsObserved.xlsx contains all of the out-of-sample predictions and actual values on the MNIST test set.
- TuningHistory.xlsx provides the exhaustive history of grid points evaluated and associated performance when tuning in asynchronous parallel.

# Data
- bs.RData provides box scores for NBA 2010-2022.
- game_logs_2010_2022.RData provides game logs for NBA 2010-2022.
- mnist_test_conv.RData provides a recipe of MNIST test data, in long form.
- mnist_train_conv.RData provides a recipe of MNIST training data, in long form.
- Recipes for ames, hnscc, and microbiome are eponymous.


