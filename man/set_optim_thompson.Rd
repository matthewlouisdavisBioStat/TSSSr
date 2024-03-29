% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_optim_thompson.R
\name{set_optim_thompson}
\alias{set_optim_thompson}
\title{Equip ModelSpecification objects in MachineShop with a Thompson sampler with smoothing splines for tuning}
\usage{
set_optim_thompson(
  object,
  times = 10,
  initial_points = 10,
  degree = 3,
  plot_predicted_performance = T,
  max_grid_rows = 2500000,
  max_knots = Inf,
  obs_per_spline_coefficient = NULL,
  random = TRUE,
  use_transf = TRUE,
  log_outcome = NULL,
  logit_outcome = NULL,
  save_model_fit = FALSE,
  include_numeric_interactions = T,
  include_high_order_interactions = F,
  include_quadratic_interactions = F,
  qualify_numeric_numeric_interactions = F,
  save_model_fit_lable = "thompson_model_fit.RData",
  generic_tolerance = sqrt(.Machine$double.eps),
  cardinal_knots = T,
  latin_grid = T,
  enforce_positive_definiteness = F,
  acquisition_function = NULL,
  print_colnames = T,
  stopping_rule_total = NULL,
  ...
)
}
\arguments{
\item{object}{MachineShop \code{\link[MachineShop]{ModelSpecification}} object}

\item{times}{how many total grid points will be evaluated (including initial points evaluated)?}

\item{initial_points}{what is the number of randomly sampled points to evaluate initially?}

\item{degree}{= 3 what is the maximum polynomial degree for numeric and integer parameters? First order interactions will be included regardless. Goes from 3 to and including 9.}

\item{plot_predicted_performance}{should the predicted performance of unevaluated grid points be plotted?}

\item{max_grid_rows}{what is the maximum design-matrix number of rows?}

\item{max_knots}{what is the maximum number of knots to include (for fitting splines)?}

\item{obs_per_spline_coefficient}{function that takes as a single argument the number of knots, and returns the number of observations required to add that knot}

\item{random}{when algorithm fails, should all grid points or a single randomly sampled one (default) be evaluated?}

\item{use_transf}{should a default variance-stabilizing transformation be used?}

\item{log_outcome}{should the metric being modeled be log-transformed as an outcome?}

\item{logit_outcome}{should the metric being modeled be logit-transformed as an outcome?}

\item{save_model_fit}{should the design matrix and beta coefficients be saved to an RData file?}

\item{include_numeric_interactions}{should numeric-numeric interactions be included?}

\item{include_high_order_interactions}{should degree 3-and-higher numeric polynomial terms and categorical interactions be included?}

\item{include_quadratic_interactions}{should degree 2 numeric polynomial terms and categorical interactions be included?}

\item{qualify_numeric_numeric_interactions}{should three-way interactions between categorical and first-order numeric pairings be included?}

\item{save_model_fit_lable}{label for the data file saving the beta coefficients and design matrix,}

\item{generic_tolerance}{tolerance for convergence on a variety of tasks}

\item{cardinal_knots}{should the knots be spaced evenly?}

\item{latin_grid}{should numerics (not integers) be latin-grid'd?}

\item{enforce_positive_definiteness}{should we ensure that each polynomial partition unconstrained is full rank?}

\item{acquisition_function}{custom acquisition function}

\item{print_colnames}{should the column names of the design matrix be printed?}

\item{stopping_rule_total}{stopping rule will be considered after this many grid points have been evaluated since improvement}
}
\value{
A model specification object ready to be fit, equipped with Thompson Sampling with Smoothing Splines.
}
\description{
Equip ModelSpecification objects in MachineShop with a Thompson sampler with smoothing splines for tuning
}
\details{
A custom set_optim_method to be paired with MachineShop for tuning models and inputs using Thompson Sampling with Smoothing Splines.
}
\examples{
\donttest{

library(recipes)
library(MachineShop)
## tune LASSO for l1-norm penalty term, with 26 grid point evaluations
tuned_model <- 
  recipe(sale_amount ~ .,data = ICHomes) \%>\%    
  step_normalize(all_numeric_predictors()) \%>\%  
  ModelSpecification(
    model = TunedModel(        
      GLMNetModel,                           
      grid = expand_params(                  
        lambda = seq(0,2500,len = 1000)))) \%>\%                     
  set_optim_thompson(times = 26)  \%>\%           
  fit 
}

}
\seealso{
\code{\link[MachineShop]{set_optim_method}}, \code{\link[MachineShop]{ModelSpecification}}
}
