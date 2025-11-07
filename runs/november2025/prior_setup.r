# prior_setup.r
#
# This file defines the set up parameters (traits and initial conditions) to
# calibrate, as well as prior distributions on these parameters. 
#

prior_list <- list()

prior_list$root_turnover_rate <- list(
  param_name = "root_turnover_rate",
  dist_name = "Beta",
  constraint = c(0,1),
  length = 1L,
  dist_params = list(shape1=2, shape2=4)
)

prior_list$leafGrowth <- list(
  param_name = "leafGrowth",
  dist_name = "Uniform",
  constraint = c(0, 252),
  length = 1L,
  dist_params = list(min=0, max=252)
)

# This puts most mass between [0,2]. Need to check that this is reasonable.
prior_list$stem_respiration_rate <- list(
  param_name = "stem_respiration_rate",
  dist_name = "Gamma",
  constraint = c(0, Inf),
  length=1L,
  dist_params = list(shape=2, rate=3)
)

# Puts most mass between [5, 35] (degrees celsius)
prior_list$psnTOpt <- list(
  param_name = "psnTOpt",
  dist_name = "Gamma",
  constraint = c(0, Inf),
  length=1L,
  dist_params = list(shape=20, rate=1)
)

prior_list$half_saturation_PAR <- list(
  param_name = "half_saturation_PAR",
  dist_name = "Uniform",
  constraint = c(4, 27),
  length=1L,
  dist_params = list(min=4, max=27)
)

# Puts most mass in [15, 40].
# See: https://onlinelibrary.wiley.com/doi/10.1111/pce.14815
prior_list$Amax <- list(
  param_name = "Amax",
  dist_name = "Gamma",
  constraint = c(0, Inf),
  length=1L,
  dist_params = list(shape=40, rate=1.5)
)


