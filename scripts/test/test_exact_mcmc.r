# 
# test_exact_mcmc.r
# Testing exact MCMC sampling. Using linear Gaussian example so that we can 
# compare to the exact closed-form posterior. These are tests of the PEcAn
# MCMC code, but no PEcAn functionality is actually used here. Separate tests
# are required in the case that the log-likelihood depends on a PEcAn model 
# run.
#
# Andrew Roberts
#

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(data.table)
library(assertthat)

# Directories. 
base_dir <- file.path("/projectnb", "dietzelab", "arober", "sipnet_calibration")
src_dir <- file.path("/projectnb", "dietzelab", "arober", "gp-calibration", "src")
pecan_base_dir <- file.path(base_dir, "runs", "sipnet_calibration_2024")
calibration_dir <- file.path(pecan_base_dir, "calibration")
pecan_src_dir <- file.path(base_dir, "src")

# Source gp-calibration files.
source(file.path(src_dir, "general_helper_functions.r"))
source(file.path(src_dir, "inv_prob_test_functions.r"))
source(file.path(src_dir, "plotting_helper_functions.r"))
source(file.path(src_dir, "seq_design.r"))
source(file.path(src_dir, "mcmc_helper_functions.r"))
source(file.path(src_dir, "ens_Kalman_inversion.r")) # For calc_lin_Gauss_cond_moments_obs()
source(file.path(src_dir, "gp_helper_functions.r")) # For get_bounds()

# Source PEcAn code that is being tested.
source(file.path(pecan_src_dir, "param_calibration_functions.r"))
source(file.path(pecan_src_dir, "mcmc_pecan.r"))

# Random number generator seed.
seed <- 96878
set.seed(seed)

# ------------------------------------------------------------------------------
# Setup statistical model
# ------------------------------------------------------------------------------

# Parameter dimension (d) and observation dimension (d).
d <- 2L
p <- 2L

# Mean and covariance of Gaussian prior on parameter.
m0 <- rep(0,d)
L0 <- matrix(rnorm(d*d), nrow=d, ncol=d)
C0 <- tcrossprod(L0)

# Covariance of noise term.
L_Sigma <- matrix(rnorm(p*p), nrow=p, ncol=p)
Sigma <- tcrossprod(L_Sigma)

# Linear forward model. 
G <- matrix(rnorm(p*d), nrow=p, ncol=d)

# Define the observed data vector. 
y <- c(1,2)

# Log likelihood function. Note that we're treating the noise covariance 
# as fixed/known here, so we don't include the log-det normalizing constant
# in the Gaussian likelihood. The `llik` argument to the MCMC functions 
# must be a function that includes `...` as an argument. This function is 
# also vectorized to accept a matrix of different parameter values.
llik <- function(par, ...) {
  if(is.null(dim(par))) par <- matrix(par, nrow=1L)
  -0.5 * colSums(forwardsolve(L_Sigma, add_vec_to_mat_cols(-y, G %*% t(par)))^2)
}


# Log prior.
lprior <- function(par) {
  if(is.null(dim(par))) par <- matrix(par, nrow=1L)
  -0.5 * colSums(forwardsolve(L0, add_vec_to_mat_cols(-m0, t(par)))^2)
}


# ------------------------------------------------------------------------------
# Exact posterior
# ------------------------------------------------------------------------------

# Get mean and covariance of the true posterior density.
true_post_moments <- calc_lin_Gauss_cond_moments_obs(G, y, m0, C0=C0, Sig=Sigma)
m_post <- true_post_moments$mean
L_post <- t(chol(true_post_moments$cov))
cov_post <- tcrossprod(L_post)

# Exact posterior density.
lpost_true <- function(par) {
  if(is.null(dim(par))) par <- matrix(par, nrow=1L)
  -0.5 * colSums(forwardsolve(L_post, add_vec_to_mat_cols(-m_post, t(par)))^2)
}

# ------------------------------------------------------------------------------
# Metropolis-Hastings Sampling
# ------------------------------------------------------------------------------

# Sample initial condition from prior.
par_init <- m0 + drop(L0 %*% matrix(rnorm(d), ncol=1L))

# Parameter names required by `mcmc_mh()`.
par_names <- paste0("par", 1:d)

# Number MCMC iterations.
n_itr <- 100000L

# Run Metropolis-Hastings, using default proposal covariance adaptation
# settings.
mcmc_output <- mcmc_mh(llik, lprior, par_names, n_itr, par_init)

# Drop first half as burn-in.
samp_dt <- format_mcmc_output(mcmc_output$samp, test_label="mcmc_mh")
samp_dt <- select_mcmc_itr(samp_dt, itr_start=round(n_itr/2))

# ------------------------------------------------------------------------------
# MCMC Diagnostics
# ------------------------------------------------------------------------------

# Trace plots.
trace_plots <- get_trace_plots(samp_dt)
for(plt in trace_plots) plot(plt)

# R-hat.
r_hat_results <- calc_R_hat(samp_dt, within_chain=TRUE)
print(r_hat_results$R_hat_vals)


# ------------------------------------------------------------------------------
# Comparing to known posterior.
#   Comparing contours of (unnormalized) posterior density.
# ------------------------------------------------------------------------------

# Compare means and covs.
samp_mat <- select_mcmc_samp_mat(samp_dt, "mcmc_mh", "par")[,par_names]
m_post_mcmc <- colMeans(samp_mat)
cov_post_mcmc <- cov(samp_mat)

print("Mean error:")
print(m_post - m_post_mcmc)

print("Cov error:")
print(cov_post - cov_post_mcmc)

# Compare 1d marginals.
J <- 100000L
samp_exact <- add_vec_to_mat_rows(m_post, t(L_post %*% matrix(rnorm(J*d), nrow=d, ncol=J)))
colnames(samp_exact) <- par_names
samp_dt <- append_samples_mat(samp_dt, samp_exact, "par", "exact")

hist_plots <- get_hist_plot_comparisons(samp_dt, test_label_baseline="exact")
for(plt in hist_plots) plot(plt)

# Compare 2d KDE.
df_mcmc_mh <- data.frame(select_mcmc_samp_mat(samp_dt, "mcmc_mh", "par"))[,par_names]

kde2d_comparison <- ggplot() + 
                    geom_density_2d(aes(x=par1, y=par2), data.frame(samp_exact)) +
                    geom_density_2d(aes(x=par1, y=par2), df_mcmc_mh, color="red")

plot(kde2d_comparison)








# Draw samples from prior. These will be used to define the plotting grid.
J <- 500L
U <- add_vec_to_mat_rows(m0, t(L0 %*% matrix(rnorm(J*d), nrow=d, ncol=J)))
colnames(U) <- par_names

# Define grid for plotting.
n_grid <- 50^2
u_bounds <- get_bounds(U) # This function is defined in `gp_helper_functions.r`
u_grid <- get_batch_design("tensor_product_grid", bounds=u_bounds, N_batch=n_grid)
colnames(u_grid) <- par_names

# Plot log density contours.
lpost_true_grid <- lpost_true(u_grid)
df_true_lpost <- data.frame(u_grid, lpost=lpost_true_grid)

plt_true_lpost <- ggplot( ) +
                  geom_contour_filled(aes(x=par1, y=par2, z=lpost), df_true_lpost,
                                      bins=20) + 
                  geom_point(aes(x=m_post[1], y=m_post[2]), color="red", shape=8) +
                  geom_point(aes(x=m0[1], y=m0[2]), color="blue", shape=8) +
                  theme_minimal() + theme(legend.position="none")
plot(plt_true_lpost)

# Overlay KDE from MCMC samples.
df_samp <- data.frame(select_mcmc_samp_mat(samp_dt, "mcmc_mh", "par"))


plt_true_lpost <- plt_true_lpost + geom_density_2d(aes(x=par1, y=par2), df_samp)
plot(plt_true_lpost)

plt_true_lpost +                   
  geom_point(aes(x=u1, y=u2), data.frame(U), color="gray", alpha=0.8) + 
  ggtitle("True Log Posterior Density vs. Prior Samples")












