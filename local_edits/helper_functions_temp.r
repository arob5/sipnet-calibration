#
# helper_functions_temp.r
# Helper functions for getting parameter calibration up an running for SIPNET. 
# These are not intended to be production code, just short-term solutions. 
#
# Andrew Roberts
# 
# Dependencies: dietzelab/arober/sipnet_calibration/local_edits/local_edits.r
# 

gen_prior_param_design <- function(settings) {
  # Constructs a set of design inputs for the calibration parameters relative 
  # to specified prior distributions. This roughly follows the current PDA 
  # code, but leverages simplified/modified functions in the `local_edits.r`
  # file. 
  #
  # Args:
  #    - settings: list, the PEcAn settings list. 
  #
  # Returns: 
  #   Returns a list with the following named elements: 
  #    - design_list_by_pft: list of length equal to the number of PFTs, with names 
  #      set to the PFT names. Each element is a matrix with dimension 
  #      (number design points, number params), where I believe number params is 
  #      the number of params (falling within that PFT with specified prior 
  #      distributions not the number of calibration parameters, and I also think not 
  #      necessarily the total number of parameters of the model). Each row of 
  #      the matrix is one parameter set; the subset of columns associated with 
  #      the calibration parameters represent the input design, while the 
  #      non-calibration parameters are set to constant values given by their 
  #      prior medians. 
  #    - n_design: Integer, the number of design points. 
  #    - n_param_by_pft: list, of length equal to the number of PFTs, with names 
  #                      set to the PFT names. The number of parameters with priors 
  #                      for each PFT (agrees with the dimensions of the matrices
  #                      in `design_list_by_pft`, see above). 
  #    - dim_param: The dimension of the parameter calibration space; i.e., 
  #                 the total number of calibration parameters. 
  #    - cal_params_by_pft: list, of length equal to the number of PFTs, with names 
  #                         set to the PFT names. Each element is a character vector 
  #                         of the calibration parameter names within each 
  #                         respective PFT. 

  # Generate list over PFTs, with each element being a data.frame specifying the 
  # prior distributions for the parameters within each PFT. 
  prior_dist_list <- pda.load.priors_test_edit(settings)
  
  # Convert the priors into their relevant distribution, density, sampling, etc. functions. 
  # This results in a list over PFTs, each element itself being a list with elements 
  # "dprior", "rprior", "qprior", etc. Each of these in turn is also a list of 
  # length equal to the number of parameters for that PFT, with each element 
  # being an R expression (not an R function). 
  prior_funcs_list <- lapply(prior_dist_list, PEcAn.assim.batch::pda.define.prior.fn)
  
  # Store indices associated with parameters marked for calibration. This results in 
  # a list over PFTs, each element being a vector of indices that can be used 
  # to subset `param_names`. I feel like it would make more sense to have all
  # parameters in a single data.frame with a column "pft" and a logical column 
  # "calibrate". 
  param_names <- lapply(prior_dist_list, rownames)
  pft_names <- sapply(settings$pfts, `[[`, "name")
  cal_param_idx <- vector("list", length(settings$pfts))
  names(cal_param_idx) <- pft_names
  
  for(i in seq_along(settings$pfts)) {
    pft_name <- settings$pfts[[i]]$name
    if(pft_name %in% names(settings$assim.batch$param.names)) {
      cal_param_idx[[i]] <- which(param_names[[i]] %in% settings$assim.batch$param.names[[pft_name]])
    }
  }
  
  cal_params_by_pft <- lapply(pft_names, function(pft) param_names[[pft]][cal_param_idx[[pft]]])
  names(cal_params_by_pft) <- pft_names

  # Dimension of the parameter space for calibration. 
  dim_param <- length(unlist(cal_param_idx))
  
  # Generates design points: calibration parameter design generated with 
  # latin hypercube sampling; non-calibration parameters fixed at their prior 
  # medians. Note that the design is generated one PFT at a time, which isn't a problem 
  # for latin hypercube sampling (which can be constructed variable-by-variable), 
  # but this would be a problem for other designs (e.g. minimax/maximin).
  n_param_by_pft <- sapply(prior_dist_list, nrow)
  n_design <- as.integer(settings$ensemble$size)
  param_design_list <- lapply(seq_along(settings$pfts),
                              function(i) PEcAn.assim.batch::pda.generate.knots(n.knot=n_design, 
                                                                                sf=NULL, probs.sf=NULL,
                                                                                n.param.all=n_param_by_pft[i],
                                                                                prior.ind=cal_param_idx[[i]],
                                                                                prior.fn=prior_funcs_list[[i]],
                                                                                pname=param_names[[i]])[["params"]])
  names(param_design_list) <- pft_names
  
  return(list(design_list_by_pft=param_design_list, n_design=n_design,
              n_param_by_pft=n_param_by_pft, dim_param=dim_param, 
              cal_params_by_pft=cal_params_by_pft))
  
}

