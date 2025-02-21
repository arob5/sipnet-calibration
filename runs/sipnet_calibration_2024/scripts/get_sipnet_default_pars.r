#
# get_sipnet_default_pars.r
#
# The main point of this script is to save a complete list of SIPNET parameters, 
# along with their defaults and default bounds. While this is already 
# provided in pecan/models/sipnet/inst/template.param, the parameters in 
# this template follow the SIPNET naming standards. This file constructs 
# a parameter list based on the valid PEcAn parameters. The associated SIPNET
# parameter name is included here as well to provide a mapping between the 
# parameter sets; this mapping between the two naming conventions is taken from
# the code in pecan/models/sipnet/write.configs.SIPNET.R. Unlike the 
# template.param file, this script does not include the initial conditions in 
# the parameters table. In some cases there is not a one-to-one mapping between
# the parameters; these cases are noted in the "note" column of the saved table.
#
# One issue with the default priors is that they are based on the default bounds
# in template.param. Sometimes these bounds are legitimate; e.g., a parameter
# cannot be negative. However, many of these bounds are largely artificial
# and an unbounded prior would probably be a better fit.
#
# Andrew Roberts
# 

# Questions:
# fineRootFrac/coarseRootFrac: are these associated with the ICs, not params?

library(data.table)

base_dir <- file.path("/projectnb", "dietzelab", "arober", "sipnet_calibration", 
                      "runs", "sipnet_calibration_2024")
out_dir <- file.path(base_dir, "param")


par_list <- list()

par_list[[1]] <- list(pecan_name="Amax", sipnet_name="aMax", default=112.0, 
                      lower=0, upper=20, note="this param gets modified by SLA")
par_list[[2]] <- list(pecan_name="AmaxFrac", sipnet_name="aMaxFrac", 
                      default=0.76, lower=0.66, upper=0.86)
par_list[[3]] <- list(pecan_name="extinction_coefficient", 
                      sipnet_name="attenuation", 
                      default=0.58, lower=0.3800, upper=0.620000)
par_list[[4]] <- list(pecan_name="leaf_respiration_rate_m2", 
                      sipnet_name="", 
                      default=NA, lower=NA, upper=NA,
                      note="used to compute baseFolRespFrac := max(min(leaf_respiration_rate_m2/Amax, 1), 0)")
par_list[[5]] <- list(pecan_name="Vm_low_temp", 
                      sipnet_name="psnTMin", 
                      default=2, lower=-8.0000, upper=8.000000)
par_list[[6]] <- list(pecan_name="growth_resp_factor", 
                      sipnet_name="growthRespFrac", 
                      default=0.2, lower=0.1, upper=0.3)
par_list[[7]] <- list(pecan_name="half_saturation_PAR", 
                      sipnet_name="halfSatPar", 
                      default=17, lower=4, upper=27)
par_list[[8]] <- list(pecan_name="stomatal_slope.BB", 
                      sipnet_name="m_ballBerry", 
                      default=3.890000, lower=0.500000, upper=15.000000)
par_list[[9]] <- list(pecan_name="dVPDSlope", 
                      sipnet_name="dVpdSlope", 
                      default=0.05, lower=0.010, upper=0.250000)
par_list[[10]] <- list(pecan_name="dVpdExp", 
                       sipnet_name="dVpdExp", 
                       default=2.0, lower=1.0000, upper=3.000000)
par_list[[11]] <- list(pecan_name="leaf_turnover_rate", 
                       sipnet_name="leafTurnoverRate", 
                       default=0.13, lower=0.001, upper=1.0)
par_list[[12]] <- list(pecan_name="wueConst", 
                       sipnet_name="wueConst", 
                       default=10.9, lower=0.01, upper=109)
par_list[[13]] <- list(pecan_name="veg_respiration_Q10", 
                       sipnet_name="vegRespQ10", 
                       default=2.0, lower=1.40000, upper=2.60000)
par_list[[14]] <- list(pecan_name="stem_respiration_rate", 
                       sipnet_name="", 
                       default=NA, lower=NA, upper=NA,
                       note="used as intermediate quantity")
par_list[[15]] <- list(pecan_name="root_turnover_rate", 
                       sipnet_name="fineRootTurnoverRate", 
                       default=0.137, lower=0.00100000, upper=1.000000)
par_list[[16]] <- list(pecan_name="fine_root_respiration_Q10", 
                       sipnet_name="fineRootQ10", 
                       default=2.600000, lower=1.400000, upper=5.000000)
par_list[[17]] <- list(pecan_name="root_respiration_rate", 
                       sipnet_name="", 
                       default=NA, lower=NA, upper=NA,
                       note="used as intermediate quantity")
par_list[[18]] <- list(pecan_name="coarse_root_respiration_Q10", 
                       sipnet_name="coarseRootQ10", 
                       default=2.600000, lower=1.400000, upper=5.000000)
par_list[[19]] <- list(pecan_name="root_allocation_fraction", 
                       sipnet_name="fineRootAllocation", 
                       default=0.400000, lower=0.000000, upper=0.600000)
par_list[[20]] <- list(pecan_name="wood_allocation_fraction", 
                       sipnet_name="woodAllocation", 
                       default=0.200000, lower=0.000000, upper=0.600000)
par_list[[21]] <- list(pecan_name="leaf_allocation_fraction", 
                       sipnet_name="leafAllocation", 
                       default=.20, lower=0, upper=1)
par_list[[22]] <- list(pecan_name="wood_turnover_rate", 
                       sipnet_name="woodTurnoverRate", 
                       default=0.014000, lower=0, upper=1.0)
par_list[[22]] <- list(pecan_name="soil_respiration_Q10", 
                       sipnet_name="soilRespQ10", 
                       default=2.9, lower=1.400000, upper=5)
par_list[[23]] <- list(pecan_name="som_respiration_rate", 
                       sipnet_name="baseSoilResp", 
                       default=0.06, lower=0.003, upper=0.6)
par_list[[24]] <- list(pecan_name="turn_over_time", 
                       sipnet_name="litterBreakdownRate", 
                       default=0.4, lower=0.13, upper=1.2)
par_list[[25]] <- list(pecan_name="frozenSoilEff", 
                       sipnet_name="frozenSoilEff", 
                       default=1.000000, lower=0, upper=1)
par_list[[26]] <- list(pecan_name="frozenSoilFolREff", 
                       sipnet_name="frozenSoilFolREff", 
                       default=0, lower=0, upper=1)
par_list[[27]] <- list(pecan_name="soilWHC", 
                       sipnet_name="soilWHC", 
                       default=12.0, lower=0.1, upper=36.000000)
par_list[[28]] <- list(pecan_name="immedEvapFrac", 
                       sipnet_name="immedEvapFrac", 
                       default=0.100000, lower=0.0, upper=0.2)
par_list[[29]] <- list(pecan_name="leafWHC", 
                       sipnet_name="leafPoolDepth", 
                       default=0.100000, lower=0.000000, upper=0.200000)
par_list[[30]] <- list(pecan_name="waterRemoveFrac", 
                       sipnet_name="waterRemoveFrac", 
                       default=0.088000, lower=0.001000, upper=0.160000)
par_list[[31]] <- list(pecan_name="fastFlowFrac", 
                       sipnet_name="fastFlowFrac", 
                       default=0, lower=0, upper=0.2)
par_list[[32]] <- list(pecan_name="rdConst", 
                       sipnet_name="rdConst", 
                       default=300.000000, lower=1.0, upper=1500)
par_list[[33]] <- list(pecan_name="GDD", 
                       sipnet_name="gddLeafOn", 
                       default=500, lower=100, upper=900)
par_list[[34]] <- list(pecan_name="fracLeafFall", 
                       sipnet_name="fracLeafFall", 
                       default=1, lower=0, upper=1)
par_list[[35]] <- list(pecan_name="leafGrowth", 
                       sipnet_name="leafGrowth", 
                       default=126, lower=0, upper=252)


for(i in seq_along(par_list)) {
  if(is.null(par_list[[i]]$note)) par_list[[i]]$note <- NA_character_
}
dt <- rbindlist(par_list, use.names=TRUE)

# Save parameter table.
out_path_table <- file.path(out_dir, "sipnet_pecan_params.csv")
fwrite(dt, out_path_table)

# Construct default priors based on default parameter bounds.
out_path_priors <- file.path(out_dir, "priors", "sipnet_params_default_priors.csv")
default_priors <- dt[!is.na(lower) & !is.na(upper)]
default_priors <- default_priors[, .(par_name=pecan_name, distn="unif",
                                     parama=lower, paramb=upper)]

fwrite(default_priors, out_path_priors)

# Split default prior file in two parts based on plant vs soil params.
soil_param <- grepl("soil", default_priors$par_name, ignore.case=TRUE)
default_priors_soil <- default_priors[soil_param,]
default_priors_plant <- default_priors[!soil_param,]

out_path_prior_soil <- file.path(out_dir, "priors", "sipnet_params_default_priors_soil.csv")
out_path_prior_plant <- file.path(out_dir, "priors", "sipnet_params_default_priors_plant.csv")

fwrite(default_priors_soil, out_path_prior_soil)
fwrite(default_priors_plant, out_path_prior_plant)

