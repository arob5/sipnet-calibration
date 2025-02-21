#
# exploring_pecan_settings_obj.r
# A script for learning the basic structure of the pecan settings object, 
# which is loaded from the XML file. 
#
# Andrew Roberts
#
# Working directory: sipnet_calibration

library(PEcAn.settings) # Required for loading pecan settings object. 

base_dir <- getwd() # Must be "sipnet_calibration"
pecan_settings_path_dongchen <- file.path(base_dir, "pecan_dongchen.xml")
pecan_settings_path_istem <- "/projectnb/dietzelab/istfer/outputs/istfer/pecan.pda1000021382.xml"


#
# Dongchen's settings 
#
# Copied from: 
# /projectnb/dietzelab/dongchen/All_NEON_SDA/NEON42/SDA.xml
#

# Open and read in settings file for PEcAn run.
settings <- PEcAn.settings::read.settings(inputfile=pecan_settings_path_dongchen)

# "MultiSettings" object, inherits from "list". 
class(settings)

# Looking at some setings. 
names(settings)
settings$outdir
settings$rundir
settings$modeloutdir
names(settings$model) # Defines which model (e.g. SIPNET) to use. 
settings$model # Not sure what all of these settings mean. 
names(settings$pfts) # Appears to be a list of PFTs to use; in this case, 3 PFTs. 

# PFT settings. 
names(settings$pfts[[1]]) # Assuming "posteriorid" points to parameter calibration posterior samples. 
settings$pfts[[1]]$name

# "Run" settings. 
names(settings$run) # One element per run. 
names(settings$run$settings.1) # A run is defined by site, time range, and inputs. 
names(settings$run$settings.1$site) # A site is defined by its ID, meterorology time range, and geographic coords. 
settings$run$settings.1$start.date # How does this date differ from/interact with the site dates?
settings$run$settings.1$end.date
names(settings$run$settings.1$inputs) # Inputs: meteorology, PFT site, pool initial conditions. 
names(settings$run$settings.1$inputs$met) # Meterology: source, output, id, path. 
length(settings$run$settings.1$inputs$met$path) # A bunch of paths. 
settings$run$settings.1$inputs$pft.site$path # How does site PFT differ from other PFT settings? 
names(settings$run$settings.1$inputs$poolinitcond) # Initial conditions: source, output, ensemble, path
settings$run$settings.1$inputs$poolinitcond$ensemble # Is this the number of initial condition ensemble members? 
length(settings$run$settings.1$inputs$poolinitcond$path) # Seems like one path per ensemble member. 


#
# Istem's settings 
#

# Open and read in settings file for PEcAn run.
settings_istem <- PEcAn.settings::read.settings(inputfile=pecan_settings_path_istem)

# Has "assim.batch" and "meta.analysis" elements. 
names(settings_istem)

# SIPNET model 
settings_istem$model

# Meta-analysis
settings_istem$meta.analysis # iter, random.effects, threshold, update 

# assim.batch: 
# Looks like essentially all of the settings are stored in here: Priors, 
# MCMC info, likelihood, parameter info, emulator info. 
settings_istem$assim.batch







