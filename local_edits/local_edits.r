#
# local_edits.r
# Copies of PEcAn functions, modified for various testing purposes. These are not 
# official PEcAn updates that should be committed. They are simply local 
# hacks for convenience in running various tests. The function naming 
# convention is "<original_function_name>_test_edit". Comments above 
# the edited functions should (1) include the file that the function was 
# copied from, (2) the Git hash, and (3) a summary of the function edits. 
# Note that the Git hash can be found by running git rev-parse HEAD; see reference: 
# https://stackoverflow.com/questions/949314/how-do-i-get-the-hash-for-the-current-commit-in-git
#
# Andrew Roberts
# 

# Copied from: pecan/modules/assim.batch/pda.utils.R
# Git hash: 05edb49
# Edits: Only change is to add `IC` argument to `pda.init.run()` so that it can 
#        be passed to `my.write.config()` to set the initial conditions.
pda.init.run_test_edit <- function(settings, con, my.write.config, workflow.id, params,
                                   n = ifelse(is.null(dim(params)), 1, nrow(params)),
                                   run.names = paste("run", 1:n, sep = "."), IC=NULL) {
  run.ids <- rep(NA, n)
  for (i in seq_len(n)) {
    ## set RUN.ID
    if (!is.null(con)) {
      paramlist <- run.names[i]
      run.ids[i] <- PEcAn.DB::db.query(
        paste0(
          "INSERT INTO runs",
          "(model_id, site_id, start_time, finish_time, outdir,",
          "ensemble_id, parameter_list) ",
          "values ('",
          settings$model$id, "','", settings$run$site$id, "','", settings$run$start.date, "','",
          settings$run$end.date, "','", settings$modeloutdir , "','",
          settings$assim.batch$ensemble.id, "','", paramlist,
          "') RETURNING id"),
        con)
      
    } else {
      run.ids[i] <- run.names[i]
    }
    dir.create(file.path(settings$rundir, run.ids[i]), recursive = TRUE)
    dir.create(file.path(settings$modeloutdir, run.ids[i]), recursive = TRUE)
    
    # # works for one PFT
    # trait.values = list(pft=as.list(params[i,]),env=NA)
    # newnames <- sapply(settings$pfts, "[[", "name")
    # names(trait.values)[which(!(names(trait.values) %in% 'env'))] <- newnames
    
    ## write config
    do.call(my.write.config, args = list(defaults = settings$pfts,
                                         trait.values = lapply(params,
                                                               function(x, n) { as.data.frame(x[n, , drop=FALSE]) },
                                                               n = i),
                                         settings = settings,
                                         run.id = run.ids[i], IC=IC)) # EDIT: added IC=IC
    
    # Identifiers for ensemble 'runtype'
    if (settings$assim.batch$method == "bruteforce" |
        settings$assim.batch$method == "bruteforce.bs" |
        settings$assim.batch$method == "bayesian.tools") {
      ensemble.type <- "pda.MCMC"
    } else if (settings$assim.batch$method == "emulator") {
      ensemble.type <- "pda.emulator"
    }
    
    ## write a README for the run
    cat("runtype     : ",
        paste("pda", settings$assim.batch$method, sep = "."), "\n", "workflow id : ",
        as.character(workflow.id), "\n",
        "ensemble id : ", as.character(settings$assim.batch$ensemble.id), "\n",
        "chain       : ", settings$assim.batch$chain, "\n",
        "run         : ", run.names[i], "\n",
        "run id      : ", as.character(run.ids[i]), "\n",
        "pft names   : ", as.character(lapply(settings$pfts, function(x) x[['name']])), "\n",
        "site        : ", settings$run$site$name, "\n",
        "site  id    : ", settings$run$site$id, "\n",
        "met data    : ", settings$run$site$met, "\n",
        "start date  : ", settings$run$start.date, "\n",
        "end date    : ", settings$run$end.date, "\n",
        "hostname    : ", settings$host$name, "\n",
        "rundir      : ", file.path(settings$host$rundir, run.ids[i]), "\n",
        "outdir      : ", file.path(settings$host$outdir, run.ids[i]), "\n",
        file = file.path(settings$rundir, run.ids[i], "README.txt"))
    
    ## add the job to the list of runs
    append <- ifelse(i == 1, FALSE, TRUE)
    cat(as.character(run.ids[i]),
        file = file.path(settings$rundir, "runs.txt"),
        sep = "\n",
        append = append)
  }  # end for
  
  return(unlist(run.ids))
} # pda.init.run


# Copied from: pecan/modules/assim.batch/pda.utils.R
# Git hash: ede6ff5
# Edits: Trimming down to only the essentials for my test run. Removing anything to 
#        to with the database or a "run extension". Setting this up to run locally, 
#        so that prior files must be saved to the paths specified in `pfts$pft$outdir`. 
pda.load.priors_test_edit <- function(settings) {
  
  # Loop over PFTs and load the prior distribution Rdata file for each PFT. These 
  # files contain priors for all parameters falling within each PFT. For SIPNET, 
  # there are typically 2 PFTs: a plant PFT and a soil PFT. This results in 
  # `prior_dists`, a list of data.frames (one element per PFT). 
  pft_prior_paths <- sapply(settings$pfts, function(pft) pft[["outdir"]])
  pft_prior_paths <- file.path(pft_prior_paths, "prior.distns.Rdata")
  prior_dists_env <- new.env()
  prior_dists <- list()
  for(i in seq_along(pft_prior_paths)) {
    pft_name <- settings$pfts[[i]]$name
    load(pft_prior_paths[[i]], envir=prior_dists_env)
    prior_dists[[pft_name]] <- prior_dists_env$prior.distns
  }
  
  # Ensure all parameters marked as calibration parameters have been provided 
  # priors. 
  calibration_params <- unlist(settings$assim.batch$param.names)
  params_with_priors <- unlist(lapply(prior_dists, rownames))
  calibration_params_no_priors <- setdiff(calibration_params, params_with_priors)
  if(length(calibration_params_no_priors) > 0) {
    stop("Calibration parameters missing priors: ", 
         paste(calibration_params_no_priors, collapse=", "))
  }
  
  return(prior_dists)
  
}


# Copied from: pecan/modules/assim.batch/pda.utils.R
# Git hash: ede6ff5
# Edits: TODO 
pda.generate.knots_test_edit <- function(n.knot, sf, probs.sf, n.param.all, prior.ind, prior.fn, pname) {
  # By default, all parameters will be fixed at their median
  probs <- matrix(0.5, nrow = n.knot, ncol = n.param.all)
  
  # if the params are going to be scaled leave them out, if sf is NULL nothing happens
  inds <- prior.ind[!prior.ind %in% which(pname %in% sf)]
  
  if(length(inds) !=0){
    # Fill in parameters to be sampled with probabilities sampled in a LHC design
    probs[, inds] <- PEcAn.emulator::lhc(t(matrix(0:1, ncol = length(inds), nrow = 2)), n.knot)
  }
  
  inds <- prior.ind[prior.ind %in% which(pname %in% sf)]
  
  if(!is.null(sf) & length(inds) !=0){
    match.ind <- sapply(pname[inds], function(x) which(sf == x))
    probs[, inds] <- probs.sf[, match.ind]
  }
  
  # Convert probabilities to parameter values
  params <- NA * probs
  for (i in seq_len(n.param.all)) {
    params[, i] <- eval(prior.fn$qprior[[i]], list(p = probs[, i]))
  }
  colnames(params) <- pname
  colnames(probs) <- pname
  
  return(list(params = params, probs = probs))
}


# Copied from: pecan/base/workflow/R/run.write.configs.R
# Git hash: 6425b3c2733a1a6ee9e373085a308f9e2b0ca771
# Edits: Significant changes. Maintaining the bare minimum to be able to execute 
# a forward model evaluation by reading in initial conditions, drivers, and parameters
# from file. In particular, all code involving the bety database has been removed. 
# `params_path` is the path to an .RDdata file containing the parameter values. 
# This function requires the run to be an "ensemble run". 
run.write.configs_test_edit <- function(settings, ens.sample.method="uniform", overwrite=FALSE, 
                                        params_path=NULL) {

  # This function is only applicable for ensemble runs. 
  assert_that("ensemble" %in% names(settings))
  
  # According to documentation, these defaults should be defined somewhere, but 
  # they're not getting filled when I run the code so adding them here. 
  if(is.null(settings$host$rundir)) settings$host$rundir <- file.path(settings$rundir)
  if(is.null(settings$host$outdir)) settings$host$outdir <- file.path(settings$modeloutdir)
  
  # If `outdir` has not been specified for a PFT, set a default directory. 
  for(i in seq_along(settings$pfts)) {
    if(is.null(settings$pfts[[i]]$outdir) || is.na(settings$pfts[[i]]$outdir)) {
      settings$pfts[[i]]$outdir <- file.path(settings$outdir, "pfts", settings$pfts[[i]]$name)
    }
  }
  
  # Select parameters to be used in forward model runs. Read parameters from local file. 
  if(is.null(params_path)) params_path <- file.path(settings$outdir, "samples.Rdata")
  assert_that(file.exists(params_path))
  
  # Read parameters. 
  samples <- new.env()
  load(params_path, envir=samples)
  trait.samples <- samples$trait.samples
  ensemble.samples <- samples$ensemble.samples
  sa.samples <- samples$sa.samples
  runs.samples <- samples$runs.samples
  
  # In this edited function, all that is needed is `ensemble.samples`. 
  assert_that(!is.null(ensemble.samples))  
  
  pft.names <- names(trait.samples)
  trait.names <- lapply(trait.samples, names)
  scipen <- getOption("scipen")
  options(scipen=12)
  
  # remove previous runs.txt
  if (overwrite && file.exists(file.path(settings$rundir, "runs.txt"))) {
    PEcAn.logger::logger.warn("Existing runs.txt file will be removed.")
    unlink(file.path(settings$rundir, "runs.txt"))
  }
  
  #
  # Forward model setup. 
  #
  
  model <- settings$model$type
  PEcAn.utils::load.modelpkg(model)
  
  # Check for model-specific write configs.
  my.write.config <- paste0("write.config.", model)
  if (!exists(my.write.config)) {
    PEcAn.logger::logger.error(my.write.config, 
                               "does not exist, please make sure that the model package contains a function called", 
                               my.write.config)
  }
  
  # Prepare for model output.  Clean up any old config files (if exists)
  my.remove.config <- paste0("remove.config.", model)
  if (exists(my.remove.config)) {
    do.call(my.remove.config, args=list(settings$rundir, settings))
  }
  
  #
  # Writing configs for a forward ensemble run. 
  #
  ens.runs <- PEcAn.uncertainty::write.ensemble.configs(defaults=settings$pfts,
                                                        ensemble.samples=ensemble.samples, 
                                                        settings=settings,
                                                        model=model, 
                                                        write.to.db=FALSE)

  # Store output in settings and output variables
  runs.samples <- list()
  runs.samples$ensemble <- ens.run.ids <- ens.runs$runs
  settings$ensemble$ensemble.id <- ens.ensemble.id <- ens.runs$ensemble.id
  ens.samples <- ensemble.samples  # rename just for consistency
    
  # Save ensemble analysis info
  fname <- PEcAn.uncertainty::ensemble.filename(settings, "ensemble.samples", "Rdata", all.var.yr = TRUE)
  save(ens.run.ids, ens.ensemble.id, ens.samples, pft.names, trait.names, file = fname)

    
  PEcAn.logger::logger.info("###### Finished writing model run config files #####")
  PEcAn.logger::logger.info("config files samples in ", file.path(settings$outdir, "run"))
  
  ### Save output from SA/Ensemble runs
  # A lot of this is duplicate with the ensemble/sa specific output above, but kept for backwards compatibility. 
  save(ensemble.samples, trait.samples, sa.samples, runs.samples, pft.names, trait.names, 
       file = file.path(settings$outdir, "samples.Rdata"))
  PEcAn.logger::logger.info("parameter values for runs in ", file.path(settings$outdir, "samples.RData"))
  options(scipen=scipen)
  
  return(invisible(settings))
}


# Copied from: pecan/modules/uncertainty/R/ensemble.R
# Git hash: 6425b3c2733a1a6ee9e373085a308f9e2b0ca771
# Edits: None so far. Just copying for use with browser(). 
write.ensemble.configs_test_edit <- function(defaults, ensemble.samples, settings, model, 
                                             clean = FALSE, write.to.db = TRUE, restart=NULL, rename = FALSE) {

  browser()
  
  my.write.config <- paste("write.config.", model, sep = "")
  my.write_restart <- paste0("write_restart.", model)
  
  # Get the workflow id
  if (!is.null(settings$workflow$id)) {
    workflow.id <- settings$workflow$id
  } else {
    workflow.id <- -1
  }
  
  # Modifying this function to only allow fresh runs - no restarts. 
  assert_that(is.null(restart))
  
  # Setup for sub-sampling parameter ensemble. Note that ensembles of other model  
  # inputs can also be created using this functionality, but my focus here is 
  # on parameter ensembles.
  ensemble.id <- NA
  required_tags <- c("met","parameters")
  samp <- settings$ensemble$samplingspace
  parents <- lapply(samp,'[[', 'parent')
  order <- names(samp)[lapply(parents, function(tr) which(names(samp) %in% tr)) %>% unlist()] 
  samp.ordered <- samp[c(order, names(samp)[!(names(samp) %in% order)])]
  samples <- list()
  
  # Perform the sub-sampling. Note that the "parent" notion here is some way to 
  # link samples across different model components (drivers, parameters, etc.). 
  # I am not using this so the parent will just be NULL. Note that `input.ens.gen` 
  # only performs subsampling for components other than parameters. Since I'm 
  # only dealing with parameter ensembles here it will return NULL. 
  for(i in seq_along(samp.ordered)){
    myparent <- samp.ordered[[i]]$parent
    parent.ids <- NULL
    if(!is.null(myparent)) parent.ids <- samples[[myparent]]
    samples[[names(samp.ordered[i])]] <- input.ens.gen(settings=settings,
                                                       input=names(samp.ordered)[i],
                                                       method=samp.ordered[[i]]$method,
                                                       parent_ids=parent.ids) 
  }

  # if there is a tag required by the model but it is not specified in the xml then I replicate n times the first element 
  required_tags%>%
    purrr::walk(function(r_tag){
      if (is.null(samples[[r_tag]]) & r_tag!="parameters") samples[[r_tag]]$samples <<- rep(settings$run$inputs[[tolower(r_tag)]]$path[1], settings$ensemble$size)
    })
    
  # Reading the site.pft specific tags from xml
  site.pfts.vec <- settings$run$site$site.pft %>% unlist %>% as.character
    
  if (!is.null(site.pfts.vec)) {
    # find the name of pfts defined in the body of pecan.xml
    defined.pfts <- settings$pfts %>% purrr::map('name') %>% unlist %>% as.character
    
    # subset ensemble samples based on the pfts that are specified in the site and they are also sampled from.
    if (length(which(site.pfts.vec %in% defined.pfts)) > 0) {
      ensemble.samples <- ensemble.samples [site.pfts.vec[which(site.pfts.vec %in% defined.pfts)]]
    }
          
    # warn if there is a pft specified in the site but it's not defined in the pecan xml.
    if (length(which(!(site.pfts.vec %in% defined.pfts))) > 0)
      PEcAn.logger::logger.warn(
        paste0(
          "The following pfts are specified for the siteid ",
          settings$run$site$id ,
          " but they are not defined as a pft in pecan.xml:",
          site.pfts.vec[which(!(site.pfts.vec %in% defined.pfts))],
          collapse = ","
        )
      )
  }
    
  # if no ensemble piece was in the xml I replicate n times the first element in params
  if(is.null(samp$parameters)) samples$parameters$samples <- ensemble.samples %>% purrr::map(~.x[rep(1, settings$ensemble$size) , ])
  # This where we handle the parameters - ensemble.samples is already generated in run.write.config and it's sent to this function as arg - 
  if(is.null(samples$parameters$samples)) samples$parameters$samples <- ensemble.samples
  #------------------------End of generating ensembles-----------------------------------
    
  # find all inputs that have an id
  inputs <- names(settings$run$inputs)
  inputs <- inputs[grepl(".id$", inputs)]
    
  # write configuration for each run of the ensemble
  runs <- data.frame()
  for(i in seq_len(settings$ensemble$size)) {
    run.id <- PEcAn.utils::get.run.id("ENS", PEcAn.utils::left.pad.zeros(i, 5), site.id=settings$run$site$id)
    runs[i, "id"] <- run.id
      
    # create folders (cleaning up old ones if needed)
    if(clean) {
      unlink(file.path(settings$rundir, run.id))
      unlink(file.path(settings$modeloutdir, run.id))
    }
    
    dir.create(file.path(settings$rundir, run.id), recursive=TRUE)
    dir.create(file.path(settings$modeloutdir, run.id), recursive=TRUE)
    
    # write run information to disk
    cat("runtype     : ensemble\n",
        "workflow id : ", format(workflow.id, scientific = FALSE), "\n",
        "ensemble id : ", format(ensemble.id, scientific = FALSE), "\n",
        "run         : ", i, "/", settings$ensemble$size, "\n",
        "run id      : ", format(run.id, scientific = FALSE), "\n",
        "pft names   : ", as.character(lapply(settings$pfts, function(x) x[["name"]])), "\n",
        "model       : ", model, "\n",
        "model id    : ", format(settings$model$id, scientific = FALSE), "\n",
        "site        : ", settings$run$site$name, "\n",
        "site  id    : ", format(settings$run$site$id, scientific = FALSE), "\n",
        "met data    : ", samples$met$samples[[i]], "\n",
        "start date  : ", settings$run$start.date, "\n",
        "end date    : ", settings$run$end.date, "\n",
        "hostname    : ", settings$host$name, "\n",
        "rundir      : ", file.path(settings$host$rundir, run.id), "\n",
        "outdir      : ", file.path(settings$host$outdir, run.id), "\n",
        file = file.path(settings$rundir, run.id, "README.txt"))
      
      #changing the structure of input tag to what the models are expecting
      for(input_i in seq_along(settings$run$inputs)) {
        input_tag <- names(settings$run$inputs)[[input_i]]
        if(!is.null(samples[[input_tag]])) {
          settings$run$inputs[[input_tag]][["path"]] <- samples[[input_tag]][["samples"]][[i]]
        }     
      }
      
      my.write.config <- "write.config.SIPNET_test_edit"
      do.call(my.write.config, args=list(defaults=defaults, 
                                         trait.values=lapply(samples$parameters$samples, function(x, n) {x[n, , drop=FALSE]}, n=i), # this is the params
                                         settings=settings, 
                                         run.id=run.id)
      
      )
      cat(format(run.id, scientific=FALSE), file=file.path(settings$rundir, "runs.txt"), sep="\n", append=TRUE)
      
    }
    return(invisible(list(runs=runs, ensemble.id=ensemble.id, samples=samples)))
   
}


write.config.SIPNET_test_edit <- function(defaults, trait.values, settings, run.id, inputs = NULL, IC = NULL,
                                          restart = NULL, spinup = NULL) {
  
  browser()
  
  ### WRITE sipnet.in
  template.in <- system.file("sipnet.in", package = "PEcAn.SIPNET")
  config.text <- readLines(con = template.in, n = -1)
  writeLines(config.text, con = file.path(settings$rundir, run.id, "sipnet.in"))
  
  ### WRITE *.clim
  template.clim <- settings$run$inputs$met$path  ## read from settings
  
  if (!is.null(inputs)) {
    ## override if specified in inputs
    if ("met" %in% names(inputs)) {
      template.clim <- inputs$met$path
    }
  }
  PEcAn.logger::logger.info(paste0("Writing SIPNET configs with input ", template.clim))
  
  # find out where to write run/ouput
  rundir <- file.path(settings$host$rundir, as.character(run.id))
  outdir <- file.path(settings$host$outdir, as.character(run.id))
  if (is.null(settings$host$qsub) && (settings$host$name == "localhost")) {
    rundir <- file.path(settings$rundir, as.character(run.id))
    outdir <- file.path(settings$modeloutdir, as.character(run.id))
  }
  
  # create launch script (which will create symlink)
  if (!is.null(settings$model$jobtemplate) && file.exists(settings$model$jobtemplate)) {
    jobsh <- readLines(con = settings$model$jobtemplate, n = -1)
  } else {
    jobsh <- readLines(con = system.file("template.job", package = "PEcAn.SIPNET"), n = -1)
  }
  
  # create host specific setttings
  hostsetup <- ""
  if (!is.null(settings$model$prerun)) {
    hostsetup <- paste(hostsetup, sep = "\n", paste(settings$model$prerun, collapse = "\n"))
  }
  if (!is.null(settings$host$prerun)) {
    hostsetup <- paste(hostsetup, sep = "\n", paste(settings$host$prerun, collapse = "\n"))
  }
  
  # create cdo specific settings
  cdosetup <- ""
  if (!is.null(settings$host$cdosetup)) {
    cdosetup <- paste(cdosetup, sep = "\n", paste(settings$host$cdosetup, collapse = "\n"))
  }
  
  hostteardown <- ""
  if (!is.null(settings$model$postrun)) {
    hostteardown <- paste(hostteardown, sep = "\n", paste(settings$model$postrun, collapse = "\n"))
  }
  if (!is.null(settings$host$postrun)) {
    hostteardown <- paste(hostteardown, sep = "\n", paste(settings$host$postrun, collapse = "\n"))
  }
  
  # create rabbitmq specific setup.
  cpruncmd <- cpoutcmd <- rmoutdircmd <- rmrundircmd <- ""
  if (!is.null(settings$host$rabbitmq)) {
    #rsync cmd from remote to local host.
    settings$host$rabbitmq$cpfcmd <- ifelse(is.null(settings$host$rabbitmq$cpfcmd), "", settings$host$rabbitmq$cpfcmd)
    cpruncmd <- gsub("@OUTDIR@", settings$host$rundir, settings$host$rabbitmq$cpfcmd)
    cpruncmd <- gsub("@OUTFOLDER@", rundir, cpruncmd)
    
    cpoutcmd <- gsub("@OUTDIR@", settings$host$outdir, settings$host$rabbitmq$cpfcmd)
    cpoutcmd <- gsub("@OUTFOLDER@", outdir, cpoutcmd)
    
    #delete files within rundir and outdir.
    rmoutdircmd <- paste("rm", file.path(outdir, "*"))
    rmrundircmd <- paste("rm", file.path(rundir, "*"))
  }
  
  # create job.sh
  jobsh <- gsub("@HOST_SETUP@", hostsetup, jobsh)
  jobsh <- gsub("@CDO_SETUP@", cdosetup, jobsh)
  jobsh <- gsub("@HOST_TEARDOWN@", hostteardown, jobsh)
  
  jobsh <- gsub("@SITE_LAT@", settings$run$site$lat, jobsh)
  jobsh <- gsub("@SITE_LON@", settings$run$site$lon, jobsh)
  jobsh <- gsub("@SITE_MET@", template.clim, jobsh)
  
  jobsh <- gsub("@OUTDIR@", outdir, jobsh)
  jobsh <- gsub("@RUNDIR@", rundir, jobsh)
  
  jobsh <- gsub("@START_DATE@", settings$run$start.date, jobsh)
  jobsh <- gsub("@END_DATE@",settings$run$end.date , jobsh)
  
  jobsh <- gsub("@BINARY@", settings$model$binary, jobsh)
  jobsh <- gsub("@REVISION@", settings$model$revision, jobsh)
  
  jobsh <- gsub("@CPRUNCMD@", cpruncmd, jobsh)
  jobsh <- gsub("@CPOUTCMD@", cpoutcmd, jobsh)
  jobsh <- gsub("@RMOUTDIRCMD@", rmoutdircmd, jobsh)
  jobsh <- gsub("@RMRUNDIRCMD@", rmrundircmd, jobsh)
  
  if(is.null(settings$state.data.assimilation$NC.Prefix)){
    settings$state.data.assimilation$NC.Prefix <- "sipnet.out"
  }
  jobsh <- gsub("@PREFIX@", settings$state.data.assimilation$NC.Prefix, jobsh)
  
  #overwrite argument
  if(is.null(settings$state.data.assimilation$NC.Overwrite)){
    settings$state.data.assimilation$NC.Overwrite <- FALSE
  }
  jobsh <- gsub("@OVERWRITE@", settings$state.data.assimilation$NC.Overwrite, jobsh)
  
  #allow conflict? meaning allow full year nc export.
  if(is.null(settings$state.data.assimilation$FullYearNC)){
    settings$state.data.assimilation$FullYearNC <- FALSE
  }
  jobsh <- gsub("@CONFLICT@", settings$state.data.assimilation$FullYearNC, jobsh)
  
  if (is.null(settings$model$delete.raw)) {
    settings$model$delete.raw <- FALSE
  }
  jobsh <- gsub("@DELETE.RAW@", settings$model$delete.raw, jobsh)
  
  writeLines(jobsh, con = file.path(settings$rundir, run.id, "job.sh"))
  Sys.chmod(file.path(settings$rundir, run.id, "job.sh"))
  
  ### WRITE *.param-spatial
  template.paramSpatial <- system.file("template.param-spatial", package = "PEcAn.SIPNET")
  file.copy(template.paramSpatial, file.path(settings$rundir, run.id, "sipnet.param-spatial"))
  
  ### WRITE *.param
  template.param <- system.file("template.param", package = "PEcAn.SIPNET")
  if ("default.param" %in% names(settings$model)) {
    template.param <- settings$model$default.param
  }
  
  param <- utils::read.table(template.param)
  
  #### write run-specific PFT parameters here #### Get parameters being handled by PEcAn
  for (pft in seq_along(trait.values)) {
    pft.traits <- unlist(trait.values[[pft]])
    pft.names <- names(pft.traits)
    
    ## Append/replace params specified as constants
    constant.traits <- unlist(defaults[[1]]$constants)
    constant.names <- names(constant.traits)
    
    # Replace matches
    for (i in seq_along(constant.traits)) {
      ind <- match(constant.names[i], pft.names)
      if (is.na(ind)) {
        # Add to list
        pft.names <- c(pft.names, constant.names[i])
        pft.traits <- c(pft.traits, constant.traits[i])
      } else {
        # Replace existing value
        pft.traits[ind] <- constant.traits[i]
      }
    }
    
    # Remove NAs. Constants may be specified as NA to request template defaults. Note that it is 'NA'
    # (character) not actual NA due to being read in as XML
    pft.names <- pft.names[pft.traits != "NA" & !is.na(pft.traits)]
    pft.traits <- pft.traits[pft.traits != "NA" & !is.na(pft.traits)]
    pft.traits <- as.numeric(pft.traits)
    
    # Leaf carbon concentration
    leafC <- 0.48  #0.5
    if ("leafC" %in% pft.names) {
      leafC <- pft.traits[which(pft.names == "leafC")]
      id <- which(param[, 1] == "cFracLeaf")
      param[id, 2] <- leafC * 0.01  # convert to percentage from 0 to 1
    }
    
    # Specific leaf area converted to SLW
    SLA <- NA
    id <- which(param[, 1] == "leafCSpWt")
    if ("SLA" %in% pft.names) {
      SLA <- pft.traits[which(pft.names == "SLA")]
      param[id, 2] <- 1000 * leafC * 0.01 / SLA
    } else {
      SLA <- 1000 * leafC / param[id, 2]
    }
    
    # Maximum photosynthesis
    Amax <- NA
    id <- which(param[, 1] == "aMax")
    if ("Amax" %in% pft.names) {
      Amax <- pft.traits[which(pft.names == "Amax")]
      param[id, 2] <- Amax * SLA
    } else {
      Amax <- param[id, 2] * SLA
    }
    # Daily fraction of maximum photosynthesis
    if ("AmaxFrac" %in% pft.names) {
      param[which(param[, 1] == "aMaxFrac"), 2] <- pft.traits[which(pft.names == "AmaxFrac")]
    }
    
    ### Canopy extinction coefficiet (k)
    if ("extinction_coefficient" %in% pft.names) {
      param[which(param[, 1] == "attenuation"), 2] <- pft.traits[which(pft.names == "extinction_coefficient")]
    }
    
    # Leaf respiration rate converted to baseFolRespFrac
    if ("leaf_respiration_rate_m2" %in% pft.names) {
      Rd <- pft.traits[which(pft.names == "leaf_respiration_rate_m2")]
      id <- which(param[, 1] == "baseFolRespFrac")
      param[id, 2] <- max(min(Rd/Amax, 1), 0)
    }
    
    # Low temp threshold for photosynethsis
    if ("Vm_low_temp" %in% pft.names) {
      param[which(param[, 1] == "psnTMin"), 2] <- pft.traits[which(pft.names == "Vm_low_temp")]
    }
    
    # Opt. temp for photosynthesis
    if ("psnTOpt" %in% pft.names) {
      param[which(param[, 1] == "psnTOpt"), 2] <- pft.traits[which(pft.names == "psnTOpt")]
    }
    
    # Growth respiration factor (fraction of GPP)
    if ("growth_resp_factor" %in% pft.names) {
      param[which(param[, 1] == "growthRespFrac"), 2] <- pft.traits[which(pft.names == "growth_resp_factor")]
    }
    ### !!! NOT YET USED
    #Jmax = NA
    #if("Jmax" %in% pft.names){
    #  Jmax = pft.traits[which(pft.names == 'Jmax')]
    ### Using Jmax scaled to 25 degC. Maybe not be the best approach
    #}
    
    #alpha = NA
    #if("quantum_efficiency" %in% pft.names){
    #  alpha = pft.traits[which(pft.names == 'quantum_efficiency')]
    #}
    
    # Half saturation of PAR.  PAR at which photosynthesis occurs at 1/2 theoretical maximum (Einsteins * m^-2 ground area * day^-1).
    #if(!is.na(Jmax) & !is.na(alpha)){
    # param[which(param[,1] == "halfSatPar"),2] = Jmax/(2*alpha)
    ### WARNING: this is a very coarse linear approximation and needs improvement *****
    ### Yes, we also need to work on doing a paired query where we have both data together.
    ### Once halfSatPar is calculated, need to remove Jmax and quantum_efficiency from param list so they are not included in SA
    #}
    ### !!!
    
    # Half saturation of PAR.  PAR at which photosynthesis occurs at 1/2 theoretical maximum (Einsteins * m^-2 ground area * day^-1).
    # Temporary implementation until above is working.
    if ("half_saturation_PAR" %in% pft.names) {
      param[which(param[, 1] == "halfSatPar"), 2] <- pft.traits[which(pft.names == "half_saturation_PAR")]
    }
    
    # Ball-berry slomatal slope parameter m
    if ("stomatal_slope.BB" %in% pft.names) {
      id <- which(param[, 1] == "m_ballBerry")
      param[id, 2] <- pft.traits[which(pft.names == "stomatal_slope.BB")]
    }
    
    # Slope of VPD–photosynthesis relationship. dVpd = 1 - dVpdSlope * vpd^dVpdExp
    if ("dVPDSlope" %in% pft.names) {
      param[which(param[, 1] == "dVpdSlope"), 2] <- pft.traits[which(pft.names == "dVPDSlope")]
    }
    
    # VPD–water use efficiency relationship.  dVpd = 1 - dVpdSlope * vpd^dVpdExp
    if ("dVpdExp" %in% pft.names) {
      param[which(param[, 1] == "dVpdExp"), 2] <- pft.traits[which(pft.names == "dVpdExp")]
    }
    
    # Leaf turnover rate average turnover rate of leaves, in fraction per day NOTE: read in as
    # per-year rate!
    if ("leaf_turnover_rate" %in% pft.names) {
      param[which(param[, 1] == "leafTurnoverRate"), 2] <- pft.traits[which(pft.names == "leaf_turnover_rate")]
    }
    
    if ("wueConst" %in% pft.names) {
      param[which(param[, 1] == "wueConst"), 2] <- pft.traits[which(pft.names == "wueConst")]
    }
    # vegetation respiration Q10.
    if ("veg_respiration_Q10" %in% pft.names) {
      param[which(param[, 1] == "vegRespQ10"), 2] <- pft.traits[which(pft.names == "veg_respiration_Q10")]
    }
    
    # Base vegetation respiration. vegetation maintenance respiration at 0 degrees C (g C respired * g^-1 plant C * day^-1)
    # NOTE: only counts plant wood C - leaves handled elsewhere (both above and below-ground: assumed for now to have same resp. rate)
    # NOTE: read in as per-year rate!
    if ("stem_respiration_rate" %in% pft.names) {
      vegRespQ10 <- param[which(param[, 1] == "vegRespQ10"), 2]
      id <- which(param[, 1] == "baseVegResp")
      ## Convert from umols CO2 kg s-1 to gC g day-1
      stem_resp_g <- (((pft.traits[which(pft.names == "stem_respiration_rate")]) *
                         (44.0096 / 1e+06) * (12.01 / 44.0096)) / 1000) * 86400
      ## use Q10 to convert stem resp from reference of 25C to 0C param[id,2] =
      ## pft.traits[which(pft.names=='stem_respiration_rate')]*vegRespQ10^(-25/10)
      param[id, 2] <- stem_resp_g * vegRespQ10^(-25/10)
    }
    
    # turnover of fine roots (per year rate)
    if ("root_turnover_rate" %in% pft.names) {
      id <- which(param[, 1] == "fineRootTurnoverRate")
      param[id, 2] <- pft.traits[which(pft.names == "root_turnover_rate")]
    }
    
    # fine root respiration Q10
    if ("fine_root_respiration_Q10" %in% pft.names) {
      param[which(param[, 1] == "fineRootQ10"), 2] <- pft.traits[which(pft.names == "fine_root_respiration_Q10")]
    }
    
    # base respiration rate of fine roots (per year rate)
    if ("root_respiration_rate" %in% pft.names) {
      fineRootQ10 <- param[which(param[, 1] == "fineRootQ10"), 2]
      id <- which(param[, 1] == "baseFineRootResp")
      ## Convert from umols CO2 kg s-1 to gC g day-1
      root_resp_rate_g <- (((pft.traits[which(pft.names == "root_respiration_rate")]) *
                              (44.0096/1e+06) * (12.01 / 44.0096)) / 1000) * 86400
      ## use Q10 to convert stem resp from reference of 25C to 0C param[id,2] =
      ## pft.traits[which(pft.names=='root_respiration_rate')]*fineRootQ10^(-25/10)
      param[id, 2] <- root_resp_rate_g * fineRootQ10 ^ (-25 / 10)
    }
    
    # coarse root respiration Q10
    if ("coarse_root_respiration_Q10" %in% pft.names) {
      param[which(param[, 1] == "coarseRootQ10"), 2] <- pft.traits[which(pft.names == "coarse_root_respiration_Q10")]
    }
    # WARNING: fineRootAllocation + woodAllocation + leafAllocation isn't supposed to exceed 1
    # see sipnet.c code L2005 :
    # fluxes.coarseRootCreation=(1-params.leafAllocation-params.fineRootAllocation-params.woodAllocation)*npp;
    # priors can be chosen accordingly, and SIPNET doesn't really crash when sum>1 but better keep an eye
    alloc_params <- c("root_allocation_fraction", "wood_allocation_fraction", "leaf_allocation_fraction")
    if (all(alloc_params %in% pft.names)) {
      sum_alloc <- pft.traits[which(pft.names == "root_allocation_fraction")] +
        pft.traits[which(pft.names == "wood_allocation_fraction")] +
        pft.traits[which(pft.names == "leaf_allocation_fraction")]
      if(sum_alloc > 1){
        # I want this to be a severe for now, lateer can be changed back to warning
        PEcAn.logger::logger.warn("Sum of allocation parameters exceeds 1 for runid = ", run.id,
                                  "- This won't break anything since SIPNET has internal check, but notice that such combinations might not take effect in the outputs.")
      }
    }
    
    
    # fineRootAllocation
    if ("root_allocation_fraction" %in% pft.names) {
      param[which(param[, 1] == "fineRootAllocation"), 2] <- pft.traits[which(pft.names == "root_allocation_fraction")]
    }
    
    # woodAllocation
    if ("wood_allocation_fraction" %in% pft.names) {
      param[which(param[, 1] == "woodAllocation"), 2] <- pft.traits[which(pft.names == "wood_allocation_fraction")]
    }
    
    # leafAllocation
    if ("leaf_allocation_fraction" %in% pft.names) {
      param[which(param[, 1] == "leafAllocation"), 2] <- pft.traits[which(pft.names == "leaf_allocation_fraction")]
    }
    
    # wood_turnover_rate
    if ("wood_turnover_rate" %in% pft.names) {
      param[which(param[, 1] == "woodTurnoverRate"), 2] <- pft.traits[which(pft.names == "wood_turnover_rate")]
    }
    
    ### ----- Soil parameters soil respiration Q10.
    if ("soil_respiration_Q10" %in% pft.names) {
      param[which(param[, 1] == "soilRespQ10"), 2] <- pft.traits[which(pft.names == "soil_respiration_Q10")]
    }
    # soil respiration rate -- units = 1/year, reference = 0C
    if ("som_respiration_rate" %in% pft.names) {
      param[which(param[, 1] == "baseSoilResp"), 2] <- pft.traits[which(pft.names == "som_respiration_rate")]
    }
    
    # litterBreakdownRate
    if ("turn_over_time" %in% pft.names) {
      id <- which(param[, 1] == "litterBreakdownRate")
      param[id, 2] <- pft.traits[which(pft.names == "turn_over_time")]
    }
    # frozenSoilEff
    if ("frozenSoilEff" %in% pft.names) {
      param[which(param[, 1] == "frozenSoilEff"), 2] <- pft.traits[which(pft.names == "frozenSoilEff")]
    }
    
    # frozenSoilFolREff
    if ("frozenSoilFolREff" %in% pft.names) {
      param[which(param[, 1] == "frozenSoilFolREff"), 2] <- pft.traits[which(pft.names == "frozenSoilFolREff")]
    }
    
    # soilWHC
    if ("soilWHC" %in% pft.names) {
      param[which(param[, 1] == "soilWHC"), 2] <- pft.traits[which(pft.names == "soilWHC")]
    }
    # 10/31/2017 IF: these were the two assumptions used in the emulator paper in order to reduce dimensionality
    # These results in improved winter soil respiration values
    # they don't affect anything when the seasonal soil respiration functionality in SIPNET is turned-off
    if(TRUE){
      # assume soil resp Q10 cold == soil resp Q10
      param[which(param[, 1] == "soilRespQ10Cold"), 2] <- param[which(param[, 1] == "soilRespQ10"), 2]
      # default SIPNET prior of baseSoilRespCold was 1/4th of baseSoilResp
      # assuming they will scale accordingly
      param[which(param[, 1] == "baseSoilRespCold"), 2] <- param[which(param[, 1] == "baseSoilResp"), 2] * 0.25
    }
    
    if ("immedEvapFrac" %in% pft.names) {
      param[which(param[, 1] == "immedEvapFrac"), 2] <- pft.traits[which(pft.names == "immedEvapFrac")]
    }
    
    if ("leafWHC" %in% pft.names) {
      param[which(param[, 1] == "leafPoolDepth"), 2] <- pft.traits[which(pft.names == "leafWHC")]
    }
    
    if ("waterRemoveFrac" %in% pft.names) {
      param[which(param[, 1] == "waterRemoveFrac"), 2] <- pft.traits[which(pft.names == "waterRemoveFrac")]
    }
    
    if ("fastFlowFrac" %in% pft.names) {
      param[which(param[, 1] == "fastFlowFrac"), 2] <- pft.traits[which(pft.names == "fastFlowFrac")]
    }
    
    if ("rdConst" %in% pft.names) {
      param[which(param[, 1] == "rdConst"), 2] <- pft.traits[which(pft.names == "rdConst")]
    }
    ### ----- Phenology parameters GDD leaf on
    if ("GDD" %in% pft.names) {
      param[which(param[, 1] == "gddLeafOn"), 2] <- pft.traits[which(pft.names == "GDD")]
    }
    
    # Fraction of leaf fall per year (should be 1 for decid)
    if ("fracLeafFall" %in% pft.names) {
      param[which(param[, 1] == "fracLeafFall"), 2] <- pft.traits[which(pft.names == "fracLeafFall")]
    }
    
    # Leaf growth.  Amount of C added to the leaf during the greenup period
    if ("leafGrowth" %in% pft.names) {
      param[which(param[, 1] == "leafGrowth"), 2] <- pft.traits[which(pft.names == "leafGrowth")]
    }
    
    #update LeafOnday and LeafOffDay
    if (!is.null(settings$run$inputs$leaf_phenology)){
      obs_year_start <- lubridate::year(settings$run$start.date)
      obs_year_end <- lubridate::year(settings$run$end.date)
      if (obs_year_start != obs_year_end) {
        PEcAn.logger::logger.info("Start.date and end.date are not in the same year. Currently start.date is used for refering phenological data")
      }
      leaf_pheno_path <- settings$run$inputs$leaf_phenology$path  ## read from settings
      if (!is.null(leaf_pheno_path)){
        ##read data
        leafphdata <- utils::read.csv(leaf_pheno_path)
        leafOnDay <-  leafphdata$leafonday[leafphdata$year == obs_year_start & leafphdata$site_id==settings$run$site$id]
        leafOffDay<-  leafphdata$leafoffday[leafphdata$year== obs_year_start & leafphdata$site_id==settings$run$site$id]
        if (!is.na(leafOnDay)){
          param[which(param[, 1] == "leafOnDay"), 2] <- leafOnDay
        }
        if (!is.na(leafOffDay)){
          param[which(param[, 1] == "leafOffDay"), 2] <- leafOffDay
        }
      } else {
        PEcAn.logger::logger.info("No phenology data were found. Please consider running `PEcAn.data.remote::extract_phenology_MODIS` to get the parameter file.")
      }
    }
  } ## end loop over PFTS
  ####### end parameter update
  #working on reading soil file (only working for 1 soil file)
  if(length(settings$run$inputs$soilinitcond$path)==1){
    soil_IC_list <- PEcAn.data.land::pool_ic_netcdf2list(settings$run$inputs$soilinitcond$path)
    #SoilWHC and LitterWHC
    if("volume_fraction_of_water_in_soil_at_saturation"%in%names(soil_IC_list$vals)){
      #SoilWHC
      param[which(param[, 1] == "soilWHC"), 2] <- mean(unlist(soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"]))*100
      
      #LitterWHC
      #param[which(param[, 1] == "litterWHC"), 2] <- unlist(soil_IC_list$vals["volume_fraction_of_water_in_soil_at_saturation"])[1]*100
    }
    if("soil_hydraulic_conductivity_at_saturation"%in%names(soil_IC_list$vals)){
      #litwaterDrainrate
      param[which(param[, 1] == "litWaterDrainRate"), 2] <- unlist(soil_IC_list$vals["soil_hydraulic_conductivity_at_saturation"])[1]*100/(3600*24)
    }
  }
  if (!is.null(IC)) {
    ic.names <- names(IC)
    ## plantWoodInit gC/m2
    plant_wood_vars <- c("AbvGrndWood", "abvGrndWoodFrac", "coarseRootFrac", "fineRootFrac")
    if (all(plant_wood_vars %in% ic.names)) {
      # reconstruct total wood C
      if(IC$abvGrndWoodFrac < 0.05){
        wood_total_C <- IC$AbvGrndWood
      }else{
        wood_total_C <- IC$AbvGrndWood / IC$abvGrndWoodFrac
      }
      
      #Sanity check
      if (is.infinite(wood_total_C) | is.nan(wood_total_C) | wood_total_C < 0) {
        wood_total_C <- 0
        if (round(IC$AbvGrndWood) > 0 & round(IC$abvGrndWoodFrac, 3) == 0)
          PEcAn.logger::logger.warn(
            paste0(
              "There is a major problem with ",
              run.id,
              " in either the model's parameters or IC.",
              "Because the ABG is estimated=",
              IC$AbvGrndWood,
              " while AGB Frac is estimated=",
              IC$abvGrndWoodFrac
            )
          )
      }
      param[which(param[, 1] == "plantWoodInit"),  2] <- wood_total_C
      param[which(param[, 1] == "coarseRootFrac"), 2] <- IC$coarseRootFrac
      param[which(param[, 1] == "fineRootFrac"),   2] <- IC$fineRootFrac
    }
    ## laiInit m2/m2
    if ("lai" %in% ic.names) {
      param[which(param[, 1] == "laiInit"), 2] <- IC$lai
    }
    ## litterInit gC/m2
    if ("litter_carbon_content" %in% ic.names) {
      param[which(param[, 1] == "litterInit"), 2] <- IC$litter_carbon_content
    }
    ## soilInit gC/m2
    if ("soil" %in% ic.names) {
      param[which(param[, 1] == "soilInit"), 2] <- IC$soil
    }
    ## litterWFracInit fraction
    if ("litter_mass_content_of_water" %in% ic.names) {
      #here we use litterWaterContent/litterWHC to calculate the litterWFracInit
      param[which(param[, 1] == "litterWFracInit"), 2] <- IC$litter_mass_content_of_water/(param[which(param[, 1] == "litterWHC"), 2]*10)
    }
    ## soilWFracInit fraction
    if ("soilWFrac" %in% ic.names) {
      param[which(param[, 1] == "soilWFracInit"), 2] <- IC$soilWFrac
    }
    ## snowInit cm water equivalent
    if ("SWE" %in% ic.names) {
      param[which(param[, 1] == "snowInit"), 2] <- IC$SWE
    }
    ## microbeInit mgC/g soil
    if ("microbe" %in% ic.names) {
      param[which(param[, 1] == "microbeInit"), 2] <- IC$microbe
    }
  }
  
  else if (length(settings$run$inputs$poolinitcond$path)>0) {
    ICs_num <- length(settings$run$inputs$poolinitcond$path)
    IC.path <- settings$run$inputs$poolinitcond$path[[sample(1:ICs_num, 1)]]
    
    IC.pools <- PEcAn.data.land::prepare_pools(IC.path, constants = list(sla = SLA))
    
    if(!is.null(IC.pools)){
      IC.nc <- ncdf4::nc_open(IC.path) #for additional variables specific to SIPNET
      ## plantWoodInit gC/m2
      if ("wood" %in% names(IC.pools)) {
        param[which(param[, 1] == "plantWoodInit"), 2] <- PEcAn.utils::ud_convert(IC.pools$wood, "kg m-2", "g m-2")
      }
      ## laiInit m2/m2
      lai <- IC.pools$LAI
      if (!is.na(lai) && is.numeric(lai)) {
        param[which(param[, 1] == "laiInit"), 2] <- lai
      }
      
      #Initial LAI is set as 0 for deciduous forests and grasslands for non-growing seasons
      if (!(lubridate::month(settings$run$start.date) %in% seq(5,9))){ #Growing seasons are coarsely defined as months from May to September for non-conifers in the US
        site_pft <- utils::read.csv(settings$run$inputs$pft.site$path)
        site.pft.name <- site_pft$pft[site_pft$site == settings$run$site$id]
        if (site.pft.name!="boreal.coniferous") {   #Currently only excluding boreal conifers. Other evergreen PFTs could be added here later.
          param[which(param[, 1] == "laiInit"), 2] <- 0       
        }
      }
      ## neeInit gC/m2
      nee <- try(ncdf4::ncvar_get(IC.nc,"nee"),silent = TRUE)
      if (!is.na(nee) && is.numeric(nee)) {
        param[which(param[, 1] == "neeInit"), 2] <- nee
      }
      ## litterInit gC/m2
      if ("litter" %in% names(IC.pools)) {
        param[which(param[, 1] == "litterInit"), 2] <- PEcAn.utils::ud_convert(IC.pools$litter, 'g m-2', 'g m-2') # BETY: kgC m-2
      }
      ## soilInit gC/m2
      if ("soil" %in% names(IC.pools)) {
        param[which(param[, 1] == "soilInit"), 2] <- PEcAn.utils::ud_convert(sum(IC.pools$soil), 'kg m-2', 'g m-2') # BETY: kgC m-2
      }
      ## soilWFracInit fraction
      soilWFrac <- try(ncdf4::ncvar_get(IC.nc,"SoilMoistFrac"),silent = TRUE)
      if (!"try-error" %in% class(soilWFrac)) {
        if (!is.na(soilWFrac) && is.numeric(soilWFrac)) {
          param[which(param[, 1] == "soilWFracInit"), 2] <- sum(soilWFrac)/100
        }
      }
      ## litterWFracInit fraction
      litterWFrac <- soilWFrac
      
      ## snowInit cm water equivalent (cm = g / cm2 because 1 g water = 1 cm3 water)
      snow = try(ncdf4::ncvar_get(IC.nc,"SWE"),silent = TRUE)
      if (!is.na(snow) && is.numeric(snow)) {
        param[which(param[, 1] == "snowInit"), 2] <- PEcAn.utils::ud_convert(snow, "kg m-2", "g cm-2")  # BETY: kg m-2
      }
      ## leafOnDay
      leafOnDay <- try(ncdf4::ncvar_get(IC.nc,"date_of_budburst"),silent = TRUE)
      if (!is.na(leafOnDay) && is.numeric(leafOnDay)) {
        param[which(param[, 1] == "leafOnDay"), 2] <- leafOnDay
      }
      ## leafOffDay
      leafOffDay <- try(ncdf4::ncvar_get(IC.nc,"date_of_senescence"),silent = TRUE)
      if (!is.na(leafOffDay) && is.numeric(leafOffDay)) {
        param[which(param[, 1] == "leafOffDay"), 2] <- leafOffDay
      }
      microbe <- try(ncdf4::ncvar_get(IC.nc,"Microbial Biomass C"),silent = TRUE)
      if (!is.na(microbe) && is.numeric(microbe)) {
        param[which(param[, 1] == "microbeInit"), 2] <- PEcAn.utils::ud_convert(microbe, "mg kg-1", "mg g-1") #BETY: mg microbial C kg-1 soil
      }
      
      ncdf4::nc_close(IC.nc)
    }else{
      PEcAn.logger::logger.error("Bad initial conditions filepath; keeping defaults")
    }
  }else{
    #some stuff about IC file that we can give in lieu of actual ICs
  }
  
  
  if (!is.null(settings$run$inputs$soilmoisture)) {
    #read soil moisture netcdf file, grab closet date to start_date, set equal to soilWFrac
    if(!is.null(settings$run$inputs$soilmoisture$path)){
      soil.path <- settings$run$inputs$soilmoisture$path
      soilWFrac <- ncdf4::ncvar_get(ncdf4::nc_open(soil.path), varid = "mass_fraction_of_unfrozen_water_in_soil_moisture")
      
      param[which(param[, 1] == "soilWFracInit"), 2] <- soilWFrac
    }
    
  }
  if(file.exists(file.path(settings$rundir, run.id, "sipnet.param"))) file.rename(file.path(settings$rundir, run.id, "sipnet.param"),file.path(settings$rundir, run.id, paste0("sipnet_",lubridate::year(settings$run$start.date),"_",lubridate::year(settings$run$end.date),".param")))
  
  
  utils::write.table(param, file.path(settings$rundir, run.id, "sipnet.param"), row.names = FALSE, col.names = FALSE,
                     quote = FALSE)
}




