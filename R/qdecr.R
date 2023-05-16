#' Run vertex-wise analyses
#' 
#' Runs vertex-wise analyses using a variety of statistical models
#' 
#' This function is the worker function of the \code{QDECR} package. It allows
#' .mgh format data as input and allows statistical analyses per vertex. A variety
#' of statistical models have been implemented, such as linear regression.
#' 
#' @param id the name of the id variable that matches the dataset to the Freesurfer output
#' @param data a required argument that contains a data frame, a list of data frames or an imputed object that is supported by the `imp2list` function (mice, mi, etc.).
#' @param vertex the preposition to the vertex measure (default: "qdecr_")
#' @param margs the arguments that should be provided to the function of interest (e.g. stats::lm)
#' @param model the function to grab the arguments from (this will be removed in a later version)
#' @param target the target template (default = "fsaverage")
#' @param hemi hemisphere to analyze ("lh" or "rh")
#' @param measure the vertex-wise measure to use ("thickness", "area", etc.) as specified in `measure_choices`
#' @param measure_choices the possible measures that can be chosen from
#' @param fwhm full width half max (default = 10 mm, for pial_lgi it is 5 mm)
#' @param mcz_thr A numeric value for the Monte Carlo simulation threshold (default: 0.001). Any of the following are accepted (equivalent values separate by `/`): 13/1.3/0.05, 20/2.0/0.01, 23/2.3/0.005, 30/3.0/0.001, 33/3.3/0.0005, 40/4.0/0.0001.
#' @param cwp_thr the cluster-wise p-value threshold on top of all correction (default = 0.025, as there are 2 hemispheres)
#' @param mgh NOT IMPLEMENTED; path to existing merged mgh file, default is NULL
#' @param mask mgh file to mask analysis; default is to use the cortex label from the target
#' @param mask_path path to the mask; default is the cortex mask that is provided with the QDECR package
#' @param project the base name you want to assign to the output files
#' @param dir_subj directory contain the surface-based maps (mgh files); defaults to SUBJECTS_DIR
#' @param dir_fshome Freesurfer directory; defaults to FREESURFER_HOME
#' @param dir_tmp directory to store the temporary big matrices; useful for shared memory; defaults to `dir_out`
#' @param dir_out the directory where to save the data to (defaults to the current directory)
#' @param dir_out_tree if TRUE, creates a dir_out/project directory. If FALSE, all output is placed directory into dir_out
#' @param file_out_tree if TRUE, adds the full project name to the output file names. By default it is the inverse of dir_out_tree
#' @param clean_up NOT IMPLEMENTED; will be used for setting cleaning of other files
#' @param clean_up_bm if TRUE, cleans all big matrices (.bk) that were generated in dir_tmp
#' @param clobber if TRUE, ignores already existing directories and writes over them; if FALSE, stops and warns user that a given directory already exists
#' @param verbose if TRUE, writes out standard log; if FALSE, no output is generated
#' @param debug NOT IMPLEMENTED; will output the maximal log to allow for easy debugging
#' @param n_cores the number of cores to be used
#' @param prep_fun Name of the function that needs to be called for the preparation step (do not touch unless you know what you are doing!)
#' @param analysis_fun Name of the function that needs to be called for the analysis step (do not touch unless you know what you are doing!)
#' @param chunk_size Integer; the desired chunk size for the chunked lm
#' 
#' @return out
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
# @importFrom methods formalArgs

qdecr <- function(id,
                  data,
                  vertex = "qdecr",
                  margs = NULL,
                  model = c("stats::lm", "QDECR::megha", "RcppEigen::fastLm", 
                            "stats::glm", "survival::coxph", "default"),
                  target = "fsaverage", 
                  hemi = c("lh", "rh"),
                  measure,
                  measure_choices = c("thickness", "area", "area.pial", "curv", 
                                      "jacobian_white", "pial", "pial_lgi", "sulc", "volume", 
                                      "w_g.pct", "white.H", "white.K"),
                  fwhm = ifelse(measure == "pial_lgi", 5, 10),
                  mcz_thr = 0.001, 
                  cwp_thr = 0.025,
                  mgh = NULL,
                  mask = NULL,
                  mask_path = system.file("extdata", paste0(hemi, ".fsaverage.cortex.mask.mgh"), package = "QDECR"),
                  project,
                  dir_subj = Sys.getenv("SUBJECTS_DIR"),
                  dir_fshome = Sys.getenv("FREESURFER_HOME"),
                  dir_tmp = dir_out,
                  dir_out,
                  dir_out_tree = TRUE,
                  file_out_tree = !dir_out_tree, 
                  clean_up = TRUE,
                  clean_up_bm = TRUE,
                  clobber = FALSE,
                  verbose = TRUE,
                  debug = FALSE,
                  n_cores = 1,
                  prep_fun = "prep_fastlm",
                  analysis_fun = "analysis_chunkedlm",
                  chunk_size = 1000
                  ){

  if (verbose) {
    message("\n")
    catprompt(paste0("QDECR v", utils::packageVersion("QDECR")))
    message(paste0("Welcome to QDECR (v", utils::packageVersion("QDECR"), ")"))
    message("Authors: Sander Lamballais & Ryan Muetzel")
    message("Website: www.qdecr.com")
    message("Repository: www.github.com/slamballais/QDECR")
  }
  
  ## Part 1: All checks and data prep
  catprompt("Integrity checks", verbose = verbose)

  vw <- list()
  vw$qdecr_call <- match.call()
  measure <- match.arg(measure, choices = measure_choices)
  hemi <- match.arg(hemi)
  model <- match.arg(model)

  # Check margs
  if (is.null(margs)) stop("`margs` is not defined, so no model can be run.")
  margs <- if (model == "default") qdecr_setnames(margs, deparse(margs[[1]])) else qdecr_setnames(margs, model)
  
  # Check data
  if(missing(data)){ 
    data <- qdecr_extract(margs, model)
    if (identical(data, "fail")) {
      stop("No input provided to the `data` argument. \n", 
           "Additionally, no data found in margs.")
    }
  }
  md <- imp2list(data)
  
  # Assemble and check input arguments
  vw$input <- qdecr_check(id, md, margs, hemi, vertex, measure, model, target, project, dir_out_tree, clobber, fwhm, n_cores, prep_fun, file_out_tree)
  vw$mask <- qdecr_check_mask(mask, mask_path)
  vw$paths <- check_paths(vw, dir_tmp, dir_subj, dir_out, dir_fshome, mask_path)
  mcz_thr <- check_mcz_thr(mcz_thr)

  # Check model
  vw$model <- qdecr_model(vw$input$model, vw$input$prep_fun, vw$input$md, vw$input$id, vw$input$vertex, vw$input$margs, vw$paths$dir_tmp2)
  
  # Check backing
  qdecr_check_backing(c(vw$model$backing, vw$model$backing_to_remove, vw$paths$backing_mgh), clobber)
  if (clean_up_bm){
    on.exit({
      for (n in c(vw$model$backing, vw$model$backing_to_remove, vw$paths$backing_mgh)) {
        if (file.exists(paste0(n, ".bk"))) file.remove(paste0(n, ".bk"))
      }
      vw$mgh <- NULL
    })
  }

  message2("Checked all input arguments.", verbose = verbose)

  ## Part 2: Load vertex data
  catprompt("Loading vertex data", verbose = verbose)
  vw$input$mgh_ids <- as.character(unlist(vw$input$md[[1]][,vw$input$id, drop = TRUE]))
  
  vw$mgh <- qdecr_prep_mgh(vw$paths$dir_subj, vw$input$mgh_ids, vw$mask, 
                           vw$input$hemi, vw$input$measure, vw$input$fwhmc, vw$input$target, 
                           vw$input$n_cores, vw$paths$dir_tmp, vw$input$project, 
                           backing = vw$paths$backing_mgh, verbose = verbose)
  message2("\n", verbose = verbose)
  
  ## Part 3: Make descriptions, print them if verbose
  catprompt("Pre-analysis overview", verbose = verbose)
  vw$describe <- qdecr_pre_describe(vw, mask, verbose = verbose)
  
  ## Part 4: Analyses
  catprompt("Vertex-wise analyses", verbose = verbose)
  class(vw) <- class(vw$model)
  vw$results <- vertexwise(vw, analysis_fun = analysis_fun, chunk = chunk_size)
  message2("\n\n", verbose = verbose)

  ## Part 5: Post-processing
  vw$post <- list()
  
  # Some stack stuff
  vw$stack <- qdecr_make_stack(vw, vw$model$so[1:4], mcz_thr)

  # Split all coefficients into separate files and save them in the final directory
  message2("Converting coefficients, p-values, etc to .mgh format", verbose = verbose)
  qdecr_move(vw$results[1:4], vw$stack, vw$model$so[1:4])

  # Save out the residuals
  message2("Converting residual data to .mgh format", verbose = verbose)
  vw$post$eres <- sub(".bk$", ".mgh", vw$results$residuals$backingfile)
  bsfbm2mgh(vw$results$residuals, vw$post$eres)
  if (clean_up_bm){
    on.exit(if (file.exists(vw$post$eres)) file.remove(vw$post$eres),
            add = TRUE)
  }

  catprompt("Estimating smoothness", verbose = verbose)

  # Estimate smoothness
  message2("Calculating smoothing for multiple testing correction.", 
           verbose = verbose)
  vw$post$fwhm_est <- calc_fwhm(vw$paths$final_path, vw$paths$final_mask_path, vw$paths$est_fwhm_path, vw$input$hemi, vw$post$eres, mask = vw$mask, verbose = verbose)[,]
  if(vw$post$fwhm_est > 30) {
    message2(paste0("Estimated smoothness is ", vw$post$fwhm_est, ", which is really high. Reduced to 30."), verbose = verbose)
    vw$post$fwhm_est <- 30
  }
  if(vw$post$fwhm_est < 1) {
    message2(paste0("Estimated smoothness is ", vw$post$fwhm_est, ", which is really low. Increased to 1."), verbose = verbose)
    vw$post$fwhm_est <- 1
  }

  catprompt("Clusterwise correction", verbose = verbose)

  # Run clusterwise and move all files to the right directory
  for ( i in seq_along(vw$stack$names)){
    message2("\n", verbose = verbose)
    run_mri_surf_cluster(vw, vw$stack$p[[i]], vw$post$fwhm_est, mask_path = vw$paths$final_mask_path, 
                      verbose = verbose, mcz_thr = mcz_thr, cwp_thr = cwp_thr,
                      stack = i, stack_name = vw$stack$names[i])
  }
  vw$post$final_mask <- load.mgh(vw$paths$final_mask_path)$x
  
  catprompt("Calculating vertex-wise values", verbose = verbose)

  # Do FBM-calculations for plotting histograms and such
  vw$post$mgh_description <- list()
  vw$post$mgh_description$vertex_mean <- fbm_row_mean(vw$mgh, n_cores, row.mask = vw$post$final_mask)
  vw$post$mgh_description$vertex_sd <- fbm_row_sd(vw$mgh, n_cores, row.mask = vw$post$final_mask)
  vw$post$mgh_description$subject_mean <- fbm_col_mean(vw$mgh, n_cores, row.mask = vw$post$final_mask)
  vw$post$mgh_description$subject_sd <- fbm_col_sd(vw$mgh, n_cores, row.mask = vw$post$final_mask)
  
  catprompt("Summarizing results", verbose = verbose)
  
  # Post-describe
  vw$describe <- qdecr_post_describe(vw, verbose = verbose)
  vw$summary <- summary(vw, verbose = verbose)
  if(clean_up_bm) vw$mgh <- NULL
  
  return(vw)
}
