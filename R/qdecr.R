#' Run vertex-wise analyses
#' 
#' Runs vertex-wise analyses using a variety of statistical models
#' 
#' This function is the worker function of the \code{QDECR} package. It allows
#' .mgh format data as input and allows statistical analyses per vertex. A variety
#' of statistical models have been implemented, such as linear regression and 
#' GCTA. 
#' 
#' @param a test
#' @param b test
#' @param c test
#' @param d test
#' @param e test
#' @param f test
#' @param g test
#' @param h test
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
                  measure = c("thickness", "area", "area.pial", "curv", 
                              "jacobian_white", "pial", "pial_lgi", "sulc", "volume", 
                              "w-g.pct.mgh", "white.H", "white.K"),
                  fwhm = ifelse(measure == "pial_lgi", 5, 10),
                  mcz_thr = 30, 
                  mgh = NULL,
                  mask = NULL,
                  mask_path = system.time("extdata", paste0(hemi, ".fsaverage.cortex.mask.mgh"), package = "QDECR"),
                  project,
                  dir_subj = Sys.getenv("SUBJECTS_DIR"),
                  dir_fshome = Sys.getenv("FREESURFER_HOME"),
                  dir_tmp = dir_out,
                  dir_out,
                  dir_out_tree = TRUE,
                  clean_up = TRUE,
                  clean_up_bm = TRUE,
                  clobber = FALSE,
                  verbose = TRUE,
                  debug = FALSE,
                  n_cores = 1,
                  prep_fun = "prep_fastlm",
                  analysis_fun = "analysis_chunkedlm",
                  chunk_size = 1000,
                  sander = FALSE
                  ){
  
  ### GENERAL STRUCTURE: All steps of qdecr are split into specific functions.
  ### Each function takes a list (of arguments, data, etc.) and outputs a modified
  ### list (of arguments, data, etc.) that can act as input for the next step.
  ### A massive further improvement is the `vw` object, which can take ANY
  ### function and can apply it to all vertices. This allows the user to apply
  ### ANY statistical model, vertex-wise. Other improvements:
  ### * qdecr calls without relying on the BigMemory package
  ### * a more accurate progress bar
  ### * more expansive contrast options
  ### * specifying vertex data as outcome or determinant by using a reserved name

  # Intro
  if (verbose) {
    message("\n")
    catprompt("QDECR v2.0")
    message("Welcome to QDECR (v2.0)")
    message("Authors: Ryan Muetzel & Sander Lamballais")
    message("Website: www.github.com/slamballais/QDECR")
  }
  
  ## Part 1: All checks and data prep
  catprompt("Integrity checks", verbose = verbose)
  
  # Create a vw object (this step should be more properly formalized with S4 objects
  vw <- list()
  
  # Check call
  vw$qdecr_call <- match.call()
  measure <- match.arg(measure)
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
  vw$input <- qdecr_check(id, md, margs, hemi, vertex, measure, model, target, project, dir_out_tree, clobber, fwhm, n_cores)
  vw$mask <- qdecr_check_mask(mask, mask_path)
  vw$paths <- check_paths(vw, dir_tmp, dir_subj, dir_out, dir_fshome, mask_path)

  # Check model
  vw$model <- qdecr_model(vw$input$model, vw$input$md, vw$input$id, vw$input$vertex, vw$input$margs, vw$paths$dir_tmp2)
  
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

  catprompt("Clusterwise correction", verbose = verbose)

  # Run clusterwise and move all files to the right directory
  for ( i in seq_along(vw$stack$names)){
    message2("\n", verbose = verbose)
    runMriSurfCluster(vw$paths$final_path, vw$paths$dir_fshome, vw$input$hemi, vw$stack$p[[i]], vw$post$fwhm_est, 
                      mask_path = vw$paths$final_mask_path, verbose = verbose, mczThr = mcz_thr, 
                      stack = i, stack_name = vw$stack$names[i])
  }
  vw$post$final_mask <- load.mgh(vw$paths$final_mask_path)$x
  
    catprompt("Calculating vertex-wise values", verbose = verbose)
  
  ### FINAL THINGS TO DO:
  # Do FBM-calculations for plotting histograms and such
  vw$post$mgh_description <- list()
  vw$post$mgh_description$vertex_mean <- fbm_row_mean(vw$mgh, n_cores, row.mask = vw$post$final_mask)
  vw$post$mgh_description$vertex_sd <- fbm_row_sd(vw$mgh, n_cores, row.mask = vw$post$final_mask)
  vw$post$mgh_description$subject_mean <- fbm_col_mean(vw$mgh, n_cores, row.mask = vw$post$final_mask)
  vw$post$mgh_description$subject_sd <- fbm_col_sd(vw$mgh, n_cores, row.mask = vw$post$final_mask)
  
    catprompt("Summarizing results", verbose = verbose)
  
  # Post-describe
  if (sander) return(vw)
  vw$describe <- qdecr_post_describe(vw, verbose = verbose)
  vw$summary <- summary(vw, verbose = verbose)
  if(clean_up_bm) vw$mgh <- NULL
  
  return(vw)
  
  
  # Save out the .rds files
  # Make JSON-like output file that contains all info
  
  # Move all files to the proper output structure
  
  
  # Clean everything
  qdecr_clean()
  
  return(qout)
  
}


  