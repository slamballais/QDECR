## COPIED FROM bigparallelr:::default_nproc_blas
default_nproc_blas <- function() {

  cl <- parallel::makePSOCKcluster(1)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, RhpcBLASctl::blas_get_num_procs())[[1]]
}

qdecr_check_backing <- function(backing, clobber){
  for (n in backing){
    n <- paste0(n, ".bk")
    if (file.exists(n)) {
      if (!clobber) {
        stop("File ", n, " already exists, and `clobber` is set to FALSE. \n",
             "Set `clobber` to TRUE to overwrite the backing file.")
      }
      file.remove(n)
    }
  }
}

qdecr_check <- function(id, md, margs, hemi, vertex, measure, model, target,
                        project, dir_out_tree, clobber, fwhm, n_cores, prep_fun,
                        file_out_tree){

  # verify that id is correct and exists within data
  check_id(id, md)

  # verify that vertex is correct, is not id and is not in data
  check_vertex(vertex, id, md)

  # Verify the integrity of fwhm
  check_fwhm(fwhm)

  # Check whether the number of cores is correct
  check_cores(n_cores)
  
  # Check whether prep_fun exists
  check_prep_fun(prep_fun)
  
  # make sure `dir_out_tree` and `clobber` do not interfere
  if (!dir_out_tree && clobber) stop("PREVENTATIVE STOP! `clobber` is set to TRUE while ", 
                                     "`dir_out_tree` is set to FALSE. This can have disastrous ",
                                     "consequences if used incorrectly. As a preventative ",
                                     "measure, we do not allow this to be run.")

  # Return all output arguments assembled
  input <- list()
  input[["id"]] <- id
  input[["md"]] <- md
  input[["margs"]] <- margs
  input[["hemi"]] <- hemi
  input[["vertex"]] <- vertex
  input[["measure"]] <- measure
  input[["model"]] <- model
  input[["target"]] <- target
  input[["project"]] <- project
  input[["project2"]] <- paste(hemi, project, measure, sep = ".")
  input[["dir_out_tree"]] <- dir_out_tree
  input[["file_out_tree"]] <- file_out_tree
  input[["clobber"]] <- clobber
  input[["fwhm"]] <- fwhm
  input[["fwhmc"]] <- if (fwhm > 0) paste0("fwhm", fwhm) else ""
  input[["n_cores"]] <- n_cores
  input[["prep_fun"]] <- prep_fun
  input
}

check_cores <- function(n_cores) {
  if (n_cores < 1) stop("You specified `n_cores` to be smaller than 1. Please choose a number that is 1 or higher.")
  rec_cores <- parallel::detectCores() - 1
  if (rec_cores == 0) rec_cores <- 1
  if (n_cores > rec_cores) stop("You specified `n_cores` to be too high (", n_cores, "). Recommended is ", rec_cores)
  if (n_cores > 1 && default_nproc_blas() > 1) stop("`n_cores` > 1, but there already seems to be a parallel BLAS library present. Either set `n_cores` to 1, or set the BLAS library to 1.")
  NULL
}

check_dir_out <- function(dir_out, project, project2, dir_out_tree, clobber){
  if (!dir.exists(dir_out)) {
      if(!dir.exists(dirname(dir_out)))
        stop("The provided `dir_out` does not exist, but neither does it parent.")
      dir.create(dir_out)
  } else {
    if(!clobber && !dir_out_tree)
      stop("`dir_out` already exists, and both `clobber` and `dir_out_tree` are FALSE. \n",
           "Either set `clobber` = TRUE to overwrite the output directory, \n",
           "or set `dir_out_tree` = TRUE and provide `project` to \n",
           "create a new subdirectory.")
    if(clobber && !dir_out_tree){
      stop("PREVENTATIVE STOP! `clobber` is set to TRUE while ", 
           "`dir_out_tree` is set to FALSE. This can have disastrous ",
           "consequences if used incorrectly. As a preventative ",
           "measure, we do not allow this to be run.")
    }
  }
  if(dir_out_tree){
    if(missing(project))
      stop("`dir_out_tree` is set to TRUE but no `project` provided.")
    p <- file.path(dir_out, project2)
    if(dir.exists(p)) {
      if (!clobber)
        stop("The output directory ", p, " already exists and `clobber` = FALSE. \n",
             "Set `clobber` = TRUE to overwrite the existing output directory.")
      unlink(p, recursive = TRUE)
    }
    dir.create(p)
    dir_out <- p
  }
  dir_out
}

check_dir_subj <- function(dir_subj, target, md, id){
  if(missing(dir_subj)) stop("No `dir_subj` argument provided. \n",
                             "Please provide a path to the directory that contains ",
                             "all the vertex input data.")
  if(!dir.exists(dir_subj))
    stop("The provided `dir_subj` does not seem to exist.")
  if(!dir.exists(target)){
    if(!dir.exists(file.path(dir_subj, target)))
      stop("The provided `target` directory does not seem to exist.")
  }
  idx <- md[[1]][[id]]
  f <- list.files(dir_subj)
  q <- !idx %in% f
  if(any(q)){
    qq <- sum(q, na.rm = TRUE)
    if (qq > 10){
      stop(qq, " ids in the `", id, "` column of the data are not present in `dir_subj`.")
    } else {
      stop("The following `id` values were not found in `dir_subj`: ",
           paste(idx[q], collapse = ", "),
           ".")
    }
  }
  NULL
}

check_fwhm <- function(fwhm){
  if(!is.numeric(fwhm)) stop("The provided `fwhm` is not numeric.")
  if(fwhm < 0) stop("The provided `fwhm` is a negative number.")
  NULL
}

check_id <- function(id, md){
  if(missing(id)) stop("No input provided to the `id` argument. \n",
                       "Please provide the column name that contains the subject names.")
  q <- sapply(md, function(x) !id %in% names(x))
  if (any(q)) {
    if (length(md) == 1) {
      stop("The provided id column name is not present in the dataset.")
    } else if (sum(q) == 1) {
      stop("The provided id column name is not present in dataset: ",
           paste(which(q), collapse = ", "),
           ".")
    } else {
      stop("The provided id column name is not present in datasets: ",
           paste(which(q), collapse = ", "),
           ".")
    }
  }
  if (length(unique(md[[1]][,id, drop = TRUE])) == 1)
    stop("The provided `id` column only has 1 unique value in your data, suggesting only 1 person is loaded in.")
  NULL
}

check_mcz_thr <- function(mcz_thr) {
  if (!is.numeric(mcz_thr)) stop("`mcz_thr` is not numeric.")
  if (mcz_thr %in% c(13, 1.3, 0.05)) 13 else 
    if (mcz_thr %in% c(20, 2.0, 0.01)) 20 else 
      if (mcz_thr %in% c(23, 2.3, 0.005)) 23 else 
        if (mcz_thr %in% c(30, 3.0, 0.001)) 30 else 
          if (mcz_thr %in% c(33, 3.3, 0.0005)) 33 else 
            if (mcz_thr %in% c(40, 4.0, 0.0001)) 40 else 
              stop("`mcz_thr` does not have an accepted value (13/1.3/0.05, 20/2.0/0.01, 23/2.3/0.005, 30/3.0/0.001, 33/3.3/0.0005, 40/4.0/0.0001).")
}

check_paths <- function(vw, dir_tmp, dir_subj, dir_out, dir_fshome, mask_path){
 check_dir_subj(dir_subj, vw$input$target, vw$input$md, vw$input$id)
 dir_out2 <- check_dir_out(dir_out, vw$input$project, vw$input$project2, vw$input$dir_out_tree, vw$input$clobber)
 if (is.null(dir_fshome) || dir_fshome == "") stop("dir_fshome is not specified. Please set the global variable FREESURFER_HOME.")

 final_path <- paste0(dir_out2, "/")
 if (vw$input$file_out_tree) final_path <- paste0(final_path, vw$input$project2)
 dat <- if (vw$input$file_out_tree) ".fwhm.dat" else "fwhm.dat"
 finalmask <- if (vw$input$file_out_tree) ".finalMask.mgh" else "finalMask.mgh"
 dir_tmp2 <- file.path(dir_tmp, vw$input$project2)
 backing_mgh <- paste0(dir_tmp2, "_mgh_backend")

 paths <- list()
 paths[["dir_tmp"]] <- dir_tmp
 paths[["dir_tmp2"]] <- dir_tmp2
 paths[["dir_subj"]] <- dir_subj
 paths[["dir_fshome"]] <- dir_fshome
 if(!identical(dir_out, dir_out2)) paths[["orig_dir_out"]] <- dir_out
 paths[["dir_out"]] <- dir_out2
 paths[["backing_mgh"]] <- backing_mgh
 paths[["mask_path"]] <- mask_path
 paths[["final_path"]] <- final_path
 paths[["est_fwhm_path"]] <- paste0(final_path, dat)
 paths[["final_mask_path"]] <- paste0(final_path, finalmask)
 paths[] <- lapply(paths, sub, pattern = "//", replacement = "/", fixed = TRUE)
 paths[] <- lapply(paths, sub, pattern = "/$", replacement = "")
 paths
}

check_prep_fun <- function(prep_fun){
  tryCatch(get2(prep_fun), error = function(e) stop("Provided `prep_fun` cannot be found."))
  invisible(NULL)
}

check_vertex <- function(vertex, id, md){
  stxt <- c("Please provide a unique, unused name to `vertex`. \n",
            "We only use it as an identifier in the model output.")
  if (vertex == id)
    stop("Provided vertex name is identical to the provided id column name. \n",
         stxt)
  q <- sapply(md, function(x) !id %in% names(x))
  if (any(q))
    stop("Provided vertex name is already present in the dataset. \n",
         stxt)
  NULL
}

qdecr_check_mask <- function(mask, mask_path) {
  if (is.null(mask)) {
    if (!file.exists(mask_path)) stop("Mask file does not exist. \n", mask_path)
    return(load.mgh(mask_path)$x)
  }
  else if (is.numeric(mask) || is.logical(mask)) {
    return(mask)
  } else {
    stop("No mask can be found.")
  }
}
