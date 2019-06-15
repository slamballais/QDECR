

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
                        project, dir_out_tree, clobber, fwhm, n_cores){

  # verify that id is correct and exists within data
  check_id(id, md)

  # verify that vertex is correct, is not id and is not in data
  check_vertex(vertex, id, md)

  # Verify the integrity of fwhm
  check_fwhm(fwhm)

  # Check whether the number of cores is correct
  check_cores(n_cores)

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
  input[["clobber"]] <- clobber
  input[["fwhm"]] <- fwhm
  input[["fwhmc"]] <- paste0("fwhm", fwhm)
  input[["n_cores"]] <- n_cores
  input
}

check_cores <- function(n_cores){
  rec_cores <- bigstatsr::nb_cores()
  if (n_cores > rec_cores) stop("You specified `", n_cores, "` to be too high. Recommended is ", rec_cores)
  NULL
}

check_dir_out <- function(dir_out, project, project2, dir_out_tree, clobber){
  if (!dir.exists(dir_out)) {
      if(!dir.exists(dirname(dir_out)))
        stop("The provided `dir_out` does not exist, but neither does it parent.")
      dir.create(dir_out)
  } else {
    if(clobber == FALSE && dir_out_tree == FALSE)
      stop("`dir_out` already exists, and both `clobber` and `dir_out_tree` are FALSE. \n",
           "Either set `clobber` = TRUE to overwrite the output directory, \n",
           "or set `dir_out_tree` = TRUE and provide `project` to \n",
           "create a new subdirectory.")
    if(clobber == TRUE && dir_out_tree == FALSE){
      unlink(dir_out, recursive = TRUE)
      dir.create(dir_out)
    }
  }
  if(dir_out_tree){
    if(missing(project))
      stop("`dir_out_tree` is set to TRUE but no `project` provided.")
    p <- file.path(dir_out, project2)
    if(dir.exists(p)) {
      if (clobber == FALSE)
        stop("The output directory ", p, " already exists and `clobber` = FALSE. \n",
             "Set `clobber` = TRUE to overwrite the existing output directory.")
      unlink(p, recursive = TRUE)
    }
    dir.create(p)
    dir_out <- p
  }
  dir_out
}

check_dir_subj <- function(dir_subj, md, id){
  if(missing(dir_subj)) stop("No `dir_subj` argument provided. \n",
                             "Please provide a path to the directory that contains ",
                             "all the vertex input data.")
  if(!dir.exists(dir_subj))
    stop("The provided `dir_subj` does not seem to exist.")
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

check_dir_target <- function(dir_target, target){
  path_target <- file.path(dir_target, target)
  if(!dir.exists(path_target)){
    stop("The provided `target` directory does not seem to exist.")
  }
  path_target
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
  if (length(unique(md[[1]][,id])) == 1)
    stop("The provided `id` column only has 1 unique value in your data, suggesting only 1 person is loaded in.")
  NULL
}

check_paths <- function(vw, dir_tmp, dir_subj, dir_out, dir_target, dir_fshome, mask_path){
 check_dir_subj(dir_subj, vw$input$md, vw$input$id)
 path_target <- check_dir_target(dir_target, vw$input$target)
 dir_out2 <- check_dir_out(dir_out, vw$input$project, vw$input$project2, vw$input$dir_out_tree, vw$input$clobber)
 if (is.null(dir_fshome) || dir_fshome == "") stop("dir_fshome is not specified. Please set the global variable FREESURFER_HOME.")

 final_path <- file.path(dir_out2, vw$input$project2)
 dir_tmp2 <- file.path(dir_tmp, vw$input$project2)
 backing_mgh <- paste0(dir_tmp2, "_mgh_backend")

 paths <- list()
 paths[["dir_tmp"]] <- dir_tmp
 paths[["dir_tmp2"]] <- dir_tmp2
 paths[["dir_subj"]] <- dir_subj
 paths[["dir_fshome"]] <- dir_fshome
 if(!identical(dir_out, dir_out2)) paths[["orig_dir_out"]] <- dir_out
 paths[["dir_out"]] <- dir_out2
 paths[["dir_target"]] <- dir_target
 paths[["path_target"]] <- path_target
 paths[["backing_mgh"]] <- backing_mgh
 paths[["mask_path"]] <- mask_path
 paths[["final_path"]] <- final_path
 paths[["est_fwhm_path"]] <- paste0(final_path, ".fwhm.dat")
 paths[["final_mask_path"]] <- paste0(final_path, ".finalMask.mgh")

 paths[] <- lapply(paths, sub, pattern = "//", replacement = "/", fixed = TRUE)
 paths[] <- lapply(paths, sub, pattern = "/$", replacement = "")
 paths
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
