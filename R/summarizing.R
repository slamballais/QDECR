#'Grabs the formula from qdecr_fastlm output
#'
#'Grabs the formula from qdecr_fastlm output
#'
#'See above.
#'
#'@param vw The output object of a qdecr_fastlm call
#'@export

formula.vw_fastlm <- function(vw) vw$model$formula

#'Grabs the number of included participants from qdecr_fastlm output
#'
#'Grabs the number of included participants from qdecr_fastlm output
#'
#'See above.
#'
#'@param vw The output object of a qdecr analysis call
#'@export

nobs.vw_fastlm <- function(vw) vw$results$residuals$nrow

#'Print for vw objects
#'
#'Prints a nicely formatted output for vw objects
#'
#'This is a standard print format for output of qdecr functions
#'
#'@param vw The output object of a qdecr analysis call
#'@export

print.vw <- function(vw) qdecr_print_describe(vw$describe, verbose = TRUE)

#'Grabs the estimated fwhm from qdecr_fastlm output
#'
#'Grabs the estimated fwhm from qdecr_fastlm output
#'
#'See above.
#'
#'@param vw The output object of a qdecr analysis call
#'@export

qdecr_fwhm <- function(vw) return(vw$post$fwhm_est)

#'Shows vw stacks
#'
#'Shows the names of the stacks of vw objects
#'
#'This function shows the stack/contrast names for all the variables
#'within the qdecr analysis. This functions commonly used for looking
#'more specifically at variables.
#'
#'@param vw The output object of a qdecr analysis call
#'@export

stacks <- function(vw) vw$stack$names

#'Summarize qdecr_fastlm clusters
#'
#'Summarizes the significant clusters obtained from qdecr_fastlm
#'
#'This function returns an overview of the significant clusters from qdecr_fastlm,
#'and provides information per cluster. Additionally, it is possible to obtain annotated
#'information, given that it is located in the target template's label directory.
#'
#'@param vw The output object of a qdecr call (e.g. qdecr_fastlm)
#'@param verbose Logical; if TRUE, it outputs some information of the steps; default is FALSE
#'@param annot Logical; if TRUE, `file` will be read in to obtain information on the top regions in that cluster
#'@param file string; the name of the file within the target template's label directory if annot = TRUE
#'@param regions integer; number of regions that should be added if annot = TRUE
#'@return NULL
#'@export
#'
summary.vw_fastlm <- function(vw, verbose = FALSE, annot = FALSE, file = "aparc.annot", regions = 3){
  
  if(!is.numeric(regions)) stop("Provided `annot_length` is not a number.")
  if(!(regions > 0)) stop("Provided `annot_length` is not a positive number.")
  
  # Retrieve the stacks to be loaded
  s <- vw$stack$names
  
  # Load the cluster/ocn data
  cl <- list()
  for (i in seq_along(s)) cl[[i]] <- qdecr_read_ocn(vw, i)
  
  # For each stack, deconstruct which indices each cluster is composed of
  for (i in seq_along(s)) cl[[i]] <- lapply(seq_len(max(cl[[i]]$x)), function(x) which(cl[[i]]$x == x))
  
  # Make an empty data frame
  nr <- sum(sapply(cl, length))
  cm <- paste0("mean_", vw$input$measure)
  nc <- c("variable", "cluster", "n_vertices", cm, "mean_coefficient", "mean_se")
  cs <- lapply(cl, function(x) {
    y <- as.data.frame(matrix(NA, nrow = length(x), ncol = length(nc)))
    names(y) <- nc
    y})
  
  # Get info per stack
  cl2 <- list()
  ct_m <- rep(FALSE, length(vw$post$final_mask))
  for (i in seq_along(s)){
    message2(paste0("Summarizing data for `", s[i], "`"), verbose = verbose)
    
    # Load coef and se
    coo <- qdecr_read_coef(vw, i)$x
    cse <- qdecr_read_se(vw, i)$x
    
    for (j in seq_along(cl[[i]])){
      ct <- cl[[i]][[j]]
      ct_mask <- ct_m
      ct_mask[ct] <- TRUE
      cs[[i]][j, "variable"] <- s[i]
      cs[[i]][j, "cluster"] <- j
      cs[[i]][j, "n_vertices"] <- length(ct)
      cs[[i]][j, cm] <- mean(vw$post$mgh_description$vertex_mean[ct_mask], na.rm = TRUE)
      cs[[i]][j, "mean_coefficient"] <- mean(coo[ct])
      cs[[i]][j, "mean_se"] <- mean(cse[ct])
    }
  }
  cs2 <- do.call("rbind", cs)
  
  # Get cluster annotation information
  if (annot) {
    ca <- qdecr_clusters(vw, name = file)
    if (is.null(ca)) return(cs2)
    ca2 <- lapply(ca, function(x) {
      nr2 <- nrow(x)
      if (nr2 > regions) nr2 <- regions
      x[seq_len(nr2), ]
    })
    ca3 <- lapply(ca2, function(y) apply(y, 1, function(x) paste0(x[1], " (", trimws(x[3]), "%, ", trimws(x[4]), "%)")))
    ca4 <- lapply(ca3, function(x) {
      lt <- length(x)
      if (lt < regions) x[(lt+1):regions] <- NA
      x
    })
    ca5 <- do.call("rbind", ca4)
    colnames(ca5) <- paste0("top_region", seq_len(regions))
    cs2 <- cbind(cs2, ca5)
  }
  cs2
}

qdecr_clusters <- function(vw, name = "aparc.annot") {
  file <- paste0(vw$input$hemi, ".", name)
  path <- file.path(vw$paths$dir_subj, vw$input$target, "label", file)
  if (!file.exists(path)) stop("Provided annotation file does not exist.")
  annot <- load.annot(path)
  ocn_annot <- lapply(seq_along(vw$stack$names), function(x) {
    ocn <- qdecr_read_ocn(vw, x)
    if (any(ocn$x > 0)) annots <- lapply(seq_len(max(ocn$x)), function(i) annot$vd_label[ocn$x == i])
  })
  ocn_annot2 <- do.call(c, ocn_annot[!sapply(ocn_annot, is.null)])
  if(is.null(ocn_annot2)) return(NULL)
  ocn_annot3 <- lapply(ocn_annot2, table)
  ta <- table(annot$vd_label)
  nta <- names(ta)
  ocn_annot4 <- lapply(ocn_annot3, function(x) {
    nx <- names(x)
    new <- nta[!nta %in% nx]
    lnew <- length(new)
    if (lnew > 0) {
      x <- setNames(c(x, rep(0, lnew)), c(nx, new))
      x <- x[match(nta, names(x))]
    }
    x
  })
  ocn_annot5 <- lapply(ocn_annot4, function(x) {
    names(x) <- annot$LUT$LUT_labelname[match(names(x), annot$LUT$LUT_value)]
    na_nx <- is.na(names(x))
    if (sum(na_nx) == 1) names(x)[na_nx] <- "missing label"
    x
  })
  ocn_p1 <- lapply(ocn_annot5, function(x) round(prop.table(x) * 100, 2))
  ocn_p2 <- lapply(ocn_annot5, function(x) round(x / ta * 100, 2))
  ocn_p <- lapply(seq_along(ocn_p1), function(i) {
    x <- data.frame(area = names(ocn_p1[[i]]),
                    size = ocn_annot5[[i]],
                    to_cluster = ocn_p1[[i]],
                    to_area = as.vector(ocn_p2[[i]])
    )
    x <- x[x$to_cluster != 0, ]
    x[order(x$to_cluster, decreasing = TRUE), ]
  })
  ocn_p
}