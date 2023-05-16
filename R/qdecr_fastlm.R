#' Vertex-wise linear regression
#'
#' @inheritParams qdecr
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. See `?lm`.
#' @param weights Optinal vector of weights for the observations. See `?lm`. 
#' @param custom_measure a string that starts with "qdecr_" followed by the name of a surface file that is not created by FreeSurfer by default (e.g. "qdecr_cc" or "qdecr_test"). Note that the surface files MUST be located in the surf subdirectory of each individual's FreeSurfer output, and the files must follow the naming conventions of the other .mgh files (e.g. "lh.test.fwhm10.fsaverage.mgh")
#' @param save if TRUE, saves the output to a .rds file
#' @param save_data if TRUE, includes the raw data + design matrices in the .rds file
#' @return returns an object of classes "vw_fastlm" and "vw".
#' @export

qdecr_fastlm <- function(formula,
                         data,
                         id,
                         hemi,
                         weights = NULL,
                         dir_out = getwd(),
                         project,
                         n_cores = 1,
                         target = "fsaverage",
                         fwhm = ifelse(measure == "pial_lgi", 5, 10),
                         mcz_thr = 0.001,
                         cwp_thr = 0.025,
                         mgh = NULL,
                         mask = NULL,
                         mask_path = system.file("extdata", paste0(hemi, ".fsaverage.cortex.mask.mgh"), package = "QDECR"),
                         dir_subj = Sys.getenv("SUBJECTS_DIR"),
                         dir_fshome = Sys.getenv("FREESURFER_HOME"),
                         dir_tmp = dir_out,
                         dir_out_tree = TRUE,
                         file_out_tree = !dir_out_tree,
                         clean_up_bm = TRUE,
                         clean_up = TRUE,
                         clobber = FALSE,
                         verbose = TRUE,
                         save = TRUE,
                         save_data = TRUE,
                         debug = FALSE,
                         custom_measure = NULL,
                         prep_fun = "prep_fastlm",
                         analysis_fun = "analysis_chunkedlm",
                         chunk_size = 1000){

# Take apart the formula
if (!inherits(formula, "formula")) stop("The supplied formula is not of class `formula`.")
terms <- attr(terms(formula), "factors")
rt <- rownames(terms)
ct <- colnames(terms)

# list of all the vertex-wise measures
qt <- c("qdecr_thickness", "qdecr_area", "qdecr_area.pial", "qdecr_curv", "qdecr_jacobian_white", "qdecr_pial", "qdecr_pial_lgi", "qdecr_sulc", "qdecr_volume", "qdecr_w_g.pct", "qdecr_white.H", "qdecr_white.K")
measure_choices <- sub("qdecr_", "", qt)


# process any custom vertex-wise measures
if (!is.null(custom_measure)) {
  if (!is.character(custom_measure)) stop ("The input to the `custom_measure` argument does not seem to be a character vector.")
  starts_with_qdecr <- grepl("^qdecr_.*", custom_measure)
  if (any(!starts_with_qdecr)) stop("The provided `custom_measure` does not start with 'qdecr_'.")
  if (any(custom_measure %in% qt)) stop("The provided `custom_measure` is one of the default maps.")
  qt <- c(qt, custom_measure)
  measure_choices <- c(measure_choices, custom_measure)
}

# Check if there is a vertex-wise measure present
n_vertex_vars <- length(intersect(rt, qt))
if (n_vertex_vars == 0) stop("The formula does not contain one of the default FreeSurfer surface measures (e.g. qdecr_thickness). If the user specified a custom surface map, use the `custom_measure` argument.")

# Find which one it is
if (sum(qt %in% rt) > 1) stop("qdecr currently cannot handle multiple vertex-wise measures in the formula.")
qid <- rt %in% qt
if (sum(terms[qid, ]) > 1) stop("qdecr currently cannot handle interactions or combinations of the vertex terms.")
qqt <- rt[qid]
measure <- sub("qdecr_", "", qqt)
ff <- stats::as.formula(paste(deparse(formula), collapse = " "))
margs <- c(qdecr_decon(RcppEigen::fastLm()), list(formula = ff, data = data, method = 2))
if (!is.null(weights)) margs <- c(margs, list(weights = weights))
vw <- qdecr(id = id,
            margs = margs,
            hemi = hemi,
            dir_out = dir_out,
            target = target,
            fwhm = fwhm,
            mcz_thr = mcz_thr,
            cwp_thr = cwp_thr,
            measure = measure,
            measure_choices = measure_choices,
            mgh = mgh,
            mask = mask,
            mask_path = mask_path,
            model = "RcppEigen::fastLm",
            project = project,
            vertex = qqt,
            dir_tmp = dir_tmp,
            dir_subj = dir_subj,
            dir_fshome = dir_fshome,
            dir_out_tree = dir_out_tree,
            file_out_tree = file_out_tree,
            clean_up = clean_up,
            clean_up_bm = clean_up_bm,
            clobber = clobber,
            verbose = verbose,
            debug = debug,
            n_cores = n_cores,
            prep_fun = prep_fun,
            analysis_fun = analysis_fun,
            chunk_size = chunk_size
            )

vw$describe$call <- rbind(vw$describe$call, c("call", "qdecr_fastlm call", paste(trimws(deparse(match.call())), collapse = "")))
stacks_df <- data.frame(stack_number = seq_along(stacks(vw)), stack_name = stacks(vw))
utils::write.table(stacks_df, file.path(vw$paths$dir_out, "stack_names.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
summary_df <- summary(vw, annot = TRUE)
utils::write.table(summary_df, file.path(vw$paths$dir_out, "significant_clusters.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
if (save) qdecr_save(vw, save_data = save_data)
vw
}
