#' Vertex-wise linear regression (based on RcppEigen)
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. See `?lm`.
#' @param data a required argument that contains a data frame, a list of data frames or an imputed object that is supported by the `imp2list` function (mice, mi, etc.).
#' @param id the name of the id variable that matches the dataset to the Freesurfer output
#' @param hemi "lh" or "rh"
#' @param dir_out the directory where to save the data to (defaults to the current directory)
#' @param project the base name you want to assign to the output files
#' @param n_cores the number of cores to be used
#' @param target the target template (usually "fsaverage")
#' @param mcz_thr the Monte Carlo simulation threshold times 10 (13 = 0.05, 20 = 0.01, 23 = 0.005 30 = 0.001, etc..)
#' @param mgh NOT IMPLEMENTED; path to existing merged mgh file, default is NULL
#' @param mask mgh file to mask analysis; default is to use the cortex label from the target
#' @param mask_path path to the mask; default is the cortex mask that is provided with the QDECR package
#' @param dir_subj directory contain the surface-based maps (mgh files); defaults to SUBJECTS_DIR
#' @param dir_fshome Freesurfer directory; defaults to FREESURFER_HOME
#' @param dir_tmp directory to store the temporary big matrices; useful for shared memory; defaults to `dir_out`
#' @param dir_target directory in which `target` is located (by default `dir_subj`)
#' @param dir_out_tree if TRUE, creates a dir_out/project directory. If FALSE, all output is placed directory into dir_out
#' @param clean_up_bm if TRUE, cleans all big matrices (.bk) that were generated in dir_tmp
#' @param clean_up NOT IMPLEMENTED; will be used for setting cleaning of other files
#' @param clobber if TRUE, ignores already existing directories and writes over them; if FALSE, stops and warns user that a given directory already exists
#' @param verbose if TRUE, writes out standard log; if FALSE, no output is generated
#' @param save if TRUE, saves the output to a .fst file
#' @param save_data if TRUE, includes the raw data + design matrices in the .fst file
#' @param debug NOT IMPLEMENTED; will output the maximal log to allow for easy debugging
#' @return returns an object of classes "vw_fastlm" and "vw".
#' @export

qdecr_fastlm <- function(formula,
                         data,
                         id,
                         hemi,
                         dir_out = getwd(),
                         project,
                         n_cores = 1,
                         target = "fsaverage",
                         mcz_thr = 30,
                         mgh = NULL,
                         mask = NULL,
                         mask_path = file.path(path.package("QDECR"), "extdata", paste0(hemi, ".fsaverage.cortex.mask.mgh")),
                         dir_subj = Sys.getenv("SUBJECTS_DIR"),
                         dir_fshome = Sys.getenv("FREESURFER_HOME"),
                         dir_tmp = dir_out,
                         dir_target = dir_subj,
                         dir_out_tree = TRUE,
                         clean_up_bm = TRUE,
                         clean_up = TRUE,
                         clobber = FALSE,
                         verbose = TRUE,
                         save = TRUE,
                         save_data = TRUE,
                         debug = FALSE){

# Take apart the formula
terms <- attr(terms(formula), "factors")
rt <- rownames(terms)
ct <- colnames(terms)

# Check if there is a vertex-wise measure present
qt <- c("qdecr_thickness", "qdecr_area", "qdecr_area.pial", "qdecr_curv", "qdecr_jacobian_white", "qdecr_pial", "qdecr_pial_lgi", "qdecr_sulc", "qdecr_volume", "qdecr_w-g.pct", "qdecr_white.H", "qdecr_white.K")
if (length(intersect(rt, qt)) == 0) stop("Please specify in the formula one of: ", paste(qt, collapse = ", "))

# Find which one it is
if (sum(qt %in% rt) > 1) stop("qdecr currently cannot handle multiple vertex-wise measures in the formula.")
qid <- rt %in% qt
if (sum(terms[qid, ]) > 1) stop("qdecr currently cannot handle interactions or combinations of the vertex terms.")
qqt <- rt[qid]
qqt2 <- sub("qdecr_", "", qqt)

ff <- as.formula(deparse(formula))

margs <- c(qdecr_decon(RcppEigen::fastLm()), list(formula = ff, data = data, method = 2))

vw <- qdecr(id = id,
            margs = margs,
            hemi = hemi,
            dir_out = dir_out,
            target = target,
            mcz_thr = mcz_thr,
            measure = qqt2,
            mgh = mgh,
            mask = mask,
            mask_path = mask_path,
            model = "RcppEigen::fastLm",
            project = project,
            vertex = qqt,
            dir_tmp = dir_tmp,
            dir_target = dir_target,
            dir_subj = dir_subj,
            dir_fshome = dir_fshome,
            dir_out_tree = dir_out_tree,
            clean_up = clean_up,
            clean_up_bm = clean_up_bm,
            clobber = clobber,
            verbose = verbose,
            debug = debug,
            n_cores = n_cores,
            )

vw$describe$call <- rbind(vw$describe$call, c("call", "qdecr_fastlm call", paste(trimws(deparse(match.call())), collapse = "")))

if (save) qdecr_save(vw, save_data = save_data)
vw
}
