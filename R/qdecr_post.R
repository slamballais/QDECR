


calc_fwhm <- function(final_path, final_mask_path, est_fwhm_path, hemi, eres, mask = NULL, path_target, verbose = FALSE) {
  cmdStr <- paste("mris_fwhm", "--i", eres, "--hemi", hemi, "--subject", basename(path_target), "--sd", dirname(path_target), "--prune", "--cortex", "--dat", est_fwhm_path, "--out-mask", final_mask_path)
  
  if (!is.null(mask)) paste(cmdStr, "--mask", mask)
  message2(verbose = verbose, cmdStr)
  system(cmdStr, ignore.stdout = !verbose)
  fwhm <- read.table(est_fwhm_path)
  fwhm <- round(fwhm)
  return(fwhm)
}

runMriSurfCluster <- function(final_path, dir_fshome, hemi, pval, fwhm, mask_path = NULL, cwpThr = 0.025, mczThr = NULL, csdSign = "abs", verbose = FALSE, stack, stack_name) {

  if (is.null(mczThr)){
    mczThr = paste0("th", as.character(30))
  }
  else {
    mczThr =  paste0("th", as.character(mczThr))
  }

  oBaseName <- paste(final_path, paste0("stack", stack), "cache", mczThr, csdSign, "sig", sep = ".")
  if (fwhm < 10) fwhm <- paste0("0", fwhm)
  csd <- file.path(dir_fshome, "average/mult-comp-cor/fsaverage", hemi, paste0("cortex/fwhm", fwhm), csdSign, mczThr, "mc-z.csd")
  cmdStr <- paste("mri_surfcluster",
                  "--in", pval,
                  "--csd", csd, 
                  "--cwsig", paste0(oBaseName, ".cluster.mgh"),
                  "--vwsig", paste0(oBaseName, ".voxel.mgh"),
                  "--sum", paste0(oBaseName, ".cluster.summary"),
                  "--ocn",  paste0(oBaseName, ".ocn.mgh"),
                  "--oannot", paste0(oBaseName, ".ocn.annot"),
                  "--annot", "aparc",
                  "--csdpdf",  paste0(oBaseName, ".pdf.dat"),
                  "--cwpvalthresh", cwpThr,
                  "--o", paste0(oBaseName, ".masked.mgh"),
                  "--no-fixmni",
                  "--surf", "white")
  cmdStr <- if (!is.null(mask_path)) paste(cmdStr, "--mask", mask_path) else paste(cmdStr, "--cortex")
  system(cmdStr, ignore.stdout = !verbose)
  NULL
}