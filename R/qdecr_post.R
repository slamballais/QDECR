calc_fwhm <- function(final_path, final_mask_path, est_fwhm_path, hemi, eres, mask = NULL, target = "fsaverage", verbose = FALSE) {
  cmdStr <- paste("mris_fwhm", "--i", eres, "--hemi", hemi, "--subject", target, "--prune", "--cortex", "--dat", est_fwhm_path, "--out-mask", final_mask_path)
  if (!is.null(mask)) paste(cmdStr, "--mask", mask)
  message2(verbose = verbose, cmdStr)
  system(cmdStr, ignore.stdout = !verbose)
  fwhm <- read.table(est_fwhm_path)
  fwhm <- round(fwhm)
  return(fwhm)
}

run_mri_surf_cluster <- function(vw, pval, fwhm, mask_path = NULL, cwp_thr = 0.025, mcz_thr = 30, csd_sign = "abs", verbose = FALSE, stack, stack_name) {
  mcz_thr2 <- paste0("th", mcz_thr)
  if (fwhm < 10) fwhm <- paste0("0", fwhm)
  csd <- file.path(vw$paths$dir_fshome, "average/mult-comp-cor/fsaverage", vw$input$hemi, paste0("cortex/fwhm", fwhm), csd_sign, mcz_thr2, "mc-z.csd")
  cmd_str <- paste("mri_surfcluster",
                   "--in", pval,
                   "--csd", csd, 
                   "--cwsig", vw$stack$cluster.mgh[[stack]],
                   "--vwsig", vw$stack$voxel.mgh[[stack]],
                   "--sum", vw$stack$cluster.summary[[stack]],
                   "--ocn",  vw$stack$ocn.mgh[[stack]],
                   "--oannot", vw$stack$ocn.annot[[stack]],
                   "--annot aparc",
                   "--cwpvalthresh", cwp_thr,
                   "--o", vw$stack$masked.mgh[[stack]],
                   "--no-fixmni",
                   "--surf", "white")
  cmd_str <- if (!is.null(mask_path)) paste(cmd_str, "--mask", mask_path) else paste(cmd_str, "--cortex")
  system(cmd_str, ignore.stdout = !verbose)
  NULL
}