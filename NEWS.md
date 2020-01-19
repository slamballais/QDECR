# QDECR 0.8.4

## Bug fixes
* Fixed that the model output files (p.mgh, t.mgh, etc) also contained the mcz threshold (e.g. cache.th30).
* Lots of comments removed/modified, and code was cleaned up a bit

# QDECR 0.8.3

## Minor tweaks
* The `mcz_thr` argument (for `qdecr` and `qdecr_fastlm`) now accepts: 13/1.3/0.05, 20/2.0/0.01, 23/2.3/0.005, 30/3.0/0.001, 33/3.3/0.0005, 40/4.0/0.0001.
* A new function `qdecr_mcz_thr` makes sure that the value is converted to 13/20/23/30/33/40. 

# QDECR 0.8.2

## Bug fixes
* The `fst` package was noted as imported package, but we never implemented functionality from it. Thus, all reference to it was removed.

## New (minor) tweaks
* The version is now properly displayed when running `qdecr`. This version is updated dynamically using `packageVersion("QDECR")`. The website (www.qdecr.com) was also added.
* Exported `load.annot` and the `qdecr_read` functions.
* Renamed the internal `runMriSurfCluster` function to `run_mri_surf_cluster`.

# QDECR 0.8.1

## Bug fixes
* If the estimated smoothness is below 1, we now increase it to 1 to avoid problems down the line.
* Fixed a bug in qdecr_clusters where it assumes that there is always at least 1 cluster significant.
* "w-g.pct" files can now be used as a measure by specifying "qdecr_w_g.pct" (underscore instead of hyphen).
* Added `fwhm` argument to `qdecr_fastlm`, which was missing before.

## New (minor) features
* Added `cwp_thr` argument to `qdecr_fastlm` and `qdecr` to set the further cluster-wise p-value adjustment (default is 0.025 due to having 2 hemispheres, thus 0.05 / 2).
* Automatically output two extra files: "significant_clusters.txt" (the output of `summary(vw, annot = TRUE)`) and "stack_names.txt" (the output of `stacks(vw)` and the corresponding stack numbers).
* Modified `freeview` and `qdecr_snap`. The `mask` argument is now called `sig`. Furthermore, the ranges for the overlay colors are determined dynamically. Finally, users can now set any arguments to Freeview for manipulating surface files (see `freeview --help` on the command line).

# QDECR 0.8.0: Momo

Version 0.8.0 is the first update after public release. It fixes a bunch of mistakes, introduces further modularization, improves the speed and also reduces the RAM load.

## New functions and features

* Added the input argument `dir_target`, so that the target can be specified flexibly. 
* Modularized `qdecr_model` and added the input argument `prep_fun`, so that users can choose and create their own prep functions.
* Modularized `qdecr_analysis` and added the input argument `analysis_fun`, so that users can choose and create their own analysis functions.
* The internal function to run vertex-wise analyses, called `vertexwise`, now processes regressions in chunks, i.e. more than 1 vertex at a time. This has led to a considerable upgrade in speed, especially for smaller datasets. Chunk size can be controlled with `chunk_size`

## Bug fixes

* We removed lots of unnecessary dependencies.
* We fixed the referencing to the default mask.

# QDECR 0.7.0: OHBM

Version 0.7.0 is the first version that is publically released. It is also the version that was presented at OHBM 2019 and was thus named "OHBM". As we do not have any formal news for versions before 0.7.0, we will keep it brief.

## New functions and features

* Creation of the QDECR package (before 0.7.0 it was a series of associated scripts)
* Creation of the vw object and associated functions
* Handling of imputed datasets (via imp2list)
* Creation of all summary and plot functions
* Creation of the `stack` concept
