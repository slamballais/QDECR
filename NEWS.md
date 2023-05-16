# QDECR 0.10.0: Stroopwafel

Version 0.10.0 introduces a number of quality-of-life changes to make life a bit sweeter.

## New features


## Minor Tweaks


## Bug fixes
* [9c101e5](https://github.com/slamballais/QDECR/commit/9c101e5fec5f7a7a4ef5d65ac243a3dc71d129cc): Added a check to make sure that the supplied `formula` is of class `formula`.

# QDECR 0.9.0: Lausanne

Version 0.9.0 is the first update after publication of the QDECR manuscript in Frontiers in Neuroinformatics. The latter is based in Lausanne, Switzerland. This version introduces weighted regression and a number of minor tweaks. 

## New features
* [#47](https://github.com/slamballais/QDECR/pull/47): `qdecr_fastlm` can now take weights via the `weights` argument, which act as observation weights in the linear regression. This argument works similar to `lm` (addresses [#36](https://github.com/slamballais/QDECR/issues/36)).

## Bug fixes
* [#43](https://github.com/slamballais/QDECR/pull/43): The specified number of cores is checked against `bigstatsr::nb_cores` (inside `QDECR:::check_cores`). However, it seems that it returns 0 cores when only 1 core is found (given that it always omits 1 core). We thus rewrote `check_cores` to use `parallel::detectCores`, and to make sure that it cannot return 0 cores (fixes [#42](https://github.com/slamballais/QDECR/issues/42)).
* [#43](https://github.com/slamballais/QDECR/pull/43): Fixed the text in the error message of `QDECR:::check_cores`.
* [#45](https://github.com/slamballais/QDECR/pull/45): Unsmoothed q-cached surface files can now be analyzed. Normally, QDECR would look for files containing `fwhmX`, where X is the FWHM in mm. However, unsmoothed files do not contain this part. The code was rewritten to check if fwhm == 0, and in those cases it will not insert that into the file names. This was fixed by setting `fwhmc` to "" in `QDECR:::qdecr_check`, and by tweaking `QDECR:::qdecr_prep_mgh` (fixes [#41](https://github.com/slamballais/QDECR/issues/41)).
* [#51](https://github.com/slamballais/QDECR/pull/51): Added an extra check to avoid multiple levels of parallelism, e.g. `n_cores` > 1 while simultaneously having a parallel BLAS library set up.
* [40f04d4](https://github.com/slamballais/QDECR/commit/40f04d44c08991b752cf4bf4fd0a9f7c525ca409): Added some functions to the NAMESPACE (stats::na.fail, stats::model.weights).
* [8f19d58](https://github.com/slamballais/QDECR/commit/8f19d585338a1bb0afc575c0e164ff6b0d90b5b5): Since 4.0.0 it is possible to get a warning message from `readChar`. The `load.annot` function was modified by replacing `readChar` with the equivalent call of `readBin`.

## Minor tweaks
* [#46](https://github.com/slamballais/QDECR/pull/46): Tweaked the vertex-wise analysis code (`QDECR:::analysis_chunkedlm`) to be faster and be more memory efficient (`se` is now calculated without the intermediate storage of `s2`) (addresses [#38](https://github.com/slamballais/QDECR/issues/38)).
* [306826c](https://github.com/slamballais/QDECR/commit/306826c0691fe6d94dd22b1ed6dfa0af96cc9aa0): Added an extra check in `QDECR:::qdecr_prep_mgh` so that it also checks whether the .mgh files actually exist. It will output which subjects are missing the surface files (if N < 20) or just say that people are missing surface files (if N >= 20). This adds a little bit of runtime, but it pays off (fixes [#34](https://github.com/slamballais/QDECR/issues/34)).
* [26ac549](https://github.com/slamballais/QDECR/commit/26ac549bf05fef9302f20cfb1d177fdfff277954): When there is missing data in the design matrix, `QDECR:::qdecr_fastlm` will throw an error (i.e., `na.action = na.fail`). This is to avoid problems downstream. We opted for `na.fail` and not `na.omit`-like behavior, because users who are unaware of missingness would then only find out after running QDECR that they had missingness.

# QDECR 0.8.5

## Bug fixes 
* [cfc9a35](https://github.com/slamballais/QDECR/commit/cfc9a35edbb52586bbd7777865577d49b73368fb): Fixed a line in `check_id` to explicitly state `drop = TRUE` in `md[, id, drop = TRUE]`. 
* [39b6173](https://github.com/slamballais/QDECR/commit/39b617309913595dc8883bf849085db27fb29912): Fixed that a call to `model.matrix` would supply `stats::contrasts` in `prep_fastlm`; this led to harmless warnings.

## Minor tweaks
* [5043753](https://github.com/slamballais/QDECR/commit/504375366a795926521f6ca640ee6cabc931b757): To let users analyze custom surface maps, we added the `custom_measure` argument to `qdecr_fastlm`. This argument lets users specify any name for a surface map that they want, provided that [1] it starts with "qdecr_" (e.g. "qdecr_test"), [2] such a file is located in the subj subdirectory of the FreeSurfer output directories, [3] the surface files follow the same naming convention as the other surface maps that FreeSurfer outputs (e.g. "lh.test.fwhm10.fsaverage.mgh").

# QDECR 0.8.4

## Bug fixes
* [20d3be2](https://github.com/slamballais/QDECR/commit/20d3be24f92121bfee8f1c114f86494dea51ae65): Fixed that the model output files (p.mgh, t.mgh, etc) also contained the mcz threshold (e.g. cache.th30).
* [18c9de8](https://github.com/slamballais/QDECR/commit/18c9de8e47bc9af151a5fc3db9ab2fa384358d63): Lots of comments removed/modified, and code was cleaned up a bit.
* [015d8a1](https://github.com/slamballais/QDECR/commit/015d8a19520aa00a996732b8027d01221f6dc076)/[e0aa3fe](https://github.com/slamballais/QDECR/commit/e0aa3fe9cf5c7ccb7945ff615d5891e77e74b58a): Wrote in a preventative stop for when `dir_out_tree = FALSE` and `clobber = TRUE` are combined. This can easily delete important directories unintentionally, as has happened to the authors.

## Minor Tweaks
* [e7adb17](https://github.com/slamballais/QDECR/commit/e7adb175575877e573723493d3e78e6aecfc8ea2): Added the `file_out_tree` argument, which controls whether output files also contain the full project name. By default, it is the inverse of `dir_out_tree`.
* [077e4fd](https://github.com/slamballais/QDECR/commit/077e4fd92bd120bad14c9a75c00f0e8a6f55d63f): Fixed a lot of the documentation.

# QDECR 0.8.3

## Minor tweaks
* [f019f81](https://github.com/slamballais/QDECR/commit/f019f819c01405bfef0f377ec50fcab03ade3718): The `mcz_thr` argument (for `qdecr` and `qdecr_fastlm`) now accepts: 13/1.3/0.05, 20/2.0/0.01, 23/2.3/0.005, 30/3.0/0.001, 33/3.3/0.0005, 40/4.0/0.0001.
* [f019f81](https://github.com/slamballais/QDECR/commit/f019f819c01405bfef0f377ec50fcab03ade3718): A new function `qdecr_mcz_thr` makes sure that the value is converted to 13/20/23/30/33/40. 

# QDECR 0.8.2

## Bug fixes
* [6638efb](https://github.com/slamballais/QDECR/commit/6638efb2a45492c13e284644b847591b3c9727be): The `fst` package was noted as imported package, but we never implemented functionality from it. Thus, all reference to it was removed.

## New (minor) tweaks
* [3d479f4](https://github.com/slamballais/QDECR/commit/3d479f4f47d5ad851ec2de3b8909bd3cd8faa603): The version is now properly displayed when running `qdecr`. This version is updated dynamically using `packageVersion("QDECR")`. The website (www.qdecr.com) was also added.
* [95a4d90](https://github.com/slamballais/QDECR/commit/95a4d90c6be352be8b7daf734aee126c966b971b): Exported `load.annot` and the `qdecr_read` functions.
* [fa3de94](https://github.com/slamballais/QDECR/commit/fa3de94b850c3e95f4ca6374f4f9ec3b3486788d): Renamed the internal `runMriSurfCluster` function to `run_mri_surf_cluster`.

# QDECR 0.8.1

## Bug fixes
* [e86d2e1](https://github.com/slamballais/QDECR/commit/e86d2e116f9b8976f7004044ed1b3dae7a0df629): If the estimated smoothness is below 1, we now increase it to 1 to avoid problems down the line.
* [04fcd1b](https://github.com/slamballais/QDECR/commit/04fcd1ba97a770e087c20ff902485122fb292683): Fixed a bug in `qdecr_clusters` where it assumes that there is always at least 1 cluster significant.
* [251e668](https://github.com/slamballais/QDECR/commit/251e668278fa945e9662675c279dc2d07877ac25): "w-g.pct" files can now be used as a measure by specifying "qdecr_w_g.pct" (underscore instead of hyphen) (fixes [#19](https://github.com/slamballais/QDECR/issues/19)).
* [dc81b89](https://github.com/slamballais/QDECR/commit/dc81b89ea6bece71831f6abc0e4eebcebc26f51e): Added `fwhm` argument to `qdecr_fastlm`, which was missing before (fixes [#18](https://github.com/slamballais/QDECR/issues/18)).

## New (minor) features
* [4f024f3](https://github.com/slamballais/QDECR/commit/4f024f38e5bf5a6277c847b9bf2371ffa58521b1): Added `cwp_thr` argument to `qdecr_fastlm` and `qdecr` to set the further cluster-wise p-value adjustment (default is 0.025 due to having 2 hemispheres, thus 0.05 / 2) (fixes [#23](https://github.com/slamballais/QDECR/issues/23)).
* [08ae63b](https://github.com/slamballais/QDECR/commit/08ae63b23033d0571dfe6509206d95e09863316c): Automatically output two extra files: "significant_clusters.txt" (the output of `summary(vw, annot = TRUE)`) and "stack_names.txt" (the output of `stacks(vw)` and the corresponding stack numbers) (fixes [#17](https://github.com/slamballais/QDECR/issues/17)).
* [11e2f55](https://github.com/slamballais/QDECR/commit/11e2f55e29a1e0aac0eb39ff5ac3fb92dbbe7f95): Modified `freeview` and `qdecr_snap`. The `mask` argument is now called `sig`. Furthermore, the ranges for the overlay colors are determined dynamically. Finally, users can now set any arguments to Freeview for manipulating surface files (see `freeview --help` on the command line) (fixes [#20](https://github.com/slamballais/QDECR/issues/20)).

# QDECR 0.8.0: Momo

Version 0.8.0 is the first update after public release. It fixes a bunch of mistakes, introduces further modularization, improves the speed and also reduces the RAM load.

## New functions and features

* [#2](https://github.com/slamballais/QDECR/pull/2): Added the input argument `dir_target`, so that the target can be specified flexibly (fixes [#1](https://github.com/slamballais/QDECR/issues/1)). 
* [#4](https://github.com/slamballais/QDECR/pull/4)/[#6](https://github.com/slamballais/QDECR/pull/6): Modularized `qdecr_model` and added the input argument `prep_fun`, so that users can choose and create their own prep functions (fixes [#3](https://github.com/slamballais/QDECR/issues/3)).
* [#10](https://github.com/slamballais/QDECR/pull/10): Modularized `qdecr_analysis` and added the input argument `analysis_fun`, so that users can choose and create their own analysis functions.
* [#10](https://github.com/slamballais/QDECR/pull/10): The internal function to run vertex-wise analyses, called `vertexwise`, now processes regressions in chunks, i.e. more than 1 vertex at a time. This has led to a considerable upgrade in speed, especially for smaller datasets. Chunk size can be controlled with `chunk_size` (fixes [#7](https://github.com/slamballais/QDECR/issues/7)).

## Bug fixes

* [#10](https://github.com/slamballais/QDECR/pull/10): We removed lots of unnecessary dependencies.
* [#10](https://github.com/slamballais/QDECR/pull/10): We fixed the referencing to the default mask.

# QDECR 0.7.0: OHBM

Version 0.7.0 is the first version that is publically released. It is also the version that was presented at OHBM 2019 and was thus named "OHBM". As we do not have any formal news for versions before 0.7.0, we will keep it brief.

## New functions and features

* Creation of the QDECR package (before 0.7.0 it was a series of associated scripts)
* Creation of the `vw` object and associated functions
* Handling of imputed datasets (via `imp2list`)
* Creation of all summary and plot functions
* Creation of the `stack` concept
