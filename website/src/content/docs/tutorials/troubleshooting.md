---
title: Troubleshooting and FAQ
description: 'Fixes for common problems, answers to frequent questions, and the known limitations.'
order: 9
---

Most problems stop QDECR with a message of its own. Search this page for the message's first words; if it is not here, [ask on GitHub](/help#where-to-ask), with the details [a bug report needs](/help#reporting-a-bug).

## Installation problems

### A package fails to compile

QDECR's dependencies, such as RcppEigen and bigstatsr, compile C++ code while they install. If the installation stops with compiler errors, R has no compiler: install `r-base-dev` on Debian or Ubuntu, run `xcode-select --install` on macOS, or load a compiler module on a cluster. Then install QDECR again.

### R can't find FreeSurfer

QDECR stops with

```text
dir_fshome is not specified. Please set the global variable FREESURFER_HOME.
```

or gets as far as the smoothness estimate and fails there, with the shell saying `mris_fwhm: not found` and R that it cannot open `fwhm.dat`:

```text
sh: 1: mris_fwhm: not found
cannot open file 'results/lh.age_sex.thickness/fwhm.dat': No such file or directory
```

QDECR runs FreeSurfer's programs by name, so R must start with FreeSurfer set up: `FREESURFER_HOME` set, and FreeSurfer's `bin` directory on the `PATH`. R started from a terminal where `SetUpFreeSurfer.sh` has run has both. R started any other way, as RStudio usually is, does not; give it the variables in the file `~/.Renviron`, with your own paths:

```text
FREESURFER_HOME=/usr/local/freesurfer
SUBJECTS_DIR=/data/subjects
PATH=${FREESURFER_HOME}/bin:${PATH}
```

Restart R afterwards, and run the checks under [Check that it works](/get-started#check-that-it-works).

### FreeSurfer can't find its licence

If `mris_fwhm` or `mri_surfcluster` complains about a missing licence, FreeSurfer's free licence file is not where it looks: `$FREESURFER_HOME/license.txt`. Put it there, or name its location with `FS_LICENSE`, in the shell or in `~/.Renviron`.

### magick won't install

The `magick` package needs ImageMagick's development files: `libmagick++-dev` on Debian and Ubuntu. Only `qdecr_snap()` needs `magick`; everything else works without it.

## Errors while running

### The `id` values were not found in `dir_subj`

```text
The following `id` values were not found in `dir_subj`: 1, 2, 3.
```

The values in the ID column must be the names of the subjects' directories in the subjects directory, exactly. Common causes: IDs stored as numbers (`1`) where the directories are named `sub-001`, a different `SUBJECTS_DIR` in R than in the shell, or subjects not yet processed. With more than ten missing, QDECR only gives the count.

### The subjects do not have the surface file

```text
The following subjects do not have lh.thickness.fwhm10.fsaverage.mgh: sub-001, sub-002
```

The subjects exist but lack the file for this measure, hemisphere and smoothing. Usually `recon-all -qcache` has not run for them, or the `fwhm` you asked for is not one `-qcache` made. Look in the subject's `surf` directory for what is there.

### The target directory does not exist

```text
The provided `target` directory does not seem to exist.
```

`fsaverage` is not in the subjects directory. [Link it in](/get-started#check-that-it-works) from FreeSurfer.

### Missing values in object

```text
Error in na.fail.default(...) : missing values in object
```

A variable in the formula has missing values, and QDECR does not drop subjects silently. Remove the incomplete rows, or [impute](/tutorials/imputed-data) the missing values.

### The output directory already exists

```text
The output directory results/lh.age_sex.thickness already exists and `clobber` = FALSE.
```

An analysis with this project name has run before. Give the new one another `project` name, or set `clobber = TRUE` to replace the old results, which deletes that directory first.

If the message names a `.bk` file instead, an earlier run was cut off and left its temporary files behind in `dir_tmp`. Delete them, or run with `clobber = TRUE`.

### Too many cores, or a parallel BLAS

```text
You specified `n_cores` to be too high (16). Recommended is 7
`n_cores` > 1, but there already seems to be a parallel BLAS library present.
```

QDECR accepts at most one core fewer than the machine has, and does not combine several processes with a multithreaded BLAS. [Performance and memory](/tutorials/performance#blas) explains the choice between the two.

### The design matrix is not full rank

```text
The design matrix is NOT full rank. Please check if you have collinear columns in your data.
```

One column of the design is a combination of others: two variables carry the same information, a variable is constant, or a factor has a level that no subject has. `model.matrix()` on your formula and data shows the columns.

### The formula has no vertex measure

```text
The formula does not contain one of the default FreeSurfer surface measures (e.g. qdecr_thickness).
```

The left-hand side is not one of the [measures QDECR knows](/tutorials/formulas-and-design#the-vertex-measure): a typo, a transformation such as `log(qdecr_thickness)`, or a map of your own that needs `custom_measure`.

### The threshold is not accepted

```text
`mcz_thr` does not have an accepted value (13/1.3/0.05, 20/2.0/0.01, 23/2.3/0.005, 30/3.0/0.001, 33/3.3/0.0005, 40/4.0/0.0001).
```

FreeSurfer's simulations exist for these six cluster-forming thresholds only. [Choosing thresholds](/tutorials/statistics#choosing-thresholds) lists them.

### Nothing to plot

```text
Stack does not contain information (e.g. because of no significant findings), aborting plot.
```

`freeview()` and `qdecr_snap()` show the significant clusters of a stack, and this one has none. `summary(out)` lists the stacks that do. To see a map without the clusters, pass `sig = FALSE`.

### R runs out of memory

R stops with `cannot allocate vector`, or is killed without a message. Each process holds a chunk of the data for every imputed dataset, and with `dir_tmp = "/dev/shm"` the temporary files take memory too. Lower `chunk_size` or `n_cores`, or move `dir_tmp` to disk; [Performance and memory](/tutorials/performance) shows how much each needs.

### The smoothness was reduced

```text
Estimated smoothness is 34, which is really high. Reduced to 30.
```

Not an error. FreeSurfer's simulations go up to 30 mm, so QDECR uses the smoothest one there is. Smoothness this high often means the data were smoothed heavily beforehand, or that something is systematically wrong with some subjects' data; look at `hist(out, qtype = "subject")` for outliers.

## Frequently asked questions

### Do I need to know R?

Yes. QDECR is used from R, and works like R's other model functions: if you have fitted a model with `lm()`, you know most of what it needs. You also need a working FreeSurfer, but no FreeSurfer commands beyond `recon-all`.

### How does QDECR compare with QDEC and mri_glmfit?

QDEC, FreeSurfer's graphical tool for group analysis, and `mri_glmfit` fit the same kind of model at every vertex, and `mri_glmfit-sim` applies the same cluster-wise correction QDECR does. QDECR differs in how you get there: an R formula and a data frame instead of an FSGD file and contrast matrices, imputed datasets, weights, and samples of thousands. Each coefficient of the model is tested on its own, rather than through a contrast you write.

### How does it compare with other tools?

- [PALM](https://github.com/andersonwinkler/PALM) tests surface data with permutations, which avoid the assumptions of FreeSurfer's simulations.
- [BrainStat](https://brainstat.readthedocs.io/), in Python and MATLAB, and its predecessor SurfStat fit models on the surface and correct with random field theory.
- [fsbrain](https://github.com/dfsp-spirit/fsbrain) and [freesurferformats](https://github.com/dfsp-spirit/freesurferformats) are R packages that read, write and plot FreeSurfer data, rather than analyse it.
- [VertexWiseR](https://github.com/CogBrainHealthLab/VertexWiseR) is another R package for vertex-wise analysis, of the whole cortex and of the hippocampus.
- [verywise](https://github.com/SereDef/verywise) fits vertex-wise linear mixed models, for longitudinal and multi-site data.

### Can I analyse longitudinal data?

Not with QDECR, which fits one observation per subject. A mixed-model prototype was shown at OHBM 2020 but never released; [verywise](https://github.com/SereDef/verywise) fits such models today.

### Can I analyse both hemispheres at once?

No: an analysis covers one hemisphere. Run it twice, with `hemi = "lh"` and `hemi = "rh"`. The default `cwp_thr` of 0.025 already accounts for the two.

### Does QDECR work with data from other software?

Not directly: it reads FreeSurfer's files from a FreeSurfer subjects directory. A map from other software can be analysed as a [custom measure](/tutorials/formulas-and-design#your-own-maps) once it is resampled to `fsaverage` and saved as one MGH file per subject, named the way FreeSurfer names its maps. Volume data, such as subcortical structures, are outside what QDECR does.

### How long does an analysis take?

It depends on the number of subjects, the number of imputed datasets, the cores, the BLAS and where `dir_tmp` lives. A few hundred subjects with the defaults take minutes; thousands of subjects with dozens of imputed datasets can take hours. [Performance and memory](/tutorials/performance) covers each factor.

### Why is the effect of sex called `sexmale`?

Each column of the design matrix is a [stack](/glossary#stack), named by R: the factor's name followed by the level compared with the reference. [Formulas and design](/tutorials/formulas-and-design#covariates-and-factors) explains how to choose the reference.

### Can I change the code?

Yes. QDECR is free software under the GPL-3 licence, and contributions are welcome: see [Contributing](/help#contributing).

### How do I cite QDECR?

Cite the paper; the [Cite](/cite) page has the reference, its BibTeX, and the FreeSurfer papers to cite alongside it.

## Known limitations

What QDECR 0.9.0 cannot do, or does not do right, and what to do instead.

- **Platforms.** Linux and macOS only; on Windows, [WSL2](/get-started#windows-through-wsl2). macOS itself is untested.
- **One hemisphere, one measure.** An analysis covers one hemisphere, and its formula has one vertex measure, untransformed, on the left-hand side.
- **Linear regression only.** `qdecr_fastlm()` is the one model. For mixed models, see [verywise](https://github.com/SereDef/verywise).
- **No missing values.** Remove incomplete rows or [impute them](/tutorials/imputed-data).
- **One subjects directory.** Every subject, and `fsaverage`, must be in the same directory.
- **The `fsaverage` target only.** The `target` argument exists, but the smoothness estimate, the simulations and the plots always use `fsaverage`. A fix is planned for 0.10.0.
- **A custom `mask` is not used for the smoothness.** The models are fitted only inside your mask, but the smoothness is estimated over the whole cortex.
- **p-values for small samples.** P-values use about 10,000 degrees of freedom whatever the sample size, which makes them too small for samples of a few dozen. See [the model at each vertex](/tutorials/statistics#the-model-at-each-vertex).
- **Unsigned p-value maps.** The p-value maps hold −log<sub>10</sub>(p) without the sign of the effect. Take the direction from the coefficient or t map.
- **Grey-to-white contrast.** `qdecr_w_g.pct` looks for files named like `lh.w_g.pct.fwhm10.fsaverage.mgh`, which FreeSurfer does not write. Link each subject's file to that name first, for both hemispheres:

```bash
for surf in "$SUBJECTS_DIR"/*/surf; do
  for file in "$surf"/?h.w-g.pct*fwhm10.fsaverage.mgh; do
    [ -e "$file" ] || continue
    name=$(basename "$file")
    ln -sf "$name" "$surf/${name:0:2}.w_g.pct.fwhm10.fsaverage.mgh"
  done
done
```

- **missForest results.** Pass the completed data, `result$ximp`, rather than the `missForest` object, which 0.9.0 reads wrongly. `aregImpute` objects are not supported at all.
- **Arguments that do nothing yet.** `mgh`, `clean_up` and `debug` are accepted but not used.
- **Relative paths.** A result stores its paths as given, so a relative `dir_out` only works from the same working directory. See [Saving and loading](/tutorials/saving-and-loading#saving-and-loading-a-result).
