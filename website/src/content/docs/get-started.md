---
title: Get started
description: 'What QDECR needs, how to install it on Linux, macOS, Windows or a computing cluster, and a first analysis.'
order: 0
---

QDECR runs a [vertex-wise analysis](/glossary#vertex-wise-analysis): the same regression model at each of the 163,842 points of FreeSurfer's cortical surface, followed by FreeSurfer's own correction for testing that many points at once. You write the model as an R formula, `qdecr_thickness ~ age + sex`, and QDECR takes care of loading the surfaces, fitting the models and finding the significant clusters.

It is for researchers who have already run their MRI scans through FreeSurfer and analyse their data in R. If that is you, this page takes you from nothing installed to a first result.

## Requirements

- **Linux or macOS**, or Windows through WSL2. QDECR runs its parallel work in forked processes, which Windows does not have, and it calls FreeSurfer's command-line tools.
- **FreeSurfer 6.0 or later**, set up so that `FREESURFER_HOME` points at it. QDECR calls `mris_fwhm` and `mri_surfcluster`, reads FreeSurfer's precomputed simulations, and opens Freeview for plots. FreeSurfer needs its free licence file, which you get when you [register](https://surfer.nmr.mgh.harvard.edu/registration.html).
- **Your subjects processed with `recon-all`, including `-qcache`**, which resamples each subject's measures to the `fsaverage` target and smooths them. QDECR reads the files it writes.
- **R 3.6 or later**: the current releases of the packages QDECR builds on need it, and they compile C++ code, so R needs a compiler too. A current R 4 release is best.
- **Optional:** the `magick` R package, for the snapshots `qdecr_snap()` makes, and a fast [BLAS library](/tutorials/performance#blas).

If your subjects were processed without `-qcache`, run it on its own for each one:

```bash
recon-all -s sub-001 -qcache
```

## Install

In every case the steps are the same: FreeSurfer, then R, then QDECR. What differs is how you get the first two.

### Linux

Install FreeSurfer with the [official instructions](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall) and put your licence file where they say. Then make FreeSurfer part of every shell by adding two lines to `~/.bashrc`, with the path where you installed it:

```bash
export FREESURFER_HOME=/usr/local/freesurfer
source "$FREESURFER_HOME/SetUpFreeSurfer.sh"
```

Install R from your distribution or from [CRAN](https://cran.r-project.org/bin/linux/). QDECR's dependencies compile C++ code, so R needs its development tools too, and pak, the installer below, needs libcurl's development files. On Debian or Ubuntu:

```bash
sudo apt install r-base r-base-dev libcurl4-openssl-dev
```

For `qdecr_snap()`, add ImageMagick's development files (`libmagick++-dev` on Debian and Ubuntu) and then `install.packages("magick")` in R.

Then install QDECR itself from GitHub, in R:

```r
install.packages("pak")
pak::pak("slamballais/QDECR")
```

That installs the latest release and every package it needs, which takes a few minutes the first time. To install one release in particular, name its tag: `pak::pak("slamballais/QDECR@0.9.0")`. If you already use `remotes`, `remotes::install_github("slamballais/QDECR")` does the same.

> [!NOTE]
> QDECR finds FreeSurfer through the environment R starts in. Started from a shell that has run `SetUpFreeSurfer.sh`, R has it. Started from a desktop launcher, as RStudio usually is, it does not: see [R can't find FreeSurfer](/tutorials/troubleshooting#r-cant-find-freesurfer).

### macOS

FreeSurfer and R both run on macOS, and so should QDECR, but we have not tested it there: if you try it, we would like to [hear how it went](/help#where-to-ask). Install FreeSurfer with the [official instructions](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall), R from [CRAN](https://cran.r-project.org/bin/macosx/), and Apple's command-line tools for the compiler (`xcode-select --install`). Then install QDECR as on Linux.

macOS has no `/dev/shm`, the shared memory that makes QDECR much faster on Linux. Leave `dir_tmp` at its default, or point it at the fastest disk you have.

### Windows, through WSL2

On Windows, QDECR runs inside the Windows Subsystem for Linux: a real Linux, Ubuntu by default, next to Windows. In PowerShell:

```powershell
wsl --install
```

After a restart, open Ubuntu from the Start menu and follow the Linux steps above inside it: FreeSurfer, R and QDECR all go into Ubuntu, not into Windows.

- Keep FreeSurfer's output inside the Linux file system, such as in your Ubuntu home directory, rather than under `/mnt/c/`. WSL reads Windows drives far more slowly.
- On Windows 11, Freeview's windows open on the Windows desktop by themselves, so `freeview()` and `qdecr_snap()` work as they do on Linux.
- An R installed on the Windows side cannot run QDECR, even pointed at the same files.

### A computing cluster

On a cluster, FreeSurfer and R are usually modules and QDECR goes into your own R library. The names differ between clusters; ask `module avail`:

```bash
module load freesurfer R
```

Then install QDECR once, from R on a login node, as on Linux. Run analyses as batch jobs with `Rscript`, and keep three things in step with what the job asked for:

- **`n_cores`** should match the cores the job was given, such as `SLURM_CPUS_PER_TASK` under Slurm. QDECR only checks it against the whole machine.
- **`dir_tmp`** should be on fast storage local to the node: `/dev/shm` if the cluster allows it, or the job's scratch directory.
- **Parallel BLAS**: many clusters' R uses a multithreaded BLAS, and QDECR refuses `n_cores` above 1 alongside it. Set `OPENBLAS_NUM_THREADS=1` in the job script before R starts; [Performance and memory](/tutorials/performance#blas) explains why.

## Check that it works

These checks take a minute and catch most problems before an analysis does. In R, from the environment you will run QDECR in:

```r
library(QDECR)
Sys.getenv(c("FREESURFER_HOME", "SUBJECTS_DIR"))
Sys.which(c("mris_fwhm", "mri_surfcluster", "freeview"))
```

Both variables should name directories, and all three programs should have a path. Then check that the `fsaverage` target sits in your subjects directory, and that a subject has the files `-qcache` writes:

```bash
ls "$SUBJECTS_DIR/fsaverage/surf/lh.inflated"
ls "$SUBJECTS_DIR/sub-001/surf/" | grep fwhm10.fsaverage
```

The second should list files like `lh.thickness.fwhm10.fsaverage.mgh`. If `fsaverage` is missing, which happens when the subjects were processed elsewhere, link it in from FreeSurfer:

```bash
ln -s "$FREESURFER_HOME/subjects/fsaverage" "$SUBJECTS_DIR/"
```

## Prepare the data

QDECR takes an ordinary data frame, with one row per subject:

- **An ID column** whose values are the names of the subjects' directories in `SUBJECTS_DIR`: the row for `sub-001` belongs to `$SUBJECTS_DIR/sub-001`. You tell QDECR its name with `id`.
- **Your variables**, with categorical ones as factors, so that R codes them the way you mean. Which level is the reference decides the names of the results: with `female` as the reference, the effect of sex is called `sexmale`.
- **No missing values** in the variables the model uses. QDECR stops rather than drop subjects without saying so. Remove incomplete rows yourself, or [impute them](/tutorials/imputed-data).

```r
pheno <- read.csv("phenotypes.csv")
pheno$sex <- factor(pheno$sex, levels = c("female", "male"))
pheno <- pheno[complete.cases(pheno[, c("id", "age", "sex")]), ]
```

## A first analysis

This fits `qdecr_thickness ~ age + sex` at every vertex of the left hemisphere: the effect of age on cortical thickness, adjusted for sex.

```r
out <- qdecr_fastlm(
  qdecr_thickness ~ age + sex,
  data = pheno,
  id = "id",
  hemi = "lh",
  project = "age_sex",
  dir_out = "results",
  dir_tmp = "/dev/shm"
)
```

- `qdecr_thickness` is the [vertex measure](/tutorials/formulas-and-design#the-vertex-measure): the name of a FreeSurfer surface map with `qdecr_` in front.
- `project` names the analysis. QDECR adds the hemisphere and the measure, so the results go to `results/lh.age_sex.thickness/`. `dir_out` must exist or be creatable; it defaults to the working directory.
- `dir_tmp = "/dev/shm"` keeps the large temporary files in memory, which makes a big difference on Linux. Leave it out on macOS.

QDECR reports each stage as it goes: checking the input, loading the vertex data of every subject, fitting the models, estimating the smoothness of the residuals, and the cluster-wise correction. For a few hundred subjects on one core that takes minutes; [Performance and memory](/tutorials/performance) covers larger studies.

When it finishes, `out` holds the [result](/glossary#result). Print it for the settings and the sample, and summarise it for the significant clusters and the regions they cover:

```r
out
summary(out, annot = TRUE)
```

A whole-brain study runs the right hemisphere as well, with `hemi = "rh"`. The default cluster-wise threshold, `cwp_thr = 0.025`, already splits 0.05 over the two hemispheres.

From here, the [quick start](/tutorials/quick-start) walks through a complete analysis of real data, with its output, and the [tutorials](/tutorials) take each part further.

## What gets written to disk

Everything an analysis produces stays in its output directory, `results/lh.age_sex.thickness/` above. The result in R knows the paths, so you rarely need to open these files yourself, but FreeSurfer's tools and other software can read all of them.

For each [stack](/glossary#stack), one per coefficient of the model and numbered in the order of `stacks(out)`:

| File | What it holds |
|---|---|
| `stack2.coef.mgh` | The coefficient at every vertex. |
| `stack2.se.mgh` | Its standard error. |
| `stack2.t.mgh` | The t-statistic. |
| `stack2.p.mgh` | The p-value, as −log<sub>10</sub>(p). |
| `stack2.cache.th30.abs.sig.cluster.summary` | FreeSurfer's table of the significant clusters. |
| `stack2.cache.th30.abs.sig.ocn.mgh` | The [cluster map](/glossary#cluster-map): each vertex numbered by its cluster. |
| `stack2.cache.th30.abs.sig.ocn.annot` | The same clusters as an annotation, for Freeview. |
| `stack2.cache.th30.abs.sig.cluster.mgh` | Each cluster's cluster-wise p-value, as −log<sub>10</sub>(p), on its vertices. |
| `stack2.cache.th30.abs.sig.masked.mgh` | The p-value map, kept only inside significant clusters. |
| `stack2.cache.th30.abs.sig.voxel.mgh` | Vertex-wise p-values, corrected for the whole surface. |

`th30` is the [cluster-forming threshold](/tutorials/statistics#choosing-thresholds): 30 stands for p < 0.001. For the analysis as a whole:

| File | What it holds |
|---|---|
| `lh.age_sex.thickness.rds` | The saved result, for [`qdecr_load()`](/tutorials/saving-and-loading#saving-and-loading-a-result). |
| `stack_names.txt` | Each stack's number and name. |
| `significant_clusters.txt` | The output of `summary(out, annot = TRUE)`. |
| `finalMask.mgh` | The vertices the analysis covered. |
| `fwhm.dat` | The estimated [smoothness](/tutorials/statistics#smoothness) of the residuals. |

Two arguments change this layout. With `dir_out_tree = FALSE`, the files go straight into `dir_out` rather than a directory of their own. `dir_out` must then be a new directory: QDECR will not write into one that exists, and refuses `clobber = TRUE` in this case, so that an analysis can never delete a directory it did not make. `file_out_tree` puts the project's full name in front of each file, as in `lh.age_sex.thickness.stack2.coef.mgh`. It is on whenever `dir_out_tree` is off, and off otherwise, unless you set it.

While it runs, QDECR also writes large temporary files to `dir_tmp`: the vertex data of every subject and the residuals, about 3.3 MB per subject together, so 3.3 GB for 1,000 subjects. They are deleted at the end unless you set `clean_up_bm = FALSE`. [Performance and memory](/tutorials/performance#shared-memory) has the details.
