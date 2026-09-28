# The example analysis

Every output the site shows, printed results, cluster tables, histograms, snapshots and
the interactive viewer, comes from one real analysis: cortical thickness against age and
sex in the typically developing participants of one ABIDE site, run with QDECR on
FreeSurfer 7.4.1 inside WSL2.
The scripts here reproduce it from nothing installed, on Ubuntu (in WSL2 or not). The
downloaded data and the FreeSurfer install stay out of the repository; what they
produce is committed.

Run them from the repository root, in this order, inside Ubuntu:

| Script                | Needs              | Does                                                                                                                                                                        |
| --------------------- | ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `setup-system.sh`     | `sudo`             | Installs R, its compilers, libcurl's and ImageMagick's headers, Xvfb, and the shared libraries FreeSurfer's binaries and Freeview's Qt load.                                 |
| `setup-freesurfer.sh` | 10 GB of disk      | Fetches the FreeSurfer 7.4.1 tarball (9.5 GB, resumable) and extracts a slim install into `~/qdecr-example/freesurfer`: the two tools QDECR calls, Freeview, fsaverage, the simulations. |
| `setup-r.sh`          |                    | Installs pak, QDECR from GitHub, magick and jsonlite into the user's R library.                                                                                              |
| `download.sh`         | R                  | Fetches ABIDE's phenotype file, picks the subjects (`select-subjects.R`), fetches their 198 thickness maps (130 MB), links fsaverage in.                                     |
| `run.sh`              | the licence        | Runs `run.R` for each hemisphere on a virtual display: the analysis, then `print`, `summary`, `hist` and `qdecr_snap`, keeping what R printed.                               |
| `export.sh`           | a run              | Runs `export.R`: writes `run.json`, the printed output, the figures and the viewer's downsampled maps into the site, and reports their sizes.                                 |

`env.sh` holds the paths they share; set `QDECR_EXAMPLE_ROOT` to put everything somewhere
other than `~/qdecr-example`.

## The licence

FreeSurfer's tools refuse to run without its licence file, which is free but personal:
[register](https://surfer.nmr.mgh.harvard.edu/registration.html), and save the
`license.txt` that arrives by email as `~/qdecr-example/freesurfer/license.txt` (or
put it at `~/license.txt` before `setup-freesurfer.sh`, which copies it). Nothing before
`run.sh` needs it.

## What lands in the site

- `src/data/example/run.json`: the run, in the shape `src/lib/example.ts` checks: the
  sample, the software, the model, and per hemisphere the stacks and every significant
  cluster with its size, cluster-wise p-value, peak and regions.
- `src/data/example/subjects.csv`: the data frame the analysis read (id, age, sex, site).
- `src/data/example/output/`: what R printed, for the site's output blocks.
- `src/assets/example/`: the histograms and the `qdecr_snap()` images.
- `public/viewer/`: fsaverage6's inflated surface and curvature, and the age stack's
  t-statistic and cluster maps on it, for the viewer.

## Credit

The data are from the [Autism Brain Imaging Data Exchange](https://fcon_1000.projects.nitrc.org/indi/abide/)
(ABIDE I), preprocessed and shared by the [Preprocessed Connectomes Project](http://preprocessed-connectomes-project.org/abide/),
under [CC BY-NC-SA 3.0](https://creativecommons.org/licenses/by-nc-sa/3.0/). Everything
derived from them here carries that licence, whatever the site's own.
