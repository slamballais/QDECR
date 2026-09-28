---
title: Saving, loading and MGH files
description: 'Saving and loading results, freeing memory, and reading and writing MGH files.'
order: 6
---

An analysis leaves two things behind: the [result](/glossary#result) in R, and the [output directory](/get-started#what-gets-written-to-disk) with every map. The result is small because the maps stay on disk; it knows where they are. This tutorial covers keeping results between sessions, getting at the vertex data and the maps, and working with FreeSurfer's MGH files directly.

## Saving and loading a result

QDECR saves every result as it finishes, as an `.rds` file in its output directory: `results/lh.age_sex.thickness/lh.age_sex.thickness.rds`. In a later session, load it by its directory or by the file:

```r
out <- qdecr_load("results/lh.age_sex.thickness")
```

The loaded result works like the original: `summary()`, `stacks()`, the plots and the functions below all read the maps from the output directory.

> [!NOTE]
> A result stores its paths as they were given. With a relative `dir_out`, such as `"results"`, load it from the same working directory; after moving the output directory, the paths no longer point at the maps. An absolute `dir_out` avoids both.

By default the saved file includes the data and the design matrices, so the analysis can be traced back completely. With many imputed datasets that makes the file large, and it holds your participants' data, which matters before you share it. `save_data = FALSE` leaves them out, and `save = FALSE` skips the file altogether:

```r
out <- qdecr_fastlm(..., save_data = FALSE)
```

To save a result again, with or without its data, use `qdecr_save()`. It writes to the output directory under the project's name, unless you give another with `file`:

```r
qdecr_save(out, save_data = FALSE)
```

## Freeing memory

During an analysis QDECR keeps the vertex data of every subject in a [file-backed matrix](/glossary#file-backed-matrix) in `dir_tmp`, and deletes it at the end. The result keeps only summaries of the vertex data, such as each vertex's mean, which is what `hist()` plots.

For analyses of your own on the vertex data, `reload()` builds that matrix again: rows are vertices, columns are subjects in the order of the data. It reads every subject's file from the subjects directory, so it takes as long as the loading step of the analysis did, and it needs the same room in `dir_tmp`.

```r
out <- reload(out)
out$mgh
```

`qdecr_load(path, reload = TRUE)` loads and reloads in one step. When you are done, `unload()` deletes the matrix's file and removes it from the result:

```r
out <- unload(out)
```

With `dir_tmp = "/dev/shm"` that file lives in memory, so unload when you are done, before R's session ends. A file left behind in `/dev/shm` holds its memory until the machine restarts or you delete it.

The `fbm_*` helpers summarise a file-backed matrix by row or by column: `fbm_row_mean()`, `fbm_row_sd()` and `fbm_row_sum()` per vertex, and `fbm_col_mean()`, `fbm_col_sd()` and `fbm_col_sum()` per subject. `row.mask` and `col.mask` limit them to some vertices or some subjects. The mean thickness of each subject inside a significant cluster, for a follow-up analysis, is a column mean over that cluster's vertices:

```r
out <- reload(out)
in_cluster <- qdecr_read_ocn(out, "age")$x == 1
pheno$cluster_thickness <- fbm_col_mean(out$mgh, row.mask = in_cluster)
out <- unload(out)
```

The columns follow the rows of the data the analysis used, so the new column lines up with `pheno` if that is the data frame you passed.

## Reading the maps

The `qdecr_read_*` functions read one map of one stack from the output directory. The stack can be named or numbered, as in `stacks(out)`:

| Function | Reads |
|---|---|
| `qdecr_read_coef(out, "age")` | The coefficient at every vertex. |
| `qdecr_read_se(out, "age")` | Its standard error. |
| `qdecr_read_t(out, "age")` | The t-statistic. |
| `qdecr_read_p(out, "age")` | The p-value, as −log<sub>10</sub>(p). |
| `qdecr_read_ocn(out, "age")` | The [cluster map](/glossary#cluster-map): each vertex's cluster number, 0 outside the significant clusters. |
| `qdecr_read_ocn_mask(out, "age")` | `TRUE` for every vertex in a significant cluster. |

Each returns an MGH object, below, whose `x` holds one value per vertex. The p-value map has no sign, unlike the maps FreeSurfer's own `mri_glmfit` writes; take the direction of an effect from the coefficient or the t-statistic.

```r
coef <- qdecr_read_coef(out, "age")$x
significant <- qdecr_read_ocn_mask(out, "age")
sum(significant)                # vertices in significant clusters
range(coef[significant])        # the effects found there
```

A few more functions describe the analysis as a whole: `qdecr_fwhm(out)` returns the estimated [smoothness](/tutorials/statistics#smoothness), `nobs(out)` the number of subjects, and `formula(out)` the model.

## MGH files

MGH is FreeSurfer's format for surface and volume data. QDECR reads subjects' maps from it and writes its results in it, and three functions let you do the same:

- `load.mgh(path)` reads a file into a list: `x` holds the values, and the other elements the header (`ndim1`, the number of vertices, and so on).
- `as_mgh(x)` makes the same list from a numeric vector, one value per vertex.
- `save.mgh(mgh, path)` writes such a list to a file.

Together they turn anything you compute per vertex into a map for Freeview or FreeSurfer's tools. The effect of age per decade rather than per year, for example:

```r
decade <- as_mgh(qdecr_read_coef(out, "age")$x * 10)
save.mgh(decade, "results/lh.age_per_decade.mgh")
```

These functions handle the kind of MGH file FreeSurfer writes for surface maps: uncompressed, with 32-bit values. They do not read compressed `.mgz` files.

Two more read and write related formats:

- `load.annot(path)` reads a FreeSurfer annotation such as `lh.aparc.annot`: `vd_label` gives each vertex's region code, and `LUT` the table from codes to region names and colours. `summary(out, annot = TRUE)` uses it to name the regions a cluster covers.
- `bsfbm2mgh(fbm, path)` writes rows of a file-backed matrix as the frames of an MGH file. QDECR uses it to split its results into one file per stack.
