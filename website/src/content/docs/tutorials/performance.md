---
title: Performance and memory
description: 'Making an analysis faster and lighter: cores, BLAS, shared memory and chunk size.'
order: 7
---

A vertex-wise analysis repeats one regression about 150,000 times per hemisphere, once for every vertex of the cortex, and once more for every imputed dataset. On a few hundred subjects that takes minutes with the defaults. For thousands of subjects, or dozens of imputed datasets, four settings decide how long it takes and how much memory it needs.

## Cores

`n_cores` sets how many processes QDECR runs at once. It uses them to load the subjects' data, to fit the models, and to summarise the results; FreeSurfer's steps, the smoothness estimate and the cluster correction, run on one core whatever you set.

```r
out <- qdecr_fastlm(..., n_cores = 8)
```

- QDECR accepts at most one core fewer than the machine reports, and stops above that. R counts logical cores, so on a processor with two threads per core, the number of physical cores is usually the better choice.
- More cores pay off most with imputed data, where every vertex is fitted once per dataset. With one dataset, much of the time goes to reading and writing, which extra cores speed up less.
- The processes are forked copies of R, which is why QDECR needs Linux or macOS. Forking from inside RStudio can be fragile, as R's own documentation warns; for long analyses, run a script with `Rscript` from a terminal.

On a computing cluster, set `n_cores` to the cores the job was given: QDECR only checks it against the whole machine.

## BLAS

The regression is matrix algebra, and R hands matrix algebra to a BLAS library. R's own is simple and single-threaded; an optimised one, such as OpenBLAS, speeds up every model QDECR fits. On Debian or Ubuntu:

```bash
sudo apt install libopenblas-dev
```

`sessionInfo()` shows which BLAS R uses, on its `BLAS:` line.

An optimised BLAS usually runs on several threads of its own, and that collides with `n_cores`: eight processes that each start eight threads fight over the same cores and run slower than either alone. Since 0.9.0 QDECR refuses the combination, with this error:

```text
`n_cores` > 1, but there already seems to be a parallel BLAS library present. Either set `n_cores` to 1, or set the BLAS library to 1.
```

So choose one of the two:

- **One process, a threaded BLAS.** Leave `n_cores = 1` and let the BLAS use the cores.
- **Several processes, a single-threaded BLAS.** Tell the BLAS to use one thread before R starts, then set `n_cores`. For OpenBLAS, start R from a shell like this, or put the line in your job script:

```bash
export OPENBLAS_NUM_THREADS=1
Rscript analysis.R
```

It has to be set before R starts: QDECR checks the BLAS in a fresh R process, and the processes it forks inherit the BLAS threads of the session they are copied from. Other libraries use other variables, such as `MKL_NUM_THREADS` for Intel's MKL. Which of the two is faster depends on the machine and the model; for a large study, time both on one hemisphere of a subset first.

## Shared memory

While it runs, QDECR keeps the vertex data of every subject, and the model's residuals, in [file-backed matrices](/glossary#file-backed-matrix): files it reads and writes constantly. Where they live is set with `dir_tmp`, and it defaults to the output directory. On a network drive or a spinning disk, that is the slowest part of the analysis.

On Linux, `/dev/shm` is shared memory: a directory whose files are held in RAM. Pointing `dir_tmp` at it makes the reading and writing as fast as the machine allows:

```r
out <- qdecr_fastlm(..., dir_tmp = "/dev/shm")
```

The files take room in RAM, so there must be enough. For the `fsaverage` target they take about 3.3 MB per subject:

| Subjects | Room needed in `dir_tmp` |
|---|---|
| 500 | 1.7 GB |
| 2,000 | 6.6 GB |
| 10,000 | 33 GB |

`/dev/shm` is usually allowed half the machine's memory. Check how much it has free with:

```bash
df -h /dev/shm
```

- QDECR deletes the files when it finishes, and after an error. If R itself is killed, they stay behind, named after the project, such as `lh.age_sex.thickness_mgh_backend.bk`. Delete them by hand; until then they hold their memory, and a new run of the same project stops because they exist.
- Docker gives a container only 64 MB of `/dev/shm` unless you ask for more with `--shm-size`.
- macOS has no `/dev/shm`. Use the fastest local disk instead.

## Chunk size

QDECR fits the vertices in chunks, 1,000 at a time by default: each process takes a chunk, fits it for every imputed dataset, writes the results and takes the next. `chunk_size` sets how many vertices a chunk holds.

A chunk's working memory grows with the subjects, the imputed datasets and the chunk size together: roughly 8 × subjects × `chunk_size` × (datasets + 2) bytes, for each process. For 5,000 subjects and 20 imputed datasets, the default chunk needs about 880 MB per process, so 7 GB for eight.

```r
# Half the memory per process, for a little more overhead.
out <- qdecr_fastlm(..., data = imp, n_cores = 8, chunk_size = 500)
```

Lower `chunk_size` when an analysis with many subjects or datasets runs out of memory. Raising it above the default gains little.
