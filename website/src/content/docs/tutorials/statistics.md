---
title: Understanding the statistics
description: 'What QDECR computes at each vertex, and how it corrects for multiple testing across the surface.'
order: 8
---

An analysis has two halves. First QDECR fits a linear regression at every vertex, which gives a coefficient, a standard error, a t-statistic and a p-value per vertex. Then it asks which of those p-values would have been expected from noise alone, given that there are about 150,000 of them and that neighbouring vertices are alike. This tutorial explains both, and the settings that change them.

## The model at each vertex

At every vertex, QDECR fits the same ordinary least squares regression, with the vertex measure as the outcome and your design as the predictors. The design matrix is the same at every vertex; only the outcome changes. So QDECR computes the part of the solution that depends on the design once, and applies it to all vertices in [chunks](/tutorials/performance#chunk-size), which is what makes it fast. With [weights](/tutorials/formulas-and-design#weights), it fits weighted least squares the same way.

For each coefficient, a [stack](/glossary#stack), QDECR keeps four maps:

- **The coefficient**, in the units of the measure per unit of the predictor: millimetres of thickness per year of age, for example.
- **The standard error**, from the residual variance with the usual n − p degrees of freedom.
- **The t-statistic**, the coefficient over its standard error. It is the same t that `lm()` reports.
- **The p-value**, two-sided, stored as −log<sub>10</sub>(p): 3 means p = 0.001.

The analysis covers the cortex of the [target](/glossary#target), without the medial wall: 149,955 vertices of the left hemisphere of `fsaverage`, and 149,926 of the right. Vertices where any subject's value is exactly zero, as at a vertex FreeSurfer left empty, are left out as well.

To analyse part of the cortex, give a mask of your own: `mask`, with one `TRUE` or `FALSE` per vertex of the target, or `mask_path`, an MGH file of ones and zeros. The models are fitted only where the mask is true, and so only there can clusters form. In 0.9.0 the [smoothness](#smoothness) is still estimated over the whole cortex.

> [!NOTE]
> QDECR 0.9.0 computes every p-value from a t-distribution with about 10,000 degrees of freedom, whatever the sample size. The t-statistic is right, but for small samples the p-value comes out smaller than `lm()` would give: with 30 subjects, a t of 2.1 gives p = 0.035 instead of 0.044. From a few hundred subjects on, the difference is negligible. The cause is in how QDECR [pools imputed datasets](#pooling-imputed-data), which it also does for a single one.

## Cluster-wise correction

Testing 150,000 vertices at p < 0.05 would find thousands of them by chance. Correcting each vertex on its own, as Bonferroni would, is far too strict, because neighbouring vertices are not independent tests: the surface data are smooth, so an effect at one vertex is shared with those around it.

QDECR uses FreeSurfer's [cluster-wise correction](/glossary#cluster-wise-correction) instead, the same one `mri_glmfit-sim` applies with its precomputed simulations. For each stack:

1. **Threshold.** Every vertex whose p-value passes the [cluster-forming threshold](/glossary#cluster-forming-threshold), p < 0.001 by default, is marked. Effects in both directions count.
2. **Cluster.** Neighbouring marked vertices join into clusters, and each cluster's size is its area in mm² on the white matter surface.
3. **Compare.** FreeSurfer has simulated Gaussian noise on the `fsaverage` surface thousands of times, at a range of smoothness levels and thresholds, and recorded the largest cluster each simulation produced (Hagler et al., 2006). The [cluster-wise p-value](/glossary#cluster-wise-p-value) of a cluster is the share of simulations whose largest cluster was at least as large.
4. **Keep.** Clusters whose cluster-wise p-value is below `cwp_thr`, 0.025 by default, are [significant](/glossary#significant-cluster). They are what `summary()` lists and the plots show.

FreeSurfer's `mri_surfcluster` does steps 1 to 4, reading the simulations from `$FREESURFER_HOME/average/mult-comp-cor/`.

A cluster is significant as a whole: the correction says that an effect lies somewhere in it, not that every vertex in it has one. A large cluster that spans several regions, which lenient thresholds produce, says little about where the effect is.

This correction assumes that the noise is smooth in the same way everywhere, like the simulated Gaussian noise. Real anatomy is not quite like that. Greve and Fischl (2018) found that, for thickness, it lets through about 10% false positives where 5% was intended, and 20 to 30% for surface area and volume; permutation testing, which QDECR does not offer, did not have the problem. Weigh results for area and volume with that in mind. The references are on the [Cite](/cite#freesurfer-and-the-correction) page.

## Smoothness

Which simulation a cluster is compared with depends on how smooth the data are: in smooth data, large clusters arise by chance more easily. What counts is not the smoothing applied beforehand, the `fwhm` argument, but the smoothness of what the model leaves unexplained. The anatomy itself is smooth, so the residuals are usually smoother than the kernel alone would make them.

After fitting, QDECR writes the residuals as a map per subject and has FreeSurfer's `mris_fwhm` estimate their [smoothness](/glossary#smoothness) over the cortex, as a FWHM in millimetres. It rounds the estimate to a whole millimetre, which is what the simulations are indexed by, and keeps it between 1 and 30 mm, the range they cover.

```r
qdecr_fwhm(out)
```

returns the estimate, and `fwhm.dat` in the output directory holds it too. With imputed data, the residuals are first averaged over the imputed datasets.

## Pooling imputed data

With [imputed data](/tutorials/imputed-data), QDECR fits the model to each imputed dataset at each vertex and pools the fits with Rubin's rules: the mean coefficient, a standard error that adds the spread between the datasets to the uncertainty within them, and a t-statistic and p-value from those. The correction then runs on the pooled maps exactly as above.

For the degrees of freedom of the pooled t, QDECR uses Rubin's large-sample formula, treating the complete data as if they had near-infinite degrees of freedom. A single dataset goes through the same code, and that is where the roughly 10,000 degrees of freedom in the note above come from.

## Choosing thresholds

Two thresholds decide what counts as a finding, and both are set before the analysis.

`mcz_thr` is the cluster-forming threshold. FreeSurfer only simulated six of them, so QDECR accepts those six, written as the p-value, as −log<sub>10</sub>(p), or as FreeSurfer's own number:

| `mcz_thr` | Vertex-wise p | Also accepted |
|---|---|---|
| `0.05` | p < 0.05 | `1.3`, `13` |
| `0.01` | p < 0.01 | `2`, `20` |
| `0.005` | p < 0.005 | `2.3`, `23` |
| `0.001` | p < 0.001, the default | `3`, `30` |
| `0.0005` | p < 0.0005 | `3.3`, `33` |
| `0.0001` | p < 0.0001 | `4`, `40` |

Anything else stops the analysis. FreeSurfer's number appears in the output files' names: `th30` is p < 0.001.

`cwp_thr` is the threshold for the cluster-wise p-value. Its default, 0.025, is 0.05 divided over the two hemispheres, because a whole-brain study is two analyses. Divide further for every further set of tests you make: 0.05 / 4 = 0.0125 for two measures in both hemispheres, for example. Use 0.05 only when your hypothesis concerns one hemisphere alone.

Which cluster-forming threshold to use is a question of what you are looking for. A strict threshold finds focal effects and places them precisely; a lenient one can find weak, widespread effects, as large clusters that locate them poorly. The default, p < 0.001, is the cautious choice. Whichever you choose, choose it before you see the results, and report it with `cwp_thr` and the estimated smoothness.
