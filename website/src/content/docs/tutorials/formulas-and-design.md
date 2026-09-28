---
title: Formulas and design
description: 'How QDECR reads a model formula: the vertex measure, covariates, factors, interactions and weights.'
order: 2
---

QDECR reads its model the way `lm()` does: as an R formula, turned into a design matrix by R itself. Everything R's formulas can say about the right-hand side, QDECR can use, and the left-hand side names the surface map to analyse.

```r
qdecr_thickness ~ age + sex
```

## The vertex measure

The left-hand side is the [vertex measure](/glossary#vertex-measure): the name of a FreeSurfer surface map with `qdecr_` in front. QDECR recognises these:

| In the formula | What it is |
|---|---|
| `qdecr_thickness` | Cortical thickness. |
| `qdecr_area` | Surface area, of the white matter surface. |
| `qdecr_area.pial` | Surface area, of the pial surface. |
| `qdecr_volume` | Grey matter volume. |
| `qdecr_curv` | Curvature of the white matter surface. |
| `qdecr_sulc` | Sulcal depth. |
| `qdecr_white.H` | Mean curvature of the white matter surface. |
| `qdecr_white.K` | Gaussian curvature of the white matter surface. |
| `qdecr_jacobian_white` | How much registration to the target stretched the white matter surface. |
| `qdecr_pial` | The map named `pial`. |
| `qdecr_pial_lgi` | Local gyrification index, which needs `recon-all -localGI` first. |
| `qdecr_w_g.pct` | Grey-to-white contrast. In 0.9.0 it needs a [workaround](/tutorials/troubleshooting#known-limitations). |

The name tells QDECR which file to read for each subject: `qdecr_thickness`, with `hemi = "lh"` and the defaults, reads `surf/lh.thickness.fwhm10.fsaverage.mgh` in the subject's directory, one of the files `recon-all -qcache` writes.

The `fwhm` argument picks the smoothing: 10 mm by default, and 5 mm for `qdecr_pial_lgi`. `-qcache` smooths at 0, 5, 10, 15, 20 and 25 mm, so any of those works. `fwhm = 0` reads the unsmoothed file, `lh.thickness.fsaverage.mgh`.

A formula has exactly one vertex measure, on its left-hand side, as it is. QDECR stops if the measure is transformed, as in `log(qdecr_thickness)`, appears twice, or appears among the predictors. To analyse a transformed map, write the transformed values to MGH files of their own and use them as a custom measure.

### Your own maps

Any other surface map works too, if it is stored like FreeSurfer's: one MGH file per subject, in the subject's `surf` directory, resampled to the target and named the same way. Give its name, with the `qdecr_` prefix, to `custom_measure` as well as to the formula:

```r
# Reads surf/lh.myelin.fwhm10.fsaverage.mgh for each subject.
out <- qdecr_fastlm(qdecr_myelin ~ age + sex, custom_measure = "qdecr_myelin", ...)
```

## Covariates and factors

The right-hand side is an ordinary R formula, read against your data frame. Continuous variables enter as they are; `+` adds one; `- 1` removes the intercept. QDECR fits the same model at every vertex, so a covariate adjusts every vertex for the same thing.

Each column of the design matrix becomes a [stack](/glossary#stack): its coefficient, standard error, t and p maps, and its own clusters. `stacks()` lists them, in order:

```r
out <- qdecr_fastlm(qdecr_thickness ~ age + sex, ...)
stacks(out)
#> [1] "(Intercept)" "age"         "sexmale"
```

A factor contributes one column per level beyond the first, with R's default treatment coding. A `sex` factor with levels `female` and `male` becomes the stack `sexmale`, the difference of men from women. A factor `status` with levels `control`, `MCI` and `AD` becomes two stacks, `statusMCI` and `statusAD`, each compared with the controls. To compare with another level, change the reference before the analysis:

```r
pheno$status <- relevel(pheno$status, ref = "MCI")
```

The stacks are numbered in the order of `stacks()`, so `stacks(out)[2]` is the effect of age, and functions that take a stack accept either: `qdecr_snap(out, "age")` and `qdecr_snap(out, 2)` draw the same map.

The design matrix must have full rank: no column may be a combination of others. QDECR checks this before loading any vertex data and stops with "The design matrix is NOT full rank" if it is not, which usually means two variables carry the same information, such as a dummy for every level of a factor.

## Interactions and transformations

Whatever R can build into a design matrix, QDECR can fit:

| Formula | Stacks it adds |
|---|---|
| `qdecr_thickness ~ age * sex` | `age`, `sexmale` and the interaction `age:sexmale` |
| `qdecr_thickness ~ poly(age, 2)` | Orthogonal linear and quadratic age: `poly(age, 2)1`, `poly(age, 2)2` |
| `qdecr_thickness ~ age + I(age^2)` | `age` and `I(age^2)`, the raw square |
| `qdecr_thickness ~ splines::ns(age, df = 3)` | Three spline columns for a smooth curve of age |
| `qdecr_thickness ~ scale(age)` | Age in standard deviations, so the coefficient is per SD |
| `qdecr_thickness ~ I(weight / height^2)` | A value computed in the formula, here the body mass index |

The stack names are the design matrix's column names, which for splines and polynomials are long. `stacks(out)` shows them, and the stack number is often easier to type.

> [!WARNING]
> Each stack is tested and corrected on its own. Every term you add is another map that can show clusters, and nothing adjusts for how many you looked at. Choose the model, and the stacks you will interpret, before you run it.

To run the same model for several exposures, build the formula in code and give each run its own project name:

```r
for (exposure in c("bmi", "income", "sleep")) {
  f <- reformulate(c(exposure, "age", "sex"), response = "qdecr_thickness")
  qdecr_fastlm(f, data = pheno, id = "id", hemi = "lh", project = exposure, dir_tmp = "/dev/shm")
}
```

## Weights

`weights` gives each subject a weight in the regression, as `lm()` does: a subject with weight 2 counts twice as much as one with weight 1. Use it for inverse probability weights, for example, to correct for who took part in the scan.

Unlike `lm()`, QDECR takes the weights themselves, not the name of a column. Pass a numeric vector with one positive value per row of the data, in the same order:

```r
out <- qdecr_fastlm(
  qdecr_thickness ~ age + sex,
  data = pheno,
  weights = pheno$ipw,
  id = "id", hemi = "lh", project = "age_sex_weighted"
)
```

With [imputed data](/tutorials/imputed-data), the same weights apply to every imputed dataset, so they cannot themselves have missing values. `print(out)` shows whether weights were used.

## Other models

`qdecr_fastlm()` is QDECR's model today: linear regression, fitted with the same algebra at every vertex. QDECR is built so that others can be added: its `prep_fun` and `analysis_fun` arguments name the functions that build the design and fit it at each vertex, and the checks, data loading, correction and output around them stay the same. The [paper](/cite) describes those modules. They are meant for developers, and the defaults are the only values that do anything useful today.
