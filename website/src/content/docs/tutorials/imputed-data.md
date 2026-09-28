---
title: Imputed data
description: 'Running QDECR on multiply imputed data, which inputs it accepts, and how the results are pooled.'
order: 3
---

Large studies rarely have every covariate for every subject. Dropping the incomplete rows costs power and can bias the result; multiple imputation fills in the missing values several times over, and the analysis is run on each [imputed dataset](/glossary#imputed-datasets) and pooled. QDECR does the running and the pooling for you: give it the imputation instead of a data frame.

Only the covariates are imputed. The vertex data come from FreeSurfer and are complete for every subject with a scan, and they are the same in every imputed dataset.

## Supported inputs

The `data` argument takes any of these:

| Input | Made by | What QDECR uses |
|---|---|---|
| A data frame | | The data, as one dataset. |
| A matrix | | The same, turned into a data frame. |
| A list of data frames | Anything | Each data frame as one imputed dataset. |
| A `mids` object | [mice](https://cran.r-project.org/package=mice) | Each of its `m` completed datasets. |
| An `amelia` object | [Amelia](https://cran.r-project.org/package=Amelia) | Its `imputations`. |
| An `mi` object | [mi](https://cran.r-project.org/package=mi) | The completed datasets from `mi::complete()`. |
| A `missForest` result | [missForest](https://cran.r-project.org/package=missForest) | Pass its `ximp` data frame rather than the object itself, which 0.9.0 reads wrongly. |
| An `aregImpute` object | [Hmisc](https://cran.r-project.org/package=Hmisc) | Not supported: it does not keep the observed data. Build the completed datasets yourself and pass them as a list. |

`imp2list()` does the conversion, and you can call it yourself to see what QDECR will receive: a list with one data frame per imputed dataset.

```r
datasets <- imp2list(imp)
length(datasets)
head(datasets[[1]])
```

## Running on imputed data

Impute as you normally would, keeping the ID column in the data but out of the model. With mice:

```r
library(mice)

predictors <- make.predictorMatrix(pheno)
predictors[, "id"] <- 0

imp <- mice(pheno, m = 20, predictorMatrix = predictors, seed = 2026, printFlag = FALSE)
```

Then pass the `mids` object as `data`:

```r
out <- qdecr_fastlm(
  qdecr_thickness ~ age + sex + income,
  data = imp,
  id = "id",
  hemi = "lh",
  project = "income",
  dir_tmp = "/dev/shm"
)
```

`print(out)` reports the number of datasets under `data`, next to the number of subjects. The rest of QDECR works as it does for one dataset: the stacks, the summary, the plots.

> [!WARNING]
> QDECR loads each subject's vertex data once, in the order of the IDs in the first dataset, and pairs it with the rows of every dataset in that order. The imputation packages keep the rows in place, but a list you build yourself must have the same subjects, in the same order, in every data frame. QDECR does not check.

Every imputed dataset adds a model fit at every vertex, so an analysis takes longer and needs more memory as `m` grows. [Performance and memory](/tutorials/performance) covers what to expect.

## How the results are pooled

At each vertex, QDECR fits the model to every imputed dataset and combines the fits with Rubin's rules, the same rules `mice::pool()` applies:

- **The coefficient** is the mean of the coefficients over the datasets.
- **Its variance** is the mean of the squared standard errors, the variance within the datasets, plus the variance of the coefficients between them, times (1 + 1/m).
- **The t-statistic** is the pooled coefficient over the pooled standard error, and its p-value uses Rubin's degrees of freedom, which grow as the datasets agree more.

The pooled coefficient, standard error, t and p maps are what QDECR writes to disk and corrects. For the correction it needs the residuals as well: it takes the mean of each subject's residuals over the datasets, and estimates the [smoothness](/tutorials/statistics#smoothness) from those.

The degrees of freedom leave out the small-sample adjustment of Barnard and Rubin, which `mice::pool()` includes: QDECR treats the complete data as if they had near-infinite degrees of freedom. With the thousands of subjects QDECR was written for this makes no practical difference; with a few dozen, the p-values come out somewhat smaller than `mice::pool()` would give. [Understanding the statistics](/tutorials/statistics#the-model-at-each-vertex) says more.
