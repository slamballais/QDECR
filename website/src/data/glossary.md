# QDECR

QDECR fits a statistical model at every vertex of FreeSurfer's cortical surfaces, from R, and corrects for multiple testing with FreeSurfer's own cluster-wise simulations. This is the project's glossary: the code, the documentation and the Glossary page on qdecr.com use these words, and the site's page is built from this file.

## Language

### Surfaces and FreeSurfer data

**Vertex**:
A point of the triangle mesh that FreeSurfer fits to the cortex. Every surface measure has one value per vertex; the `fsaverage` target has 163,842 vertices per hemisphere.
_Avoid_: node

**Hemisphere**:
The left (`lh`) or right (`rh`) half of the cortex. One analysis covers one hemisphere, so a whole-brain study is two analyses.

**Vertex measure**:
The surface map analysed at every vertex, such as cortical thickness or surface area. In a formula it carries the `qdecr_` prefix: `qdecr_thickness`.
_Avoid_: vertex variable, surface variable

**Subjects directory**:
The directory that holds one FreeSurfer output directory per subject, named by subject ID, and the target. FreeSurfer calls it `SUBJECTS_DIR`; QDECR's argument is `dir_subj`.

**Target**:
The template subject that every subject's surface data is resampled to, so that a vertex means the same place for everyone. `fsaverage` by default.
_Avoid_: template, average subject

**qcache**:
The `recon-all -qcache` step of FreeSurfer, which resamples each subject's surface measures to the target and smooths them. QDECR reads the files it writes.

**FWHM**:
Full width at half maximum: the width, in millimetres, of the smoothing kernel applied to the surface data before the analysis. 10 mm by default, 5 mm for local gyrification.
_Avoid_: smoothing kernel size

**Mask**:
The vertices an analysis covers. By default the cortex of the target, without the medial wall, less any vertex where a subject's value is exactly zero.

**MGH file**:
FreeSurfer's binary format for surface and volume data, with the extension `.mgh`. QDECR reads subject data from it and writes every result map in it.

**Annotation**:
A division of the target's surface into named regions, such as the Desikan-Killiany atlas in `aparc.annot`. QDECR uses it to say which regions a cluster covers.
_Avoid_: parcellation file

### Models

**Vertex-wise analysis**:
The same statistical model fitted separately at every vertex, with the vertex measure as its outcome.
_Avoid_: vertexwise, mass-univariate analysis

**Formula**:
The model, in R's formula notation, with the vertex measure on the left of the tilde: `qdecr_thickness ~ age + sex`.

**Design matrix**:
The matrix R builds from the formula and the data, with one column per coefficient of the model.

**Stack**:
One coefficient of the model, a column of the design matrix, together with every map QDECR writes for it. Stacks are named after their column (`age`, `sexmale`) and numbered in the matrix's order, so stack 1 is usually the intercept.
_Avoid_: contrast

**Weights**:
Observation weights for the regression, one per subject, as in R's `lm`.

**Imputed datasets**:
Several copies of a dataset, each with its missing values filled in differently by multiple imputation. QDECR fits the model to every copy and pools the results.
_Avoid_: imputations

**Pooling**:
Combining the estimates from each imputed dataset into one coefficient, standard error and p-value per vertex, with Rubin's rules.

### Correction for multiple testing

**Cluster-wise correction**:
The correction for multiple testing QDECR applies: vertices that pass a threshold are grouped into clusters, and each cluster is judged by its size against FreeSurfer's precomputed Monte Carlo simulations of pure noise.
_Avoid_: MCZ correction, Monte Carlo correction

**Cluster-forming threshold**:
The vertex-wise p-value a vertex must pass to join a cluster, set with `mcz_thr`: 0.001 by default.
_Avoid_: MCZ threshold, vertex-wise threshold

**Smoothness**:
How smooth the residuals of the model are across the surface, estimated after the fit and expressed as a FWHM in millimetres. It picks which of FreeSurfer's simulations the clusters are compared with.
_Avoid_: estimated FWHM

**Cluster**:
A connected set of vertices that all pass the cluster-forming threshold, for one stack.

**Cluster-wise p-value**:
The probability of a cluster at least as large as the one found, if there were no effect anywhere. A cluster is kept when this is below `cwp_thr`: 0.025 by default, 0.05 split over two hemispheres.
_Avoid_: CWP

**Significant cluster**:
A cluster whose cluster-wise p-value is below `cwp_thr`. Only significant clusters appear in the summary and the plots.

**Cluster map**:
The map that numbers each vertex by the significant cluster it belongs to, and is zero elsewhere. FreeSurfer calls it the output cluster number (OCN) map.
_Avoid_: OCN map

### Results

**Project**:
The name of one analysis. QDECR combines it with the hemisphere and the measure, as in `lh.age_sex.thickness`, for the output directory and file names.

**Output directory**:
Where an analysis writes its maps, cluster tables and saved result: a directory named after the project, inside `dir_out`.

**Result**:
The R object an analysis returns, of class `vw`: the settings, paths, model and summaries. The maps themselves stay in the output directory, and the result knows where.
_Avoid_: vw object, output object

**File-backed matrix**:
A matrix stored in a file rather than in memory, which QDECR uses for the vertex data of every subject and for the model's results while it runs. Built on the bigstatsr package, which calls it an FBM.
_Avoid_: big matrix, FBM

**Temporary directory**:
Where the file-backed matrices are written during an analysis, set with `dir_tmp`. Fast storage here, such as shared memory, makes the analysis much faster.
