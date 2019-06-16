# QDECR 0.7.1: OHBM (.1)

Version 0.7.1 is the first version of QDECR after the public release. It will feature several improvements as well as additional arguments that extend (but not modify) the user experience.

## New functions and features

* Added the input argument `dir_target`, so that the target can be specified flexibly. 
* Modularized `qdecr_model` and added the input argument `prep_fun`, so that users can choose and create their own prep functions.

## Major changes

* The internal function to run vertex-wise analyses, called `vertexwise`, now processes regressions in chunks, i.e. more than 1 vertex at a time.

## Minor changes

## Bug fixes

## Deprecated and defunct

## Documentation

# QDECR 0.7.0: OHBM

Version 0.7.0 is the first version that is publically released. It is also the version that was presented at OHBM 2019 and was thus named "OHBM". As we do not have any formal news for versions before 0.7.0, we will keep it brief.

## New functions and features

* Creation of the QDECR package (before 0.7.0 it was a series of associated scripts)
* Creation of the vw object and associated functions
* Handling of imputed datasets (via imp2list)
* Creation of all summary and plot functions
* Creation of the `stack` concept
