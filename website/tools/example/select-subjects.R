#!/usr/bin/env Rscript
# Picks the example's subjects from ABIDE's phenotype file and writes the data frame the
# analysis reads, phenotypes.csv, with the columns id, age, sex and site:
#
#   Rscript select-subjects.R <Phenotypic_V1_0b_preprocessed1.csv> <phenotypes.csv> <excluded.csv>
#
# The rule: the controls of one site, modelled on age and sex only. NYU
# is the site, because it has the most controls (100) over the widest span of ages (6 to
# 32). Subjects whose anatomical scan failed the Preprocessed Connectomes Project's
# quality check are left out and listed in excluded.csv, so the site can say so. Nothing
# else is filtered: no age or sex balancing, no handedness, no IQ.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: select-subjects.R <abide phenotypes> <phenotypes.csv> <excluded.csv>")

pheno <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)

# DX_GROUP 2 is a control (1 is autism). FILE_ID is the subject's directory in the
# FreeSurfer output, "no_filename" for subjects that were not preprocessed.
controls <- pheno[pheno$SITE_ID == "NYU" & pheno$DX_GROUP == 2 & pheno$FILE_ID != "no_filename", ]

# The PCP's raters marked each anatomical scan OK, maybe or fail; a blank means no
# rating. Only a fail excludes. (%in% rather than ==, so a blank is not NA.)
failed <- controls$qc_anat_rater_2 %in% "fail" | controls$qc_anat_rater_3 %in% "fail"

kept <- controls[!failed, ]
subjects <- data.frame(
  id = kept$FILE_ID,
  age = kept$AGE_AT_SCAN,
  # SEX is 1 for male and 2 for female. Words, so that the analysis can name the
  # reference level rather than remember a code.
  sex = ifelse(kept$SEX == 1, "male", "female"),
  site = kept$SITE_ID,
  stringsAsFactors = FALSE
)
subjects <- subjects[order(subjects$id), ]
write.csv(subjects, args[2], row.names = FALSE, quote = FALSE)

excluded <- data.frame(id = controls$FILE_ID[failed], reason = "anatomical scan failed the PCP quality check")
write.csv(excluded, args[3], row.names = FALSE, quote = FALSE)

cat(sprintf(
  "%d NYU controls with FreeSurfer output; %d excluded (anatomical QC failed); %d kept: %d female, %d male, aged %.1f to %.1f.\n",
  nrow(controls), sum(failed), nrow(subjects), sum(subjects$sex == "female"), sum(subjects$sex == "male"),
  min(subjects$age), max(subjects$age)
))
