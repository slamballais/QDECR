% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarizing.R
\name{summary.vw_fastlm}
\alias{summary.vw_fastlm}
\title{Summarize qdecr_fastlm clusters}
\usage{
\method{summary}{vw_fastlm}(
  object,
  verbose = FALSE,
  annot = FALSE,
  file = "aparc.annot",
  regions = 3,
  ...
)
}
\arguments{
\item{object}{The output object of a qdecr call (e.g. qdecr_fastlm)}

\item{verbose}{Logical; if TRUE, it outputs some information of the steps; default is FALSE}

\item{annot}{Logical; if TRUE, `file` will be read in to obtain information on the top regions in that cluster}

\item{file}{string; the name of the file within the target template's label directory if annot = TRUE}

\item{regions}{integer; number of regions that should be added if annot = TRUE}

\item{...}{Further arguments for `summary`}
}
\description{
Summarizes the significant clusters obtained from qdecr_fastlm
}
\details{
This function returns an overview of the significant clusters from qdecr_fastlm,
and provides information per cluster. Additionally, it is possible to obtain annotated
information, given that it is located in the target template's label directory.
}
