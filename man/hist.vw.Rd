% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{hist.vw}
\alias{hist.vw}
\title{Histograms of the vertex measures}
\usage{
\method{hist}{vw}(x, qtype = c("vertex", "subject"), xlab = NULL, main = NULL, ...)
}
\arguments{
\item{x}{The output object of a qdecr call (e.g. qdecr_fastlm)}

\item{qtype}{A string; either "vertex" for vertex-wise means, or "subject" for subject-wise means}

\item{xlab}{See `?hist`}

\item{main}{see `?hist`}

\item{...}{Further arguments for `hist`}
}
\value{
see `?hist`
}
\description{
Plots a histogram of the mean vertex measure, either by vertex or by subject
}
\details{
This function works on top of the base `hist` function. It simply makes a histogram
of the mean vertex values either per vertex or per subject (if you specify `qtype = "subject"`)
}
