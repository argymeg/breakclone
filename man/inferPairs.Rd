% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readFunctions.R
\name{inferPairs}
\alias{inferPairs}
\title{Infer tumour pairs from segment table}
\usage{
inferPairs(segmentTable, sep = "_")
}
\arguments{
\item{segmentTable}{A segment or mutation table generated by the appropriate read function}

\item{sep}{The separator between patient and tumour identifiers in the sample names}
}
\value{
A data.frame with column names Sample1 and Sample2
}
\description{
If sample IDs comprise a patient ID and a tumour type identifier, this function can attempt to infer the tumour pairs from the cohort.
}
