% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeSAMeStr_utils.R
\name{testEnrichment_DMLs}
\alias{testEnrichment_DMLs}
\title{Plot Test Enrichments}
\usage{
testEnrichment_DMLs(test_result, CONDITION, LEVEL, out_dir)
}
\arguments{
\item{test_result}{Test results of DML summary stats}

\item{CONDITION}{Metavalue realted to colname of sample sheet}

\item{LEVEL}{Which comparison of condition levels to use DML results from}

\item{out_dir}{Path to output directory}
}
\value{
Pdf of test enrichment results for CONDITION-LEVEL DML CpGs
}
\description{
Function to calculate enrichment for different SeSAMe CpG databases from DML CpGs
}
\examples{
testEnrichment_DMLs(test_result, "Group", "Treatment", "path/out_dir")
}
