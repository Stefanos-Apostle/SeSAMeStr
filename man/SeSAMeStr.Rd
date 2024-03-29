% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeSAMeStr_utils.R
\name{SeSAMeStr}
\alias{SeSAMeStr}
\title{Final Assembly of SeSAMe STREET Subassemblies/Utilities}
\usage{
SeSAMeStr(
  Idat_dir,
  out_dir,
  sample_sheet,
  prep,
  formula,
  cores,
  subsample = NA
)
}
\arguments{
\item{Idat_dir}{Path to directory containing raw Idat files}

\item{out_dir}{Path to output directory}

\item{sample_sheet}{Path to SeSAMe STREET sample sheet}

\item{prep}{Flag code for SeSAMe Preprocessing}

\item{formula}{Formula for DML calculation}

\item{cores}{Number of cores to run DML summary stats calculation}

\item{subsample}{Number of beta values to use for DML calc and analysis for test runs (leave blank for full scale run)}
}
\value{
All outputs in or a sub directory of the out_dir
}
\description{
Function to run the SeSAMe Street Pipeline
}
\examples{
SeSAMeStr("path/Idat_dir", "path/out_dir", "path/sample_sheet.xlsx", "TQCDPB", ~ Group, 16, NA)
}
