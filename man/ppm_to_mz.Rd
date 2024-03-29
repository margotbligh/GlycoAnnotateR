% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ppm_to_mz}
\alias{ppm_to_mz}
\title{Convert from ppm to mz}
\usage{
ppm_to_mz(mz, ppm)
}
\description{
Simple function to calculate from ppm to mz, for a given mz and
amount of ppm.
}
\section{Slots}{

\describe{
\item{\code{mz}}{\emph{m/z} value to calculate for}

\item{\code{ppm}}{the ppm amount to calculate}
}}

\examples{

#calculate 2 ppm of m/z value of 259.0129

ppm_to_mz(259.0129, 2) 

}
