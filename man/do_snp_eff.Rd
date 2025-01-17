% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{do_snp_eff}
\alias{do_snp_eff}
\title{runs SNPeff with a candiSNP input file and a preinstalled genome}
\usage{
do_snp_eff(file, genome)
}
\arguments{
\item{file}{path to .csv file (Chr, Pos, Ref, Alt, Allele_freq)}

\item{genome}{genome name (annotation key) for annotation to use}
}
\value{
dataframe
}
\description{
runs SNPeff with a candiSNP input file and a preinstalled genome
}
\examples{
\dontrun{ #cant work without installing a large db
 do_snp_eff("athal_sample.csv", "athalianaTair10")
 }
}
