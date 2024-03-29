% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_multiallelic_vcf.R
\name{filter_geno_multi}
\alias{filter_geno_multi}
\title{Filter by parents alleles. Every genotype that has alleles absent 
in the parents will be considered missing data (NA).}
\usage{
filter_geno_multi(gt.multi, P1, P2)
}
\arguments{
\item{gt.multi}{genotypes matrix (markers x individuals) with multiallelic markers}

\item{P1}{character with one of the parents ID}

\item{P2}{character with other parent ID}
}
\description{
Filter by parents alleles. Every genotype that has alleles absent 
in the parents will be considered missing data (NA).
}
