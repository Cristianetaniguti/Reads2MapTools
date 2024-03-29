% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/updog_genotype.R
\name{updog_genotype}
\alias{updog_genotype}
\title{OneMap interface with updog package}
\usage{
updog_genotype(
  vcf = NULL,
  vcf.par = c("AD", "DPR"),
  out_vcf = NULL,
  parent1 = "P1",
  parent2 = "P2",
  f1 = NULL,
  crosstype = NULL,
  recovering = FALSE,
  cores = 2,
  depths = NULL,
  global_error = NULL,
  use_genotypes_errors = TRUE,
  use_genotypes_probs = FALSE,
  rm_multiallelic = TRUE
)
}
\arguments{
\item{vcf}{path and name of the vcf file}

\item{vcf.par}{Field of VCF that informs the depth of alleles}

\item{out_vcf}{path and name of the vcf file to be outputed from the generated onemap object. 
It is important to notice that only onemap informative markers are kept.}

\item{parent1}{parent 1 identification in vcfR object}

\item{parent2}{parent 2 identification in vcfR objetc}

\item{f1}{f1 individual identification if f2 cross type}

\item{crosstype}{string defining the cross type, by now it supports only 
outcross and f2 intercross}

\item{recovering}{logical defining if markers should be recovered from VCF (TRUE) or just keep the ones in the onemap object (FALSE).}

\item{cores}{number of threads}

\item{depths}{list containing a matrix for ref and other for alt allele counts, samples ID in colum and markers ID in rows}

\item{global_error}{number from 0 to 1 defining the global error to be considered together 
with the genotype errors or the genotype probabilities or NULL to not considered any global error}

\item{use_genotypes_errors}{if \code{TRUE} the error probability of each genotype will be considered in emission function of HMM}

\item{rm_multiallelic}{if \code{TRUE} multiallelic markers will be removed from the output onemap object}

\item{use_genotype_probs}{if \code{TRUE} the probability of each possible genotype will be considered in emission function of HMM}

\item{output_info_file}{define a name for the file with alleles information.}
}
\value{
onemap object with genotypes updated
}
\description{
Uses alelle counts to reestimate genotypes with updog approach and 
stores the genotypes probabilities or for further multipoint 
analysis
}
\references{
Gerard, D., Ferrão L.F.V., Garcia, A.A.F., & Stephens, M. (2018). Harnessing 
Empirical Bayes and Mendelian Segregation for Genotyping Autopolyploids from 
Messy Sequencing Data. bioRxiv. doi: 10.1101/281550.
}
\seealso{
\code{\link[onemap]{extract_depth}} 
    \code{\link[onemap]{binom_genotype}} and 
    \url{https://github.com/dcgerard/updog}.
}
\author{
Cristiane Taniguti, \email{chtaniguti@usp.br}
}
