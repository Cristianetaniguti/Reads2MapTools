#' Write VCF files from onemap objects
#'
#' Only available for biallelic markers. OneMap don't keep the reference and alternative 
#' allele information. Do not consider the outputted information of the alleles for further 
#' investigations. We can't properly relate homozygous progeny from both parents heterozygous 
#' with the reference (0/0) or alternative (1/1) alleles. 
#' 
#' @param onemap.object object of onemap class
#' @param out_vcf path and name of the vcf file to be outputed from the generated onemap object. 
#' It is important to notice that only onemap informative markers are kept.
#' @param input_info_rds output file from onemap_read_vcfR to store alleles information
#' @param probs genotyope probabilities
#' @param parent1 parent 1 identification 
#' @param parent2 parent 2 identification
#' 
#' @return vcf file
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br} 
#' 
#' @import vcfR
#'  
#' @export
onemap_write_vcfR <- function(onemap.object, 
                              out_vcf, 
                              input_info_rds,
                              parent1 = NULL,
                              parent2 = NULL,
                              probs, cores=1){
  
  vcf.template <- read.vcfR(system.file("extdata/vcf_example_out.vcf.gz", package = "onemap"), verbose = F)
  
  info <- readRDS(input_info_rds)
  info$POS <- as.character(as.numeric(info$POS))
  vcf.template@fix <- as.matrix(cbind(info[,-c(6:9)], QUAL = NA, FILTER = "PASS", INFO = NA))
  
  # GT
  GT <- t(onemap.object$geno)
  idx.m <- match(colnames(onemap.object$geno), vcf.template@fix[,3])
  info <- info[idx.m,]
  
  GT[which(onemap.object$geno == 0)] <- "./."
  GT[which(onemap.object$geno == 2)] <- "0/1"
  
  # B3.7
  idx <- which(info[,6] == 0 & info[,7] == 1 & info[,8] == 0 & info[,9] == 1)
  GT[idx, ][which(GT[idx, ] == 1)] <- "0/0" # Here we can't be sure
  GT[idx, ][which(GT[idx, ] == 3)] <- "1/1" # Here we can't be sure
  
  # D1.10
  idx <- which(info[,6] == 0 & info[,7] == 1 & info[,8] == 0 & info[,9] == 0)
  GT[idx, ][which(GT[idx, ] == 1)] <- "0/0"
  idx <- which(info[,6] == 0 & info[,7] == 1 & info[,8] == 1 & info[,9] == 1)
  GT[idx, ][which(GT[idx, ] == 1)] <- "1/1" 
  
  # D2.15
  idx <- which(info[,6] == 0 & info[,7] == 0 & info[,8] == 0 & info[,9] == 1)
  GT[idx, ][which(GT[idx, ] == 1)] <- "0/0" 
  idx <- which(info[,6] == 1 & info[,7] == 1 & info[,8] == 0 & info[,9] == 1)
  GT[idx, ][which(GT[idx, ] == 1)] <- "1/1" 
  
  # PL
  PL <- probs
  PL <- round(-10*log(PL, base = 10),0)
  PL[which(PL == "Inf")] <- 99
  
  GQ <- apply(PL, 1, function(x) {
    temp <- x[-which.min(x)]
    temp[which.min(temp)]
  })
  
  PL <- apply(PL, 1, function(x) paste(x, collapse = ","))
  PL <- paste0(GQ, ":",PL)
  PL <- split(PL, rep(1:onemap.object$n.mar, each = onemap.object$n.ind))
  PL <- do.call(rbind, PL)
  
  gt <- matrix(paste0(GT, ":", PL), nrow = dim(GT)[1])
  colnames(gt) <- rownames(onemap.object$geno)
  
  ## Parents - parents are deterministic in maps, then we output here a very low error probability
  P1 <- apply(info[,6:7], 1, function(x) paste0(x, collapse = "/"))
  P2 <- apply(info[,8:9], 1, function(x) paste0(x, collapse = "/"))
  parents <- cbind(P1, P2)
  
  parents[which(parents == "1/1")] <- "1/1:99:99,99,0"
  parents[which(parents == "0/0")] <- "0/0:99:0,99,99"
  parents[which(parents == "0/1")] <- "0/1:99:99,0,99"
  
  gt <- cbind(parents,gt)
  colnames(gt)[1:2] <- c(parent1, parent2)
  gt <- cbind(FORMAT="GT:GQ:PL", gt)
  
  vcf.template@gt <- gt
  vcf.template@meta[2] <- "##source=OneMap"
  vcf.template@meta <- vcf.template@meta[-c(4:5,8)]
  
  write.vcf(x = vcf.template, file = out_vcf)
}
