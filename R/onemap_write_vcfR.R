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
#' @param parent1.geno parent 1 genotypes (vcf notation: 0/0, 0/1, 1/1, ./.)
#' @param parent2.geno parent 2 genotypes (vcf notation: 0/0, 0/1, 1/1, ./.)
#' @param parent1.id parent 1 ID
#' @param parent2.id parent 2 ID
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
                              parent1.geno = NULL,
                              parent1.id = NULL,
                              parent2.geno = NULL,
                              parent2.id = NULL,
                              probs, 
                              ad_matrix = NULL){
  
  vcf.template <- read.vcfR(system.file("extdata/vcf_example_out.vcf.gz", package = "onemap"), verbose = F)
  
  info <- readRDS(input_info_rds)
  info$POS <- as.character(as.numeric(info$POS))
  
  # GT
  GT <- t(onemap.object$geno)
  idx.m <- match(colnames(onemap.object$geno), info$ID)
  info <- info[idx.m,]
  vcf.template@fix <- as.matrix(cbind(as.data.frame(info), QUAL = NA, FILTER = "PASS", INFO = NA))
  
  GT[which(GT == 0)] <- "./."
  GT[which(GT == 2)] <- "0/1"
  
  parent1.geno <- parent1.geno[which(names(parent1.geno) %in% rownames(GT))]
  parent2.geno <- parent2.geno[which(names(parent2.geno) %in% rownames(GT))]
  
  # B3.7
  idx <- which(parent1.geno == "0/1" & parent2.geno == "0/1")
  GT[idx, ][which(GT[idx, ] == "1")] <- "0/0" # Here we can't be sure
  GT[idx, ][which(GT[idx, ] == "3")] <- "1/1" # Here we can't be sure
  
  # D1.10
  idx <- which(parent1.geno == "0/1"& parent2.geno == "0/0")
  GT[idx, ][which(GT[idx, ] == "1")] <- "0/0"
  GT[idx, ][which(GT[idx, ] == "3")] <- "./."
  idx <- which(parent1.geno == "0/1"  & parent2.geno == "1/1")
  GT[idx, ][which(GT[idx, ] == "1")] <- "1/1" 
  GT[idx, ][which(GT[idx, ] == "3")] <- "./." 
  
  # D2.15
  idx <- which(parent1.geno == "0/0" & parent2.geno == "0/1")
  GT[idx, ][which(GT[idx, ] == "1")] <- "0/0" 
  GT[idx, ][which(GT[idx, ] == "3")] <- "./."
  idx <- which(parent1.geno == "1/1" & parent2.geno == "0/1")
  GT[idx, ][which(GT[idx, ] == "1")] <- "1/1" 
  GT[idx, ][which(GT[idx, ] == "3")] <- "./." 
  
  # non-informative
  idx <- which((parent1.geno == "0/0" & parent2.geno == "0/0") | 
                 (parent1.geno == "1/1" & parent2.geno == "1/1") |
                 (parent1.geno == "0/0" & parent2.geno == "1/1") |
                 (parent1.geno == "1/1" & parent2.geno == "0/0") |
                 any(parent1.geno == "./." | parent2.geno == "./.")
  )
  GT[idx, ][which(GT[idx, ] == "1")] <- "0/0" 
  GT[idx, ][which(GT[idx, ] == "3")] <- "1/1"
  
  # PL
  PL <- probs
  PL <- -10*log(PL, base = 10)
  PL <- apply(PL, 2, function(x) {
    x[which(x == "Inf" | x == "-Inf" | x > 99)] <- 99
    return(x)
  })
  
  if(length(which(is.na(PL))) > 0)
    PL[which(is.na(PL))] <- 0
  PL <- t(apply(PL, 1, function(x) x-x[which.min(x)]))
  PL <- floor(PL)
  
  GQ <- apply(PL, 1, function(x) {
    temp <- x[-which.min(x)]
    temp[which.min(temp)]
  })
  
  PL[is.na(probs)] <- "."
  GQ[which(GQ==0)] <- "."
  
  PL <- apply(PL, 1, function(x) paste(rev(x), collapse = ","))
  PL <- paste0(GQ, ":",PL)
  PL <- split(PL, rep(1:onemap.object$n.mar, each = onemap.object$n.ind))
  PL <- do.call(rbind, PL)
  
  gt <- matrix(paste0(GT, ":", PL), nrow = dim(GT)[1])
  colnames(gt) <- rownames(onemap.object$geno)
  
  ## Parents - parents are deterministic in maps, then we output here a very low error probability
  parents <- cbind(parent1.geno, parent2.geno)
  
  parents[which(parents == "1/1")] <- "1/1:99:99,99,0"
  parents[which(parents == "0/0")] <- "0/0:99:0,99,99"
  parents[which(parents == "0/1")] <- "0/1:99:99,0,99"
  
  gt <- cbind(parents,gt)
  colnames(gt)[1:2] <- c(parent1.id, parent2.id)
  
  FORMAT <- "GT:GQ:PL"
  
  if(!is.null(ad_matrix)){
    idx.c <- match(colnames(gt), colnames(ad_matrix))
    idx.r <- match(rownames(gt), rownames(ad_matrix))
    
    ad_matrix <- ad_matrix[idx.r, idx.c]
    gt1 <- matrix(paste0(gt, ":", ad_matrix), nrow = nrow(gt))
    colnames(gt1) <- colnames(gt)
    rownames(gt1) <- rownames(gt)
    gt <- gt1
    FORMAT <- "GT:GQ:PL:AD"
  }
  
  gt <- cbind(FORMAT, gt)
  
  vcf.template@gt <- gt
  vcf.template@meta[2] <- "##source=OneMap"
  vcf.template@meta <- vcf.template@meta[-c(5,8)]
  
  write.vcf(x = vcf.template, file = out_vcf)
}
