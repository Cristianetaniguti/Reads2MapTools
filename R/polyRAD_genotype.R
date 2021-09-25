#' OneMap interface with polyRAD package
#'
#' Uses alelle counts to reestimate genotypes with updog approach and 
#' stores the genotypes probabilities or for further multipoint 
#' analysis
#' 
#' @param vcf path and name of the vcf file
#' @param out_vcf path and name of the vcf file to be outputed from the generated onemap object. 
#' It is important to notice that only onemap informative markers are kept.
#' @param parent1 parent 1 identification in vcfR object
#' @param parent2 parent 2 identification in vcfR objetc
#' @param f1 f1 individual identification if F2 cross type
#' @param crosstype string defining the cross type, by now it supports only 
#' outcross and f2 intercross
#' @param global_error number from 0 to 1 defining the global error to be considered together 
#' with the genotype errors or the genotype probabilities or NULL to not considered any global error
#' @param use_genotypes_errors if \code{TRUE} the error probability of each genotype will be considered in emission function of HMM
#' @param use_genotype_probs if \code{TRUE} the probability of each possible genotype will be considered in emission function of HMM
#' @param rm_multiallelic if \code{TRUE} multiallelic markers will be removed from the output onemap object
#' 
#' 
#' @return onemap object with genotypes updated
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br} 
#' @seealso \code{\link[onemap]{extract_depth}} 
#'     \code{\link[onemap]{binom_error}} and 
#'     \url{https://github.com/dcgerard/updog}.
#'
#' @references 
#'
#' Clark LV, Lipka AE, and Sacks EJ (2019) Improvements to 
#' Genotype Calling in Polyploids Using the polyRAD R Package. 
#' Plant and Animal Genome Conference XXVII, January 12-16, 
#' San Diego, California, USA. doi:10.13140/RG.2.2.18358.75847
#' 
#' Clark LV, Lipka AE, and Sacks EJ (2018) polyRAD: Genotype 
#' Calling with Uncertainty from Sequencing Data in Polyploids 
#' and Diploids. Plant and Animal Genome Conference XXVI, 
#' January 13-17, San Diego, California, USA. doi:10.13140/RG.2.2.27134.08001
#'
#' @import polyRAD
#' @importFrom vcfR read.vcfR 
#' @import onemap 
#'   
#' @export
polyRAD_genotype <- function(vcf=NULL, 
                             out_vcf = NULL,
                             parent1=NULL,
                             parent2=NULL,
                             f1=NULL,
                             crosstype=NULL,
                             global_error = NULL,
                             use_genotypes_errors = TRUE,
                             use_genotypes_probs = FALSE,
                             rm_multiallelic = TRUE,
                             info_file_name = NULL,
                             recovering = FALSE){
  
  vcfR.object <- read.vcfR(vcf, verbose = F)
  
  info_file_name <- tempfile() 
  onemap.obj <- onemap_read_vcfR(vcfR.object = vcfR.object,
                                 parent1=parent1,
                                 parent2=parent2,
                                 f1=f1,
                                 cross=crosstype, 
                                 only_biallelic = F)
  
  MKS <- vcfR.object@fix[,3]
  if (any(MKS == "." | is.na(MKS))) {
    MKS <- paste0(vcfR.object@fix[,1],"_", vcfR.object@fix[,2])
    # Add tag if is duplicated positions (split form of mnps)
    for(i in 2:length(MKS)) {
      if(MKS[i] == paste0(strsplit(MKS[i-1], "_")[[1]][1:2], collapse = "_")) {
        z <- z + 1
        MKS[i] <- paste0(MKS[i], "_",z)
      } else {
        z <- 0
      }
    }
  }
  
  info <- data.frame(CHROM = vcfR.object@fix[,1], 
                     POS = as.numeric(vcfR.object@fix[,2]), 
                     ID = MKS, 
                     REF = vcfR.object@fix[,4], 
                     ALT = vcfR.object@fix[,5])
  
  saveRDS(info, file = info_file_name)
  
  # Do the checks
  poly.test <- VCF2RADdata(vcf, phaseSNPs = FALSE, 
                           min.ind.with.reads = 0,
                           min.ind.with.minor.allele = 0)
  if(crosstype=="f2 intercross"){
    poly.test <- SetDonorParent(poly.test, f1)
    poly.test <- SetRecurrentParent(poly.test, f1)
  } else if(crosstype=="outcross"){
    poly.test <- SetDonorParent(poly.test, parent1)
    poly.test <- SetRecurrentParent(poly.test, parent2)
  }
  
  mydata2 <- PipelineMapping2Parents(poly.test, 
                                     freqAllowedDeviation = 0.06,
                                     useLinkage = FALSE)
  
  
  
  seed.uniq <- sample(100000, 1)
  Export_MAPpoly(mydata2, paste0("temp.file.", seed.uniq)) 
  
  genotypes <- read.table(paste0("temp.file.",seed.uniq), skip=12)
  
  file.remove(paste0("temp.file.",seed.uniq))
  
  # this will change according to the vcf - bug!! Need attention!
  if(any(grepl(":", as.character(genotypes$V1)))){
    temp_list <- strsplit(as.character(genotypes$V1), split = "_")
    temp <- sapply(temp_list, function(x) paste0(x[-c(length(x)-1,length(x))], collapse = "_"))
    pos <- gsub(":", "_", temp)
  } else {
    temp_list <- strsplit(as.character(genotypes$V1), split = "_")
    pos <- sapply(temp_list, function(x) if(length(x) > 2) paste0(x[1:2], collapse = "_") else x[1])
  }
  
  multi <- names(which(table(pos) > onemap.obj$n.ind))
  if(length(multi) != 0){
    ## From vcf
    multi.idx <- which(pos %in% multi)
    genotypes <- droplevels(genotypes[-multi.idx,])
    pos <- pos[-multi.idx]
    
    # removing the multiallelics from the onemap object
    if(crosstype == "outcross"){
      idx.multi <- which(!(onemap.obj$segr.type %in% c("B3.7", "D1.10", "D2.15")))
      mult.obj <- split_onemap(onemap.obj, mks= idx.multi)
      idx.bi <- which(onemap.obj$segr.type %in% c("B3.7", "D1.10", "D2.15"))
      onemap.obj <- split_onemap(onemap.obj, mks= idx.bi)
    } else if(crosstype == "f2 intercross"){
      idx.multi <- which(!(onemap.obj$segr.type %in% c("A.H.B")))
      mult.obj <- split_onemap(onemap.obj, mks= idx.multi)
      idx.bi <- which(onemap.obj$segr.type %in% c("A.H.B"))
      onemap.obj <- split_onemap(onemap.obj, mks= idx.bi)
    }
  }
  
  parent1.geno <-t(mydata2$likelyGeno_donor)
  
  mks <- sapply(strsplit(rownames(parent1.geno), "_"), function(x) {
    paste0(x[-c(length(x)-1, length(x))], collapse = "_")
  })
  mks <- gsub(":", "_", mks)
  geno <- split(parent1.geno, mks)
  P1 <- sapply(geno, function(x){
    if(anyNA(x)){
      return("./.")
    }else if(all(x == c(1,1))){
      return("0/1")
    } else if(all(x == c(2,0))){
      return("0/0")
    } else if(all(x == c(0,2))){
      return("1/1")
    }
  })
  names(P1) <- unique(mks)
  
  parent2.geno <- t(mydata2$likelyGeno_recurrent)
  mks <- sapply(strsplit(rownames(parent2.geno), "_"), function(x) {
    paste0(x[-c(length(x)-1, length(x))], collapse = "_")
  })
  mks <- gsub(":", "_", mks)
  geno <- split(parent2.geno, mks)
  P2 <- sapply(geno, function(x){
    if(anyNA(x)){
      return("./.")
    }else if(all(x == c(1,1))){
      return("0/1")
    } else if(all(x == c(2,0))){
      return("0/0")
    } else if(all(x == c(0,2))){
      return("1/1")
    }
  })
  names(P2) <- unique(mks)

  unq.pos <- unique(pos)
  P1 <- P1[which(names(P1) %in% unq.pos)]
  P2 <- P2[which(names(P2) %in% unq.pos)]

  if(!recovering){
    pos.onemap <- colnames(onemap.obj$geno)
    genotypes <- genotypes[which(pos%in%pos.onemap),]
    pos <- pos[which(pos%in%pos.onemap)]
    P1 <- P1[which(names(P1) %in% pos.onemap)]
    P2 <- P2[which(names(P2) %in% pos.onemap)]
  }
  
  # Remove parents
  if(crosstype=="f2 intercross"){
    genotypes <- genotypes[-which(genotypes[,2]%in%parent1),]
    genotypes <- genotypes[-which(genotypes[,2] %in%parent2),]
  }
  
  # Updating geno matrix
  new.geno <- apply(genotypes[,3:5], 1, which.max)
  maxpostprob <- apply(genotypes[,3:5], 1, function(x) x[which.max(x)])
  
  new.geno <- matrix(new.geno,nrow = onemap.obj$n.ind, ncol = length(unique(pos)))
  maxpostprob <- matrix(maxpostprob,nrow = onemap.obj$n.ind, ncol = length(unique(pos)))
  colnames(new.geno) <- colnames(maxpostprob) <- unique(pos)
  rownames(new.geno) <- rownames(maxpostprob) <- rownames(onemap.obj$geno)
  
  # Same order priority: individuals, markers
  genotypes$V2 <- factor(genotypes$V2, levels = genotypes$V2[1:onemap.obj$n.ind])
  
  genotypes <- genotypes[order(genotypes$V2),]
  
  # marker type
  mk.type <- rep(NA, length(P1))
  mk.type.num <- rep(NA, length(P1))
  
  idx <- which(P1 == "0/1" & P2 == "0/1")
  mk.type[idx] <- "B3.7"
  mk.type.num[idx] <- 4
  idx <- which((P1 == "0/1" & P2 == "0/0") | (P1 == "0/1" & P2 == "1/1"))
  mk.type[idx] <- "D1.10"
  mk.type.num[idx] <- 6
  idx <- which((P1 == "0/0" & P2 == "0/1") | (P1 == "1/1" & P2 == "0/1"))
  mk.type[idx] <- "D2.15"
  mk.type.num[idx] <- 7
  keep.mks <- which(!is.na(mk.type.num))
  
  # Removing markers
  onemap.obj$geno <- new.geno[,keep.mks]
  onemap.obj$n.mar <- dim(onemap.obj$geno)[2]
  onemap.obj$segr.type <- mk.type[keep.mks]
  onemap.obj$segr.type.num <- mk.type.num[keep.mks]
  onemap.obj$CHROM <- info$CHROM[info$ID %in% colnames(onemap.obj$geno)]
  onemap.obj$POS <- info$POS[info$ID %in% colnames(onemap.obj$geno)]

  # Avoiding geno 3 in D1.10 and D2.15
  idx <- which(onemap.obj$segr.type == "D1.10")
  onemap.obj$geno[,idx][onemap.obj$geno[,idx] == 3] <- 1
  idx <- which(onemap.obj$segr.type == "D2.15")
  onemap.obj$geno[,idx][onemap.obj$geno[,idx] == 3] <- 1
  
  probs <- as.matrix(genotypes[,3:5])
  probs <- probs[keep.mks + rep(c(0:(onemap.obj$n.ind-1))*onemap.obj$n.mar, each=length(keep.mks)),]
  
  maxpostprob <- maxpostprob[,keep.mks]
  colnames(onemap.obj$geno) <- unq.pos[keep.mks]
  
  if(use_genotypes_probs){
    onemap.obj.new <- create_probs(input.obj = onemap.obj,
                                   genotypes_probs = probs,
                                   global_error = global_error)
  } else if(use_genotypes_errors){
    onemap.obj.new <- create_probs(input.obj = onemap.obj,
                                   genotypes_errors = 1- maxpostprob,
                                   global_error = global_error)
  } else if(!is.null(global_error)){
    onemap.obj.new <- create_probs(input.obj = onemap.obj,
                                   global_error = global_error)
  }
  
  if(!rm_multiallelic){
    if(length(multi) > 0)
      onemap.obj.new <- combine_onemap(onemap.obj.new, mult.obj)
  }
  
  if(!is.null(out_vcf)){
    onemap_write_vcfR(onemap.object = onemap.obj.new, 
                      out_vcf = out_vcf, 
                      input_info_rds = info_file_name,
                      probs = probs, 
                      parent1.id = parent1, 
                      parent2.id = parent2, 
                      parent1.geno = P1,
                      parent2.geno = P2)
  }
  
  return(onemap.obj.new)
}
