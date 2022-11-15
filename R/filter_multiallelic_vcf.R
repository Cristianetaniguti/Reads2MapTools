#' Filter vcf file keeping only alleles found in parents haplotypes
#' 
#' @param vcf path to vcf file
#' @param P1 character with one of the parents ID
#' @param P2 character with other parent ID
#' @param max.missing maximum fraction of missing data allowed by marker
#' @param vcf.out character defining output file name
#' @param ploidy integer defining the samples ploidy
#' @param min.gpm set minimum value for GPM parameter. Genotypes bellow this value will be replaced by missing data.
#' @param min.phpm set minimum value for PHPM parameter. Genotypes bellow this value will be replaced by missing data.
#' 
#' @import vcfR
#' 
#' @export
filter_multi_vcf <- function(vcf.file, P1, P2, ploidy, max.missing = NULL, min.gpm = NULL,
                             min.phpm = NULL, vcf.out = "filtered.vcf.gz"){
  vcf <- read.vcfR(vcf.file, verbose = FALSE)
  
  gt <- extract.gt(vcf)
  gq <- extract.gt(vcf, element = "GQ")
  
  filt.gt <- filter_geno_multi(gt.multi = gt, P1, P2)
  fix.up <- vcf@fix
  idx <- match(rownames(filt.gt), paste0(fix.up[,1],"_",fix.up[,2]))
  old.gt <- gt[idx,]
  fix.up <- fix.up[idx,]
  gq <- gq[idx,]

  up.fix <- get_alternatives(fix = fix.up, old.gt = old.gt, P1, P2)
  
  if(!is.null(min.gpm)){
    gpm <- extract.gt(vcf, element = "GPM")
    gpm <- gpm[idx,]
    gpm <- apply(gpm, 2, as.numeric)
    filt.gt[which(gpm < min.gpm)] <- NA
  }
  
  if(!is.null(min.phpm)){
    phpm <- extract.gt(vcf, element = "PHPM")
    phpm <- phpm[idx,]
    phpm <- apply(phpm, 2, as.numeric)
    filt.gt[which(phpm < min.phpm)] <- NA
  }
  
  if(!is.null(max.missing)){
    mis <- apply(filt.gt, 1, function(x) sum(is.na(x))/length(x))
    idx <- which(mis > max.missing)
    
    filt.gt <- filt.gt[-idx,]
    gq<- gq[-idx,]
    phpm <- phpm[-idx,]
    gpm <- gpm[-idx,]
    up.fix <- up.fix[-idx,]
  }
  
  filt.gt[is.na(filt.gt)] <- paste0(rep(".", ploidy), collapse = "/")
  format <- matrix(paste0(filt.gt, ":", gq), nrow = nrow(filt.gt))
  
  colnames(format) <- colnames(filt.gt)
  vcf.new <- vcf
  vcf.new@fix <- up.fix
  vcf.new@gt <- cbind(FORMAT="GT:GQ", format)
  
  write.vcf(vcf.new, file = vcf.out)
}

#' Filter by parents alleles. Every genotype that has alleles absent 
#' in the parents will be considered missing data (NA).
#' 
#' @param gt.multi genotypes matrix (markers x individuals) with multiallelic markers 
#' @param P1 character with one of the parents ID
#' @param P2 character with other parent ID
#' 
filter_geno_multi <- function(gt.multi, P1, P2){
  parents.id <- which(colnames(gt.multi) %in% c(P1, P2))
  
  diff.geno <- apply(gt.multi, 1, function(x) length(unique(x)))
  gt.multi <- gt.multi[-which(diff.geno == 1),]
  
  mk.split <- apply(gt.multi, 1, function(x) strsplit(x, "/"))
  
  # Filter non-informative
  parents1 <- lapply(mk.split, "[[", parents.id[1])
  parents2 <- lapply(mk.split, "[[", parents.id[2])
  parents1 <- lapply(parents1, function(x) length(unique(x)))
  parents2 <- lapply(parents2, function(x) length(unique(x)))
  parents1 <- do.call(rbind, parents1)
  parents2 <- do.call(rbind, parents2)
  rm.mk <- which(parents1 == 1 & parents2 == 1)
  mk.split <- mk.split[-rm.mk]
  
  mk.filt <- lapply(mk.split, function(y) {
    
    paren <- y[parents.id]
    
    # All genotypes with "." will be considered missing data
    idx <- lapply(paren, function(x) which(x == "."))
    if(length(idx[[1]]) > 0) paren[[1]] <- NA
    if(length(idx[[2]]) > 0)  paren[[2]] <- NA
    
    # The parents alleles will be re codified as 0 to 4 
    un.p <- unique(unlist(paren))
    un.p <- sort(un.p)
    un.p.up <- 0:(length(un.p) - 1)
    
    progeny <- sapply(y[-parents.id], function(x) {
      if(all(is.na(x)) | !all(x %in% un.p)) {
        return(NA)
      } else {
        res <- un.p.up[match(x, un.p)]
        paste0(res, collapse = "/")
      }
    })
    
    # Avoiding NA as character
    p1 <- if(!any(is.na(paren[[1]]))) paste0(un.p.up[match(paren[[1]], un.p)], collapse = "/") else NA
    p2 <- if(!any(is.na(paren[[2]]))) paste0(un.p.up[match(paren[[2]], un.p)], collapse = "/") else NA
    parents <- c(p1, p2)
    names(parents) <- c(P1, P2)
    
    return(c(progeny, parents))
  })
  
  new.gt <- do.call("rbind", mk.filt)
  return(new.gt)
}

#' Replace the reference and alternative alleles in fix to 
#' keep only the ones found in the parents
#' 
#' @param fix fix part of vcfR object
#' @param old.gt genotypes matrix (markers x individuals) with multiallelic markers 
#' @param P1 character with one of the parents ID
#' @param P2 character with other parent ID
#' 
get_alternatives <- function(fix, old.gt, P1, P2){
  p <- old.gt[,which(colnames(old.gt) %in% c(P1, P2))] 
  p <- paste0(p[,1],"/",p[,2])
  p <- strsplit(p, "/")
  p <- lapply(p, unique)
  p <- lapply(p, function(x) if(length(which(x=="NA" | x == ".")) > 0) as.numeric(x[-which(x=="NA" | x == ".")]) else as.numeric(x))
  p <- lapply(p, function(x) if(length(x) == 1)  c(x, x+1) else x)
  p <- lapply(p, function(x) sort(x + 1))
  
  alleles <- paste0(fix[,4], ",", fix[,5])
  alleles <- strsplit(alleles, ",")
  
  alleles.sele <- Map(function(x, y) {
    y[x]
  }, p, alleles)
  
  # Change reference allele
  fix[,4] <- sapply(alleles.sele, "[[", 1)
  
  # Change alternatives alleles
  fix[,5] <- sapply(alleles.sele, function(x) paste0(x[-1], collapse = ","))
  return(fix)
}

