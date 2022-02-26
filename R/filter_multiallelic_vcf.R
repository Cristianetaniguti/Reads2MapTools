#' Filter vcf file keeping only alleles found in parents haplotypes
#' 
#' @param vcf path to vcf file
#' @param P1 character with one of the parents ID
#' @param P2 character with other parent ID
#' @param max.missing maximum fraction of missing data allowed by marker
#' @param vcf.out character defining output file name
#' 
#' @import vcfR
#' 
#' @export
filter_multi_vcf <- function(vcf.file, P1, P2, max.missing = NULL, vcf.out = "filtered.vcf.gz"){
  vcf <- read.vcfR(vcf.file)

  gt <- extract.gt(vcf)
  
  filt.gt <- filter_multi(gt, P1, P2)
  up.fix <- get_alternatives(vcf@fix, gt, P1, P2)
  
  vcf.new <- vcf
  vcf.new@fix <- up.fix
  vcf.new@gt <- cbind(FORMAT="GT", filt.gt)
  
  if(!is.null(max.missing)){
    mis <- apply(vcf.new@gt, 1, function(x) sum(is.na(x))/length(x))
    idx <- which(mis > 0.25)
    
    vcf.new@gt <- vcf.new@gt[-idx,]
    vcf.new@fix <- vcf.new@fix[-idx,]
  }
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
  
  mk.split <- apply(gt.multi, 1, function(x) strsplit(x, "/"))
  
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
#' @param gt genotypes matrix (markers x individuals) with multiallelic markers 
#' @param P1 character with one of the parents ID
#' @param P2 character with other parent ID
#' 
get_alternatives <- function(fix, gt, P1, P2){
  p <- gt[,which(colnames(gt) %in% c(P1, P2))]
  p <- paste0(p[,1],"/",p[2])
  p <- strsplit(p, "/")
  p <- lapply(p, unique)
  p <- lapply(p, function(x) if(length(which(x=="NA" | x == ".")) > 0) x[-which(x=="NA" | x == ".")] else x)
  p <- lapply(p, function(x) sort(as.numeric(x) + 1))
  
  alleles <- paste0(fix[,4], ",", fix[,5])
  alleles <- strsplit(alleles, ",")
  
  alleles.sele <- mapply(function(x, y) {
    y[x]
  }, p, alleles)
  
  # Change reference allele
  fix[,4] <- sapply(alleles.sele, "[", 1)
  
  # Change alternatives alleles
  fix[,5] <- sapply(alleles.sele, function(x) paste0(x[-1], collapse = ","))
  return(fix)
}
