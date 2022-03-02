#' Filter markers with segregation distortion. Accepts only biallelic markers.
#' 
#' @param vcf path to the vcf file
#' @param threshold alpha value to be considered
#' @param P1 parent 1 ID
#' @param P2 parent 2 ID
#' @param out.vcf output VCF file name
#' 
#' @import vcfR
#' @export 
segregation_test_vcf <- function(vcf, threshold, P1, P2, out.vcf = "filtered.vcf.gz"){
  obj <- read.vcfR(vcf)
  
  gt <- extract.gt(obj)
  
  gt <- gsub("[|]", "/", gt)
  gt <- gsub("1/0", "0/1", gt)
  
  parents <- which(colnames(gt) %in% c(P1, P2))
  
  mis <- which(is.na(gt[,parents[1]]) | is.na(gt[,parents[2]]))
  non_inf <- which(gt[,parents[1]] %in% c("0/0", "1/1") & gt[,parents[2]] %in% c("0/0", "1/1"))
  names(non_inf) <- rownames(gt)[non_inf]
  
  b3 <- which(gt[,parents[1]] == "0/1" & gt[,parents[2]] == "0/1")
  gt.b3 <- gt[b3, -parents]
  b3.p <- apply(gt.b3, 1, function(gt.marker) {
    c1 <- table(gt.marker)[1]
    c2 <- table(gt.marker)[2]
    c3 <- table(gt.marker)[3]
    if(is.na(c1) | is.na(c2) | is.na(c3)){
      return(NA)
    } else {
      qui <- chisq.test(as.vector(c(c1,c2,c3)), p=c(1/4, 1/2, 1/4), correct = FALSE, simulate.p.value = TRUE)
      return(qui$p.value)
    }
  })
  
  d1.10 <- which(gt[,parents[1]] == "0/1" & gt[,parents[2]] %in% c("0/0", "1/1"))
  d2.15 <- which(gt[,parents[1]] %in% c("0/0", "1/1") & gt[,parents[2]] == "0/1")
  gt.d <- gt[c(d1.10, d2.15), -parents]
  d.p <- apply(gt.d, 1, function(gt.marker) {
    c1 <- table(gt.marker)[1]
    c2 <- table(gt.marker)[2]
    if(is.na(c1) | is.na(c2)){
      return(NA)
    } else {
      qui <- chisq.test(as.vector(c(c1,c2)), p=c(1/2, 1/2), correct = FALSE, simulate.p.value = TRUE)
      return(qui$p.value)
    } 
  })
  
  p_values <- c(b3.p, d.p)
  
  threshold <- threshold/length(p_values)
  
  rm.mks <- names(which(p_values < threshold))
  rm.mks <- c(rm.mks, names(which(is.na(p_values))))
  rm.mks <- c(rm.mks, names(mis), names(non_inf))
  rm.mks.id <- which(rownames(gt) %in% rm.mks)
  
  new.obj <- obj
  new.obj@fix <- obj@fix[-rm.mks.id,]
  new.obj@gt <- obj@gt[-rm.mks.id,]
  
  write.vcf(new.obj, file = out.vcf)
}