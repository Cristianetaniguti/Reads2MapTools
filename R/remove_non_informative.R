#' Filter non-informative markers. Accepts only biallelic markers.
#' 
#' @param vcf path to the vcf file
#' @param P1 parent 1 ID
#' @param P2 parent 2 ID
#' @param out.vcf output VCF file name
#' @param replaceAD logical, if true it replaces AD by missing data when GT is missing
#' 
#' @import vcfR
#' @export 
remove_non_informative <- function(vcf, P1, P2, out.vcf = "filtered.vcf.gz", replaceAD = TRUE){
  obj <- read.vcfR(vcf)
  gt <- extract.gt(obj)
  
  if(replaceAD){
    dp <- extract.gt(obj, element = "DP")
    idx <- which(dp==0)
    if(length(idx) >0)  dp[idx] <- NA else idx <- is.na(dp)
    if(length(idx) > 0){
      # Replace AD and DP by 0 when GT is missing
      obj@gt[,-1][which(is.na(gt))] <- obj@gt[,-1][idx][1]
    }
  }
  
  # Remove phase information
  idx <- grep("[|]", gt)
  if(length(idx) >0){
    gt <- gsub("[|]", "/", gt)
    gt[idx] <- sapply(strsplit(gt[idx], "/"), function(x) paste0(sort(x), collapse = "/"))
  }
  
  parents <- which(colnames(gt) %in% c(P1, P2))
  gt.p1 <- gt[,parents[1]]
  gt.p2 <- gt[,parents[2]]
  
  mis <- which(is.na(gt.p1) | is.na(gt.p2))
  
  non_inf <- which(((grepl("0",gt.p1) & !grepl("1",gt.p1)) | ((grepl("1",gt.p1) & !grepl("0",gt.p1)))) &
                     ((grepl("0",gt.p2) & !grepl("1",gt.p2)) | ((grepl("1",gt.p2) & !grepl("0",gt.p2)))))
  
  new.obj <- obj
  new.obj@fix <- obj@fix[-non_inf,]
  new.obj@gt <- obj@gt[-non_inf,]
  
  write.vcf(new.obj, file = out.vcf)
}
