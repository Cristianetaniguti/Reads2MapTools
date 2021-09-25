#' OneMap interface with updog package
#'
#' Uses alelle counts to reestimate genotypes with updog approach and 
#' stores the genotypes probabilities or for further multipoint 
#' analysis
#' 
#' @param vcf path and name of the vcf file
#' @param vcf.par Field of VCF that informs the depth of alleles
#' @param out_vcf path and name of the vcf file to be outputed from the generated onemap object. 
#' It is important to notice that only onemap informative markers are kept.
#' @param f1 f1 individual identification if f2 cross type
#' @param recovering logical defining if markers should be recovered from VCF
#' @param mean_phred the mean phred score of the sequencer technology
#' @param cores number of threads 
#' @param crosstype string defining the cross type, by now it supports only 
#' outcross and f2 intercross
#' @param depths list containing a matrix for ref and other for alt allele counts, samples ID in colum and markers ID in rows
#' @param parent1 parent 1 identification in vcfR object
#' @param parent2 parent 2 identification in vcfR objetc
#' @param global_error number from 0 to 1 defining the global error to be considered together 
#' with the genotype errors or the genotype probabilities or NULL to not considered any global error
#' @param use_genotypes_errors if \code{TRUE} the error probability of each genotype will be considered in emission function of HMM
#' @param use_genotype_probs if \code{TRUE} the probability of each possible genotype will be considered in emission function of HMM
#' @param rm_multiallelic if \code{TRUE} multiallelic markers will be removed from the output onemap object 
#' @param output_info_file define a name for the file with alleles information.
#' 
#' @return onemap object with genotypes updated 
#' 
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br} 
#' @seealso \code{\link[onemap]{extract_depth}} 
#'     \code{\link[onemap]{binom_genotype}} and 
#'     \url{https://github.com/dcgerard/updog}.
#'     
#'     
#' @references 
#'
#' Gerard, D., Ferrão L.F.V., Garcia, A.A.F., & Stephens, M. (2018). Harnessing 
#' Empirical Bayes and Mendelian Segregation for Genotyping Autopolyploids from 
#' Messy Sequencing Data. bioRxiv. doi: 10.1101/281550.
#'
#' @import foreach doParallel updog
#' @importFrom vcfR read.vcfR
#' @import onemap
#'   
#' @export
updog_genotype <- function(vcf=NULL,
                           vcf.par = c("AD", "DPR"),
                           out_vcf = NULL,
                           parent1="P1",
                           parent2="P2",
                           f1=NULL,
                           crosstype=NULL,
                           recovering = FALSE,
                           mean_phred = 20, 
                           cores = 2,
                           depths = NULL,
                           global_error = NULL,
                           use_genotypes_errors = TRUE,
                           use_genotypes_probs = FALSE,
                           rm_multiallelic = TRUE){
  
  # checks
  if(use_genotypes_errors & use_genotypes_probs){
    stop("You must choose only one of the offered approaches to be considered in emission function of HMM. `use_genotype_errors` or `use_genotype_probs`")
  } else if (!(use_genotypes_errors | use_genotypes_probs)){
    stop("You should choose one approach to be considered in emission function of HMM.`use_genotype_errors` or `use_genotype_probs`")
  }
  
  vcfR.object <- read.vcfR(vcf, verbose = F) 
  
  info_file_name <- tempfile()
  if(recovering){
    onemap.object <- onemap_read_vcfR(vcfR.object = vcfR.object,
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
    
  } else {
    onemap.object <- onemap_read_vcfR(vcfR.object = vcfR.object,
                                      parent1=parent1,
                                      parent2=parent2,
                                      f1=f1,
                                      cross=crosstype, 
                                      only_biallelic = F, 
                                      output_info_rds = info_file_name)
  }
  
  
  if(is.null(depths)){
    depth_matrix <- extract_depth(vcfR.object=vcfR.object,
                                  onemap.object=onemap.object,
                                  vcf.par=vcf.par,
                                  parent1=parent1,
                                  parent2=parent2,
                                  f1=f1,
                                  recovering=recovering)
  } else {
    p1 <- which(colnames(depths[[1]]) == parent1)
    p2 <- which(colnames(depths[[1]]) == parent2)
    if(is.null(f1)){
      palt <- depths[[2]][,c(p1,p2)]
      pref <- depths[[1]][,c(p1,p2)]
      oalt <- depths[[2]][,-c(p1,p2)]
      oref <- depths[[1]][,-c(p1,p2)]
      
      oalt <- oalt[,match(rownames(onemap.object$geno),colnames(oalt))]
      oref <- oref[,match(rownames(onemap.object$geno),colnames(oref))]
      
      psize <- palt + pref
      osize <- oalt + oref
      
      rownames(palt) <- rownames(pref) <- rownames(oalt)
      
    } else {
      f1i <- which(colnames(depths[[1]]) == f1)
      palt <- as.numeric(depths[[2]][,f1i])
      pref <- as.numeric(depths[[1]][,f1i])
      oalt <- depths[[2]][,-c(p1,p2,f1i)]
      oref <- depths[[1]][,-c(p1,p2,f1i)]
      
      oalt <- oalt[,match(rownames(onemap.object$geno),colnames(oalt))]
      oref <- oref[,match(rownames(onemap.object$geno),colnames(oref))]
      
      psize <- palt + pref
      osize <- oalt + oref
      
      names(palt) <- names(pref) <- rownames(oalt)
      
    }
    
    if(!recovering){
      palt <- palt[which(rownames(palt) %in% colnames(onemap.object$geno)),]
      pref <- pref[which(rownames(pref) %in% colnames(onemap.object$geno)),]
      psize <- psize[which(rownames(psize) %in% colnames(onemap.object$geno)),]
      oalt <- oalt[which(rownames(oalt) %in% colnames(onemap.object$geno)),]
      oref <- oref[which(rownames(oref) %in% colnames(onemap.object$geno)),]
      osize <- osize[which(rownames(osize) %in% colnames(onemap.object$geno)),]
    }
    
    depth_matrix <- list("palt"=palt, "pref"=pref, "psize"=psize, 
                         "oalt"=as.matrix(oalt), "oref"=as.matrix(oref), "osize"=as.matrix(osize), 
                         n.mks = dim(oalt)[1], n.ind = dim(oalt)[2],
                         inds = colnames(oalt), mks = rownames(oalt),
                         CHROM = sapply(strsplit(rownames(oalt), split="_"), "[",1),
                         POS = as.numeric(as.character(sapply(strsplit(rownames(oalt), split="_"), "[",2))),
                         onemap.object = onemap.object)
  }
  
  # Extract list objects                                                                                         
  for(i in 1:length(depth_matrix)){
    tempobj=depth_matrix[[i]]
    eval(parse(text=paste(names(depth_matrix)[[i]],"= tempobj")))
  }
  
  onemap_updog <- onemap.object
  
  # removing the multiallelics from the vcf
  multi <- grepl(",",vcfR.object@fix[,5])
  multi.mks <- vcfR.object@fix[,3][multi]
  
  if(any(multi)){ 
    idx <- which(rownames(palt) %in% multi.mks)
    palt <- palt[-idx,]
    pref <- pref[-idx,]
    psize <- psize[-idx,]
    oalt <- oalt[-idx,]
    oref <- oref[-idx,]
    osize <- osize[-idx,]
    n.mks = dim(oalt)[1]
    mks = rownames(oalt)
    CHROM = CHROM[-idx]
    POS = POS[-idx]
    
    # removing the multiallelics from the onemap object
    if(is(onemap.object, "outcross")){
      idx.multi <- which(!(onemap.object$segr.type %in% c("B3.7", "D1.10", "D2.15")))
      mult.obj <- split_onemap(onemap.object, mks= idx.multi)
      idx.bi <- which(onemap.object$segr.type %in% c("B3.7", "D1.10", "D2.15"))
      onemap.object <- split_onemap(onemap.object, mks= idx.bi)
    } else if(is(onemap.object,"f2")){
      idx.multi <- which(!(onemap.object$segr.type %in% c("A.H.B")))
      mult.obj <- split_onemap(onemap.object, mks= idx.multi)
      idx.bi <- which(onemap.object$segr.type %in% c("A.H.B"))
      onemap.object <- split_onemap(onemap.object, mks= idx.bi)
    }
  }
  
  # Missing data                                                                                                 
  if(is(onemap.object, "f2")){
    rm.mks <- which(psize ==0 | is.na(psize))
  }else{
    rm.mks <- which(psize[,1] ==0 | psize[,2] ==0 | is.na(psize[,2]) | is.na(psize[,1]) )
  }
  
  rm.mks1 <- vector()
  for(i in 1:n.mks){
    if(length(table(osize[i,]))==1)
      rm.mks1 <- c(rm.mks1,i)
  }
  rm.mks <- sort(unique(c(rm.mks, rm.mks1)))
  
  if(length(rm.mks) > 0){
    if(is(onemap.object, "f2")){
      psize <- psize[-rm.mks]
      pref <- pref[-rm.mks]
    } else {
      psize <- psize[-rm.mks,]
      pref <- pref[-rm.mks,]
    }
    osize <- osize[-rm.mks,]
    oref <- oref[-rm.mks,]
    n.mks <- n.mks - length(rm.mks)
    mks <- mks[-rm.mks]
    onemap_updog$CHROM <- CHROM[-rm.mks]
    onemap_updog$POS <- POS[-rm.mks]
  }
  
  idx <- which(osize==0)
  osize[idx] <- NA
  oref[idx] <- NA
  
  geno_matrix <- maxpostprob <- matrix(rep(NA,dim(osize)[2]*dim(osize)[1]),nrow=dim(osize)[1])
  P1 <- P2 <- rep(NA, n.mks)
  if(is(onemap.object, "f2")){
    cl <- parallel::makeCluster(as.numeric(cores))
    doParallel::registerDoParallel(cl = cl)
    gene_est <- foreach(i = 1:n.mks,
                        .combine = 'c',
                        .multicombine=TRUE,
                        .export = "flexdog") %dopar% {
                          fout <- flexdog(refvec  = oref[i,],
                                          sizevec = osize[i,],
                                          ploidy  = 2,
                                          p1ref = pref[i],
                                          p1size = psize[i],
                                          model = "s1")
                          list(fout)
                        }
    parallel::stopCluster(cl)
    
    P1 <- unlist(sapply(sapply(gene_est, "[", 8), "[", 1))
  } else if(is(onemap.object, "outcross")){
    cl <- parallel::makeCluster(as.numeric(cores))
    doParallel::registerDoParallel(cl = cl)
    gene_est <- foreach(i = 1:n.mks) %dopar% {
      ## fit flexdog                                                                                                                             
      fout <- updog::flexdog(refvec  = oref[i,],
                             sizevec = osize[i,],
                             ploidy  = 2,
                             p1ref = pref[i,2],
                             p1size = psize[i,2],
                             p2ref = pref[i,1],
                             p2size = psize[i,1],
                             model = "f1")
      fout
    }
    parallel::stopCluster(cl)
    P1 <- unlist(sapply(sapply(gene_est, "[", 8), "[", 2))
    P2 <- unlist(sapply(sapply(gene_est, "[", 8), "[", 1))
  }
  
  temp1 <- lapply(gene_est, "[[", 10)
  temp2 <- lapply(gene_est, "[[", 11)
  
  for(i in 1:n.mks){
    geno_matrix[i,] <- temp1[[i]]
    maxpostprob[i,] <- temp2[[i]]
  }
  
  postmat <- lapply(gene_est, "[", 6)
  
  genotypes_probs <- postmat[[1]]$postmat
  for(i in 2:length(postmat)) # Order: mk 1 1 1 ind 1 2 3
    genotypes_probs <- rbind(genotypes_probs, postmat[[i]]$postmat)
  
  # sort - order mk 1 2 3 ind 1 1 1 
  idx <- rep(1:dim(osize)[2], dim(osize)[1])
  genotypes_probs <- genotypes_probs[order(idx), ]
  
  if(is(onemap.object, "outcross")){
    rm.mk <- which(P1==0 & P2 == 2 | P2 == 0 & P1 == 0 | P1==2 & P2 == 2 | P1==2 & P2 ==0)
    
    if(length(rm.mk) > 0){
      cat("markers", length(mks[rm.mk]), "were reestimated as non-informative and removed of analysis \n")
      geno_matrix <- geno_matrix[-rm.mk,]
      maxpostprob <- maxpostprob[-rm.mk,]
      P1 <- P1[-rm.mk]
      P2 <- P2[-rm.mk]
      n.mks <- n.mks - length(rm.mk)
      mks <- mks[-rm.mk]
      genotypes_probs <- genotypes_probs[-c(rm.mk + rep(c(0:(dim(osize)[2]-1))*dim(osize)[1], each=length(rm.mk))),]
    }
    conv_geno <- matrix(rep(NA,dim(geno_matrix)[2]*dim(geno_matrix)[1]),nrow=dim(geno_matrix)[1])
    if(dim(conv_geno)[1] == 0) stop("All markers were filtered.\n")
    segr.type <- segr.type.num <- rep(NA, n.mks)
    conv_geno[which(geno_matrix==1)] <- 2
    # B3.7                                                                                                         
    idx <- which(P1==1 & P2==1)
    conv_geno[idx,][which(geno_matrix[idx,]==0)] <- 3
    conv_geno[idx,][which(geno_matrix[idx,]==2)] <- 1
    segr.type[idx] <- "B3.7"
    segr.type.num[idx] <- 4
    
    # D2.15                                                                                                        
    idx <- which(P1==0 & P2==1)
    conv_geno[idx,][which(geno_matrix[idx,]==0)] <- 1
    conv_geno[idx,][which(geno_matrix[idx,]==2)] <- 0
    segr.type[idx] <- "D2.15"
    segr.type.num[idx] <- 7
    
    idx <- which(P1==2 & P2==1)
    conv_geno[idx,][which(geno_matrix[idx,]==0)] <- 0
    conv_geno[idx,][which(geno_matrix[idx,]==2)] <- 1
    segr.type[idx] <- "D2.15"
    segr.type.num[idx] <- 7
    # D1.10                                                                                                        
    idx <- which(P1==1 & P2==0)
    conv_geno[idx,][which(geno_matrix[idx,]==0)] <- 1
    conv_geno[idx,][which(geno_matrix[idx,]==2)] <- 0
    segr.type[idx] <- "D1.10"
    segr.type.num[idx] <- 6
    
    idx <- which(P1==1 & P2==2)
    conv_geno[idx,][which(geno_matrix[idx,]==0)] <- 0
    conv_geno[idx,][which(geno_matrix[idx,]==2)] <- 1
    segr.type[idx] <- "D1.10"
    segr.type.num[idx] <- 6
    
  } else {
    rm.mk <- which(P1!=1)
    if(length(rm.mk) > 0){
      cat("markers", length(mks[rm.mk]), "were estimated as non-informative and removed of analysis \n")
      geno_matrix <- geno_matrix[-rm.mk,]
      maxpostprob <- maxpostprob[-rm.mk,]
      P1 <- P1[-rm.mk]
      P2 <- P2[-rm.mk]
      n.mks <- n.mks - length(rm.mk)
      mks <- mks[-rm.mk]
      genotypes_probs <- genotypes_probs[-c(rm.mk + rep(c(0:(dim(osize)[2]-1))*dim(osize)[1], each=length(rm.mk))),]
    }
    conv_geno <- matrix(rep(NA,dim(geno_matrix)[2]*dim(geno_matrix)[1]),nrow=dim(geno_matrix)[1])
    if(dim(conv_geno)[1] == 0) stop("All markers were filtered.\n")
    segr.type <- segr.type.num <- rep(NA, n.mks)
    conv_geno[which(geno_matrix==1)] <- 2
    
    # A.H.B                                                                                                        
    idx <- which(P1==1)
    conv_geno[idx,][which(geno_matrix[idx,]==0)] <- 3
    conv_geno[idx,][which(geno_matrix[idx,]==2)] <- 1
    segr.type[idx] <- "A.H.B"
    segr.type.num[idx] <- 1
  }
  
  P1 <- recode_parents(P1)
  P2 <- recode_parents(P2)
  
  conv_geno <-  t(conv_geno)
  conv_geno[which(is.na(conv_geno))] <- 0
  
  comp <- which(colnames(onemap.object$geno) %in% mks)
  onemap_updog$CHROM <- onemap.object$CHROM[comp]
  onemap_updog$POS <- onemap.object$POS[comp]
  
  comp1 <- which(depth_matrix$mks %in% mks)
  onemap_updog$CHROM <- depth_matrix$CHROM[comp1]
  onemap_updog$POS <- as.numeric(as.character(depth_matrix$POS))[comp1]
  
  cat("New onemap object contains", length(mks), "biallelic markers\n")
  
  maxpostprob <- t(1- maxpostprob)
  colnames(conv_geno)  <- colnames(maxpostprob) <- mks
  inds <- rownames(conv_geno)  <- rownames(maxpostprob) <- rownames(onemap.object$geno)
  
  onemap_updog$geno <- conv_geno
  onemap_updog$n.ind <- length(inds)
  onemap_updog$n.mar <- length(mks)
  onemap_updog$segr.type.num <- segr.type.num
  onemap_updog$segr.type <- segr.type
  
  if(use_genotypes_probs){
    onemap_updog.new <- create_probs(input.obj = onemap_updog,
                                     genotypes_probs = genotypes_probs,
                                     global_error = global_error)
  } else if(use_genotypes_errors){
    onemap_updog.new <- create_probs(input.obj = onemap_updog,
                                     genotypes_errors = maxpostprob,
                                     global_error = global_error)
  } else if(!is.null(global_error)){
    onemap_updog.new <- create_probs(input.obj = onemap_updog,
                                     global_error = global_error)
  }
  
  if(!rm_multiallelic){
    if(length(multi.mks) > 0)
      onemap_updog.new <- combine_onemap(onemap_updog.new, mult.obj)
  }
  
  if(!is.null(out_vcf)){
    onemap_write_vcfR(onemap.object = onemap_updog.new, 
                      out_vcf = out_vcf, 
                      input_info_rds = info_file_name,
                      probs = genotypes_probs, 
                      parent1.id = parent1, 
                      parent2.id = parent2, 
                      parent1.geno = P1, 
                      parent2.geno = P2)
  }
  
  structure(onemap_updog.new)
}

##' Plot a barplot with error probabilities values
##' 
##' @param onemap.obj an object of class \code{onemap} coming from \code{read_onemap}, 
##' \code{read_mapmaker}, \code{onemap_read_vcfR}, \code{updog_genotype} functions
##' 
##' @param mk.type a TRUE/FALSE value to define if genotypes will colored by marker type
##' 
##' @param select.ind a string defining specific individuals. The graphic will only contains 
##' error probability information of this individuals.
##' 
##' @param select.mk a string defining specific markers. The graphic will only contains 
##' error probability information of this markers.
##' 
##' @param n.int
##' 
##' @examples 
##' 
##' \dontrun{
##' data(onemap_example_out)
##' p <- plot_error_dist(onemap_example_out)
##' ggplot2::ggsave("errors_out.jpg", p width=7, height=4, dpi=300)
##' }
##' @import ggplot2
##' 
##' @export
plot_error_dist <- function(onemap.obj = NULL, mk.type = TRUE, 
                            select.ind = NULL, select.mk = NULL, 
                            n.int= NULL){
  
  #Do checks
  M <- reshape2::melt(onemap.obj$error)
  M <- cbind(M, type = rep(onemap.obj$segr.type, each = onemap.obj$n.ind))
  
  if(!is.null(select.ind)){
    idx <- which(as.character(M$Var1) %in% select.ind)
    if(length(idx)==0){
      stop("This individual does not exist in the given onemap object\n")
    } else{ M <- M[idx,] }
  }
  
  if(!is.null(select.mk)){
    idx <- which(as.character(M$Var2) %in% select.mk)
    if(length(idx)==0){
      stop("This marker does not exist in the given onemap object\n")
    }
    M <- M[idx,]
  }
  
  if(is.null(n.int)){
    if(nclass.FD(M$value) <= 10000000){
      breaks <- pretty(range(M$value), n = nclass.FD(M$value), min.n = 1)
      binwidth <- breaks[2]-breaks[1]
    } else {
      stop("Choose a n.int value for your graphic\n")
    }
  } else {
    breaks <- pretty(range(M$value), n = n.int, min.n = 1)
    binwidth <- breaks[2]-breaks[1]
  }
  
  if(mk.type){
    p <- ggplot(M, aes(value, fill = as.factor(type))) +
      geom_histogram(binwidth = binwidth)
  } else {
    p <- ggplot(M, aes(value)) +
      geom_histogram(binwidth = binwidth)
  }
  
  p <- p +  xlab("error probability") + ylab("count") + labs(fill = "Marker type") + scale_x_continuous(limits = c(0, max(M$value)))
  return(p)
}

recode_parents <- function(x) {
  x[which(x==1)] <- "0/1"
  x[which(x==0)] <- "1/1"
  x[which(x==2)] <- "0/0"
  return(x)
}
