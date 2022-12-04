# Contains supermassa_error, depth_prepare and supermassa_parallel

##' Runs SuperMassa for multiple markers using doParallel packages
##' Updates OneMap object in genotypes, error probabilitie and marker type
##' Also removes non-informative markers according to SuperMassa genotyping.
##' @param vcf path and name of the vcf file
##' @param vcf.par Field of VCF that informs the depth of alleles
#' @param out_vcf path and name of the vcf file to be outputed from the generated onemap object. 
#' It is important to notice that only onemap informative markers are kept.
##' @param parent1 parent 1 identification in vcfR object
##' @param parent2 parent 2 identification in vcfR objetc
##' @param f1 f1 individual identification if f2 cross type
##' @param crosstype string defining the cross type, by now it supports only 
##' outcross and f2 intercross
##' @param recovering logical defining if markers should be recovered from VCF
##' @param cores number of threads 
##' @param depths list containing a matrix for ref and other for alt allele counts, samples ID in colum and markers ID in rows
##' @param global_error number from 0 to 1 defining the global error to be considered together 
##' with the genotype errors or the genotype probabilities or NULL to not considered any global error
##' @param use_genotypes_errors if \code{TRUE} the error probability of each genotype will be considered in emission function of HMM
##' @param use_genotype_probs if \code{TRUE} the probability of each possible genotype will be considered in emission function of HMM
##' @param output_info_file define a name for the file with alleles information.
##' 
##' @import doParallel parallel
##' 
##' @importFrom matrixStats logSumExp
##' @importFrom vcfR read.vcfR
##' @import onemap
##' 
##' @export
supermassa_genotype <- function(vcf=NULL,
                                vcf.par = c("AD", "DPR"),
                                out_vcf = NULL,
                                parent1="P1",
                                parent2="P2",
                                crosstype=NULL,
                                f1=NULL,
                                recovering = FALSE,
                                cores = 2,
                                depths = NULL,
                                global_error = NULL,
                                use_genotypes_errors = TRUE,
                                use_genotypes_probs = FALSE,
                                rm_multiallelic = TRUE,
                                info_file_name=NULL){
  
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
    extracted_depth <- extract_depth(vcfR.object=vcfR.object,
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
    
    if(recovering==FALSE){
      palt <- palt[which(rownames(palt) %in% colnames(onemap.object$geno)),]
      pref <- pref[which(rownames(pref) %in% colnames(onemap.object$geno)),]
      psize <- psize[which(rownames(psize) %in% colnames(onemap.object$geno)),]
      oalt <- oalt[which(rownames(oalt) %in% colnames(onemap.object$geno)),]
      oref <- oref[which(rownames(oref) %in% colnames(onemap.object$geno)),]
      osize <- osize[which(rownames(osize) %in% colnames(onemap.object$geno)),]
    }
    
    extracted_depth <- list("palt"=palt, "pref"=pref, "psize"=psize, 
                            "oalt"=as.matrix(oalt), "oref"=as.matrix(oref), "osize"=as.matrix(osize), 
                            n.mks = dim(oalt)[1], n.ind = dim(oalt)[2],
                            inds = colnames(oalt), mks = rownames(oalt),
                            CHROM = sapply(strsplit(rownames(oalt), split="_"), "[",1),
                            POS = sapply(strsplit(rownames(oalt), split="_"), "[",2),
                            onemap.object = onemap.object)
  }
  
  # removing the multiallelics from the vcf
  multi <- grepl(",",vcfR.object@fix[,5])
  multi.mks <- vcfR.object@fix[,3][multi]
  
  if(any(multi)){ 
    idx <- which(rownames(extracted_depth$palt) %in% multi.mks)
    extracted_depth$palt <- extracted_depth$palt[-idx,]
    extracted_depth$pref <- extracted_depth$pref[-idx,]
    extracted_depth$psize <- extracted_depth$psize[-idx,]
    extracted_depth$oalt <- extracted_depth$oalt[-idx,]
    extracted_depth$oref <- extracted_depth$oref[-idx,]
    extracted_depth$osize <- extracted_depth$osize[-idx,]
    extracted_depth$n.mks = dim(extracted_depth$oalt)[1]
    extracted_depth$mks = rownames(extracted_depth$oalt)
    extracted_depth$CHROM = extracted_depth$CHROM[-idx]
    extracted_depth$POS = extracted_depth$POS[-idx]
    
    # removing the multiallelics from the onemap object
    if(is(onemap.object, "outcross")){
      idx.multi <- which(!(onemap.object$segr.type %in% c("B3.7", "D1.10", "D2.15")))
      mult.obj <- split_onemap(onemap.object, mks= idx.multi)
      idx.bi <- which(onemap.object$segr.type %in% c("B3.7", "D1.10", "D2.15"))
      onemap.object <- split_onemap(onemap.object, mks= idx.bi)
    } else if(is(onemap.object,"f2 intercross")){
      idx.multi <- which(!(onemap.object$segr.type %in% c("A.H.B")))
      mult.obj <- split_onemap(onemap.object, mks= idx.multi)
      idx.bi <- which(onemap.object$segr.type %in% c("A.H.B"))
      onemap.object <- split_onemap(onemap.object, mks= idx.bi)
    }
    extracted_depth$onemap.object <- onemap.object
  }
  
  # NA replaced by 0 - not ideal - Bugfix!
  for(i in 1:6)
    extracted_depth[[i]][is.na(extracted_depth[[i]])] <- 0
  
  prepared_depth <- depth_prepare(extracted_depth)
  obj.class <- class(extracted_depth$onemap.object)[2]
  
  cl <- parallel::makeCluster(as.numeric(cores))
  registerDoParallel(cl)
  clusterExport(cl, c("supermassa_parallel"))
  result <- parLapply(cl, prepared_depth, function(x) supermassa_parallel(supermassa_4parallel = x, class=obj.class))
  parallel::stopCluster(cl)
  mks <- unlist(lapply(result, "[", 1))
  rm.mk <- unlist(lapply(result, "[", 2))
  mks.type <- unlist(lapply(result, "[", 3))
  geno <- matrix(unlist(lapply(result, "[", 4)), ncol = extracted_depth$n.mks, byrow = FALSE)
  error <- lapply(result, "[[", 5)
  error <- as.matrix(do.call(rbind, error[!as.logical(rm.mk)]))
  pgeno <- lapply(result, "[[", 6)
  pgeno <- as.matrix(do.call(rbind, pgeno[!as.logical(rm.mk)]))
  
  # parents geno
  p1 <- apply(pgeno[,1:2], 1, recode_parents_sm)
  p2 <- apply(pgeno[,3:4], 1, recode_parents_sm)
  names(p1) <- names(p2) <- mks[!rm.mk]
  
  geno[which(is.na(geno))] <- 0
  error[which(is.na(error))] <- 1
  
  colnames(geno)  <- extracted_depth$mks
  rownames(geno)  <- extracted_depth$inds
  
  n.mar <- length(rm.mk) - sum(rm.mk)
  
  onemap_supermassa <- extracted_depth$onemap.object
  
  #error[,2] <- factor(error[,2], levels = extracted_depth$inds)
  
  # set ord mk 1 2 3 ind 1 1 1
  idx <- rep(1:length(unique(error[,2])), length(unique(error[,1])))
  error <- error[order(idx),]
  onemap_supermassa$error <- as.matrix(apply(error[,3:5],2,as.numeric))
  colnames(onemap_supermassa$error) <- NULL
  
  if (any(rm.mk == 1)){
    onemap_supermassa$geno <- geno[,-which(rm.mk==1)]
    onemap_supermassa$segr.type <- unname(mks.type[-which(rm.mk==1)])
    onemap_supermassa$CHROM <- extracted_depth$CHROM[-which(rm.mk==1)]
    onemap_supermassa$POS <- as.numeric(extracted_depth$POS[-which(rm.mk==1)])
  } else{
    onemap_supermassa$geno <- geno
    onemap_supermassa$segr.type <- unname(mks.type)
    onemap_supermassa$CHROM <- extracted_depth$CHROM
    onemap_supermassa$POS <- as.numeric(extracted_depth$POS)
  }
  
  onemap_supermassa$n.mar <- n.mar
  
  mks.type.num <- onemap_supermassa$segr.type
  mks.type.num[which(onemap_supermassa$segr.type == "B3.7")] <- 4
  mks.type.num[which(onemap_supermassa$segr.type == "D2.15")] <- 7
  mks.type.num[which(onemap_supermassa$segr.type == "D1.10")] <- 6
  mks.type.num[which(onemap_supermassa$segr.type == "A.H.B")] <- 1
  
  onemap_supermassa$segr.type.num <- mks.type.num
  
  # Genotypes percent changed - markers have repeated names (splitted multiallelics)
  # idx <- match(colnames(onemap_supermassa$geno), colnames(extracted_depth$onemap.object$geno))
  # idx2 <- match(colnames(extracted_depth$onemap.object$geno), colnames(onemap_supermassa$geno))
  # if(length(idx) > 0){
  #   diffe <- sum(!extracted_depth$onemap.object$geno[,idx] == onemap_supermassa$geno[,idx2])
  # } else {diffe <- 0}
  # perc <- (diffe/length(extracted_depth$onemap.object$geno))*100
  # cat(perc, "% of the genotypes were changed by this approach\n")
  # 
  # # Recovery markers
  # cat(sum(!colnames(onemap_supermassa$geno) %in% colnames(extracted_depth$onemap.object$geno)),"were recovered from vcf\n")
  # 
  # # Removed markers
  # cat(sum(!colnames(extracted_depth$onemap.object$geno) %in% colnames(onemap_supermassa$geno)),
  #     "markers from old onemap object were considered non-informative and removed of analysis\n")
  
  maxpostprob <- unlist(apply(onemap_supermassa$error, 1, function(x) x[which.max(x)]))
  maxpostprob <- matrix(maxpostprob, nrow = onemap_supermassa$n.ind, ncol = onemap_supermassa$n.mar, byrow = F)
  colnames(maxpostprob) <- colnames(onemap_supermassa$geno)
  rownames(maxpostprob) <- rownames(onemap_supermassa$geno)
  
  if(use_genotypes_probs){
    onemap_supermassa.new <- create_probs(input.obj = onemap_supermassa,
                                          genotypes_probs = onemap_supermassa$error,
                                          global_error = global_error)
  } else if(use_genotypes_errors){
    onemap_supermassa.new <- create_probs(input.obj = onemap_supermassa,
                                          genotypes_errors = 1- maxpostprob,
                                          global_error = global_error)
  } else if(!is.null(global_error)){
    onemap_supermassa.new <- create_probs(input.obj = onemap_supermassa,
                                          global_error = global_error)
  }
  
  if(!rm_multiallelic){ ## BUGFIX! Fix probabilities if rm_multiallelics = FALSE
    if(length(multi.mks) > 0)
      onemap_supermassa.new <- combine_onemap(onemap_supermassa.new, mult.obj)
  }
  
  # AD matrix
  ref <- cbind(extracted_depth$pref, extracted_depth$oref)
  ref[is.na(ref)] <- "."
  alt <- cbind(extracted_depth$psize - extracted_depth$pref, extracted_depth$osize - extracted_depth$oref)
  alt[is.na(alt)] <- "."
  ad_matrix <- matrix(paste0(ref, ",", alt), nrow = nrow(ref))
  colnames(ad_matrix) <- colnames(ref)
  rownames(ad_matrix) <- rownames(ref)
  
  # sort - order mk 1 1 1 ind 1 2 3 
  idx <- rep(1:onemap_supermassa$n.mar, onemap_supermassa$n.ind)
  genotypes_probs <- onemap_supermassa$error[order(idx), ]
  
  if(!is.null(out_vcf)){
    onemap_write_vcfR(onemap.object = onemap_supermassa.new, 
                      out_vcf = out_vcf, 
                      input_info_rds = info_file_name,
                      probs = genotypes_probs, 
                      parent1.geno = p1, 
                      parent2.geno = p2,
                      parent1.id = parent1,
                      parent2.id = parent2,
                      ad_matrix = ad_matrix)
  }
  
  return(onemap_supermassa.new)
}

##' Creates list with first level refering to markers and second level with
##' marker name, offspring names and alternative and reference allele counts.
##' This function was created to parallelize the SuperMassa process. Its output
##' is the input of supermassa_error
##' @export
depth_prepare <- function(extracted_depth){
  mylist <- rep(list(NA),extracted_depth$n.mks)
  for(i in 1:extracted_depth$n.mks){
    mk.name <- extracted_depth$mks[i]
    inds <- extracted_depth$inds
    odepth <- data.frame(extracted_depth$oref[i,], extracted_depth$oalt[i,])
    if(class(extracted_depth$onemap.object)[2]=="f2"){
      pdepth <- data.frame(extracted_depth$pref[i], extracted_depth$palt[i])
    } else {
      pdepth <- data.frame(extracted_depth$pref[i,], extracted_depth$palt[i,])
    }
    mylist[[i]] <- list(mk.name, inds,odepth,pdepth)
  }
  return(mylist)
}


##' Use output list from depth_prepare and run Supermassa for each marker (first level of the list)
##' Extract from the SuperMassa output, the offspring genotype, the marker type, and the error
##' probabilities
##' @export
supermassa_parallel <- function(supermassa_4parallel, class=NULL){
  n.ind <- length(supermassa_4parallel[[2]])
  all.ind.names <- supermassa_4parallel[[2]]
  if(length(which(supermassa_4parallel[[4]][,1] == 0 & supermassa_4parallel[[4]][,2]== 0)) > 0) {
    # Marker with parents missing data
    rm.mk <- 1 #index to markers to be removed
    geno_one <- prob_error <- rep(NA, n.ind)
    mk.type <- NA
    probs <- NA
    pgeno <- NA
  }else if(all(supermassa_4parallel[[3]]==0)){
    # Marker with missing data in all progeny
    rm.mk <- 1 #index to markers to be removed
    geno_one <- prob_error <- rep(NA, n.ind)
    mk.type <- NA
    probs <- NA
    pgeno <- NA
  } else {
    rm.mk <- 0
    write.table(supermassa_4parallel[[3]], file = paste0("odepth_temp", supermassa_4parallel[[1]],".txt"), quote = FALSE, col.names = FALSE)
    if(class=="f2"){
      write.table(rbind(supermassa_4parallel[[4]],supermassa_4parallel[[4]]), file = paste0("pdepth_temp", supermassa_4parallel[[1]],".txt"), quote = FALSE, col.names = FALSE)
    } else {
      write.table(supermassa_4parallel[[4]], file = paste0("pdepth_temp", supermassa_4parallel[[1]],".txt"), quote = FALSE, col.names = FALSE)
    }
    odepth_temp <- paste0("odepth_temp", supermassa_4parallel[[1]],".txt")
    pdepth_temp <- paste0("pdepth_temp", supermassa_4parallel[[1]],".txt")
    out_file <- paste0("out_prob",supermassa_4parallel[[1]],".txt")
    
    command_mass <- paste("python", paste0(system.file(package = "Reads2MapTools"),"/python_scripts/","SuperMASSA.py"),
                          "--inference f1 --file", paste0(getwd(),"/",odepth_temp), "--ploidy_range 2",
                          " --f1_parent_data", paste0(getwd(),"/",pdepth_temp), 
                          " --print_genotypes --naive_posterior_reporting_threshold 0",
                          "--save_geno_prob_dist", paste0(getwd(),"/",out_file))
    
    SuperMASSA.output<-system(command_mass, intern=TRUE)
    
    file.remove(odepth_temp,pdepth_temp)
    
    ############ Start code from Molli
    
    ## Assessing estimated ploidy level
    est.ploidy<-as.numeric(strsplit(SuperMASSA.output[4], split = " |,")[[1]][6])
    #cat(" (ploidy ", est.ploidy, ") ", sep = "")
    
    ## Assessing estimated parental dosage
    dpdq <- as.numeric(strsplit(x = SuperMASSA.output[4], split = "\\(|,")[[1]][c(5,8)])
    
    ## Reading SuperMASSA output
    d<-scan(file=out_file, what="character", nlines=1, quiet = TRUE)
    P<-matrix(as.numeric(gsub("[^0-9]", "", unlist(d))), ncol=2, byrow=TRUE)
    A<-read.table(file=out_file, skip=1)
    file.remove(out_file)
    M<-matrix(NA, nrow = length(all.ind.names), ncol = est.ploidy+1,
              dimnames = list(all.ind.names, c(0:est.ploidy)))
    M.temp<-apply(A[,2:(est.ploidy+2)], 2, parse.geno)
    if(is(M.temp, "numeric")) M.temp <- t(as.matrix(M.temp))
    M.temp<-M.temp[,order(P[,2], decreasing=TRUE)]
    if(is(M.temp,"numeric")) M.temp <- t(as.matrix(M.temp))
    dimnames(M.temp)<-list(as.character(A[,1]), c(0:est.ploidy))
    M[rownames(M.temp),]<-M.temp
    mrk1<-apply(M, 1, function(x) exp(x-matrixStats::logSumExp(x)))
    
    ## Filling NAs with mendelian expectations
    ## z<-which(apply(mrk1, 2, function(x) any(is.na(x))))
    ## if(length(z) > 0)
    ## {
    ##   if((dpdq[1] == 1 & dpdq[2] == 2) |(dpdq[1] == 2 & dpdq[2] == 1)){
    ##     a <- c(aa=0,Aa=0.5,AA=0.5)
    ##   } else if((dpdq[1] == 1 & dpdq[2] == 0) | (dpdq[1] == 0 & dpdq[2] == 1)){
    ##     a <- c(aa=0.5,Aa=0.5,AA=0)
    ##   } else if((dpdq[1] == 1 & dpdq[2] == 1)){
    ##     a <- c(aa=0.25, Aa=0.5, AA=0.25)
    ##   } else if((dpdq[1] == 2 & dpdq[2] == 0) | (dpdq[1] == 0 & dpdq[2] == 2)){
    ##     a <- c(aa=0, Aa=1, AA=0)
    ##   } else if((dpdq[1] == 2 & dpdq[2] == 2)){
    ##     a <- c(aa=0, Aa=0, AA=1)
    ##   } else if((dpdq[1] == 0 & dpdq[2] == 0)){
    ##     a <- c(aa=1, Aa=0, AA=0)  
    ##   }
    ##   for(k in names(z))
    ##     mrk1[,z]<-a
    
    ##   rownames(mrk1)<-names(a)<-0:est.ploidy
    ##   mrk1[names(which(a==0)),][]<-0
    ##   for(k1 in 1:ncol(mrk1))
    ##   {
    ##     if(any(mrk1[names(which(a!=0)),k1] > 0.5))
    ##       mrk1[,k1] <- mrk1[,k1]/sum(mrk1[,k1])
    ##     else mrk1[,k1] <- a
    ##   }
    ## }
    
    probs <- t(mrk1)
    probs <- data.frame(mrk = supermassa_4parallel[[1]],
                        ind = rownames(probs),
                        probs)
    rownames(probs)<-NULL
    colnames(probs)<-c("mrk", "ind", 0:est.ploidy)
    
    pgeno_temp <- strsplit(SuperMASSA.output[grep("ploidy", SuperMASSA.output)], split = ":")
    pgeno_temp <- unlist(strsplit(pgeno_temp[[1]][length(pgeno_temp[[1]])], split = ","))
    
    # Code for OneMap
    pgeno <- rep(NA, 4)
    pgeno[grep("0", pgeno_temp)] <- 0
    pgeno[grep("1", pgeno_temp)] <- 1
    pgeno[grep("2", pgeno_temp)] <- 2
    
    if(class=="outcross"){
      # only biallelic markers and diploid
      if((pgeno[1] == 0 | pgeno[1] == 2) & (pgeno[3] == 2 | pgeno[3] == 0)){
        # non-informative markers
        rm.mk <- 1
        geno_one <- prob_error <- rep(NA, n.ind)
        mk.type <- NA
      } else {
        rm.mk <- 0
        ind_geno <- unlist(lapply(strsplit(SuperMASSA.output[grep("\t", SuperMASSA.output)], split = "\t"), "[", 1))
        ind_geno <- gsub('.{1}$', '', ind_geno)
        geno <- unlist(lapply(strsplit(SuperMASSA.output[grep("\t", SuperMASSA.output)], split = "\t"), "[", 2))
        
        ind_geno_tot <- cbind(supermassa_4parallel[[2]], rep(NA, n.ind))
        ind_geno_tot[,2][match(ind_geno,ind_geno_tot[,1])] <- geno
        geno_one <- rep(NA, n.ind)
        
        # Code for OneMap
        idx <- which(ind_geno_tot[,2] == "(1, 1) <br>")
        geno_one[idx] <- 2
        
        mk.type <- vector()
        if(pgeno[1] == 0){
          mk.type <- c(mk.type,"D2.15")
          idx <- which(ind_geno_tot[,2] == "(2, 0) <br>")
          geno_one[idx] <- 0
          idx <- which(ind_geno_tot[,2] == "(0, 2) <br>")
          geno_one[idx] <- 1
        } else if(pgeno[2] == 0){
          mk.type <- c(mk.type,"D2.15")
          idx <- which(ind_geno_tot[,2] == "(2, 0) <br>")
          geno_one[idx] <- 1
          idx <- which(ind_geno_tot[,2] == "(0, 2) <br>")
          geno_one[idx] <- 0
        } else if(pgeno[3] == 0){
          mk.type <- c(mk.type,"D1.10")
          idx <- which(ind_geno_tot[,2] == "(2, 0) <br>")
          geno_one[idx] <- 0
          idx <- which(ind_geno_tot[,2] == "(0, 2) <br>")
          geno_one[idx] <- 1
        } else if(pgeno[3] == 2){
          mk.type <- c(mk.type,"D1.10")
          idx <- which(ind_geno_tot[,2] == "(2, 0) <br>")
          geno_one[idx] <- 1
          idx <- which(ind_geno_tot[,2] == "(0, 2) <br>")
          geno_one[idx] <- 0
        } else{
          mk.type <- c(mk.type,"B3.7")
          idx <- which(ind_geno_tot[,2] == "(2, 0) <br>")
          geno_one[idx] <- 1
          idx <- which(ind_geno_tot[,2] == "(0, 2) <br>")
          geno_one[idx] <- 3
        }
      }
    } else {
      
      # only biallelic markers and diploid
      if(pgeno[1] != 1 | pgeno[3] != 1){
        # non-informative markers
        rm.mk <- 1
        geno_one <- prob_error <- rep(NA, n.ind)
        mk.type <- NA
      } else {
        rm.mk <- 0
        ind_geno <- unlist(lapply(strsplit(SuperMASSA.output[grep("\t", SuperMASSA.output)], split = "\t"), "[", 1))
        ind_geno <- gsub('.{1}$', '', ind_geno)
        geno <- unlist(lapply(strsplit(SuperMASSA.output[grep("\t", SuperMASSA.output)], split = "\t"), "[", 2))
        
        ind_geno_tot <- cbind(supermassa_4parallel[[2]], rep(NA, n.ind))
        ind_geno_tot[,2][match(ind_geno,ind_geno_tot[,1])] <- geno
        geno_one <- rep(NA, n.ind)
        
        mk.type <- "A.H.B"
        # Code for OneMap
        idx <- which(ind_geno_tot[,2] == "(1, 1) <br>")
        geno_one[idx] <- 2
        idx <- which(ind_geno_tot[,2] == "(2, 0) <br>")
        geno_one[idx] <- 1
        idx <- which(ind_geno_tot[,2] == "(0, 2) <br>")
        geno_one[idx] <- 3
      }
    }
  }
  
  return( list(mk= supermassa_4parallel[[1]],rm.mk, mk.type, geno_one, probs, pgeno))
}

## function to parse the genotype probabilities
parse.geno<-function(x)
{
  y<-strsplit(as.character(x), split="\\[|\\,|\\]")
  sapply(y, function(x) as.numeric(x[which.max(nchar(x))]))
}

recode_parents_sm <- function(x) {
  if(anyNA(x)) {
    "./."
  } else if(x[1] == 0 & x[2] == 2){
    "1/1"
  } else if(x[1] == 2 & x[2] == 0){
    "0/0"
  } else if(x[1] == 1 & x[2] == 1){
    "0/1"
  }
}

# Contains supermassa_error, depth_prepare and supermassa_parallel

##' Runs SuperMassa for multiple markers using doParallel packages
##' Updates OneMap object in genotypes, error probabilitie and marker type
##' Also removes non-informative markers according to SuperMassa genotyping.
##' @param vcf path and name of the vcf file
##' @param vcf.par Field of VCF that informs the depth of alleles
#' @param out_vcf path and name of the vcf file to be outputed from the generated onemap object. 
#' It is important to notice that only onemap informative markers are kept.
##' @param parent1 parent 1 identification in vcfR object
##' @param parent2 parent 2 identification in vcfR objetc
##' @param f1 f1 individual identification if f2 cross type
##' @param crosstype string defining the cross type, by now it supports only 
##' outcross and f2 intercross
##' @param cores number of threads 
##' 
##' @import doParallel parallel
##' 
##' @importFrom matrixStats logSumExp
##' @importFrom stringr str_count
##' @import vcfR 
##' 
##' @export
supermassa_genotype_vcf <- function(vcf=NULL,
                                    vcf.par = c("AD", "DPR"),
                                    out_vcf = NULL,
                                    parent1="P1",
                                    parent2="P2",
                                    crosstype=NULL,
                                    cores = 2,
                                    ploidy = NULL){
  
  if(is.null(ploidy)) stop("Define a ploidy number.")
  if(crosstype != "outcross") stop("Invalid crosstype.")
  
  vcfR.object <- read.vcfR(vcf, verbose = F) 
  input_gt <- extract.gt(vcfR.object)
  
  depths <- extract.gt(vcfR.object, vcf.par)
  parents.id <- which(colnames(depths) %in% c(parent1, parent2))
  
  # Remove missing data
  rm.mks <- which(apply(depths[,parents.id],1, function(x) any(is.na(x))))
  depths <- depths[-rm.mks,]
  oref <- sapply(strsplit(depths, ","), "[[",1)
  oref <- matrix(oref, nrow = nrow(depths))
  oref <- apply(oref, 2, as.numeric)
  osize <- extract.gt(vcfR.object, "DP")
  osize <- osize[-rm.mks,]
  osize <- apply(osize, 2, as.numeric)
  oalt <- osize - oref
  colnames(oref) <- colnames(oalt) <- colnames(depths)
  rownames(oref) <- rownames(oalt) <- rownames(depths)
  
  # Replace NA by 0
  oref[is.na(oref)] <- 0
  oalt[is.na(oalt)] <- 0
  
  depths_prepared <- list()
  for(i in 1:nrow(depths)){
    prog.geno = data.frame(oref[i,-parents.id], oalt[i,-parents.id])
    inds <- colnames(depths)[-parents.id]
    depths_prepared[[i]] <- list(mk=rownames(depths)[i], 
                                 inds= inds,
                                 prog.geno = prog.geno,
                                 paren.geno = data.frame(oref[i,parents.id], oalt[i,parents.id]),
                                 ploidy = ploidy)
  }
  
  cl <- parallel::makeCluster(as.numeric(cores))
  registerDoParallel(cl)
  clusterExport(cl, c("supermassa_parallel_poly"))
  result <- parLapply(cl, depths_prepared, function(x) supermassa_parallel_poly(supermassa_4parallel = x))
  parallel::stopCluster(cl)
  mks <- unlist(lapply(result, "[", 1)) 
  geno <- t(sapply(result, "[[", 2))
  error <- lapply(result, "[[", 3)
  error <- do.call(rbind, error)
  pgeno <-t(sapply(result, "[[", 4))
  
  geno <- cbind(geno, pgeno)

  geno_recode <- recode_geno_vcf(geno, ploidy, oref_cov = oref, osize_cov = osize)

  diffe <- sum(geno_recode != input_gt[-rm.mks,], na.rm = T)/length(geno_recode)
  cat(paste("The approach changed", round(diffe*100,2), "% of the genotypes"))
  
  # set ord mk 1 2 3 ind 1 1 1
  idx <- rep(1:length(unique(error[,2])), length(unique(error[,1])))
  error <- error[order(idx),]
  probs <- as.matrix(apply(error[,-c(1,2)],2,as.numeric))
  
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
  PL <- split(PL, rep(1:length(mks), each = length(colnames(depths))-2))
  PL <- do.call(rbind, PL)
  
  # Parents probs
  zeros <- (ploidy + 1) - str_count(geno_recode[,parents.id], "0")
  template_all <- vector()
  for(i in 1:(ploidy + 1)){
    template <- rep(99, (ploidy + 1))
    template[i] <- "0"
    template_all[i] <- paste0("99:", paste0(template, collapse = ","))
  }
  names(template_all) <- 1:(ploidy + 1)
  prob.parents <- template_all[match(zeros, names(template_all))]
  prob.parents <- matrix(prob.parents, ncol = 2)

  PL <- cbind(PL, prob.parents)
  
  depths_temp <- depths[,-parents.id]
  depths <- cbind(depths_temp, depths[,parents.id])
  
  gt <- matrix(paste0(geno_recode, ":", 
                      depths, ":", 
                      PL),  
               nrow = dim(geno_recode)[1])
  
  
  colnames(gt)  <- colnames(depths)
  rownames(gt)  <- mks
  
  FORMAT <- "GT:AD:GQ:PL"
  
  gt <- cbind(FORMAT, gt)
  
  new.vcfR.object <- vcfR.object
  new.vcfR.object@gt <- gt
  new.vcfR.object@meta[2] <- "##source=Reads2MapTools"
  
  keep <- c(grep("=GT,", new.vcfR.object@meta),
            grep(vcf.par, new.vcfR.object@meta),
            grep("=GQ,", new.vcfR.object@meta),
            grep("=PL,", new.vcfR.object@meta))
  
  new.vcfR.object@meta <- new.vcfR.object@meta[c(1,2, keep)]
  new.vcfR.object@fix <- new.vcfR.object@fix[-rm.mks,]
  new.vcfR.object@fix[,3] <- paste0(new.vcfR.object@fix[,1],"_",new.vcfR.object@fix[,2])
  new.vcfR.object@fix[,"INFO"] <- "."
  
  write.vcf(new.vcfR.object, file =  out_vcf)
}

##' Use output list from depth_prepare and run Supermassa for each marker (first level of the list)
##' Extract from the SuperMassa output, the offspring genotype, the marker type, and the error
##' probabilities
##' @export
supermassa_parallel_poly <- function(supermassa_4parallel){
  ploidy <- supermassa_4parallel[[5]]
  n.ind <- length(supermassa_4parallel[[2]])
  all.ind.names <- supermassa_4parallel[[2]]
  write.table(supermassa_4parallel[[3]], file = paste0("odepth_temp", supermassa_4parallel[[1]],".txt"), quote = FALSE, col.names = FALSE)
  write.table(supermassa_4parallel[[4]], file = paste0("pdepth_temp", supermassa_4parallel[[1]],".txt"), quote = FALSE, col.names = FALSE)
  odepth_temp <- paste0("odepth_temp", supermassa_4parallel[[1]],".txt")
  pdepth_temp <- paste0("pdepth_temp", supermassa_4parallel[[1]],".txt")
  out_file <- paste0("out_prob",supermassa_4parallel[[1]],".txt")
  
  command_mass <- paste("python", paste0(system.file(package = "Reads2MapTools"),"/python_scripts/","SuperMASSA.py"),
                        "--inference f1 --file", paste0(getwd(),"/",odepth_temp), "--ploidy_range ", ploidy,
                        " --f1_parent_data", paste0(getwd(),"/",pdepth_temp), 
                        " --print_genotypes --naive_posterior_reporting_threshold 0",
                        "--save_geno_prob_dist", paste0(getwd(),"/",out_file))
  
  SuperMASSA.output<-system(command_mass, intern=TRUE)
  
  file.remove(odepth_temp,pdepth_temp)
  
  ############ Start code from Molli
  
  ## Assessing estimated ploidy level
  est.ploidy<-as.numeric(strsplit(SuperMASSA.output[4], split = " |,")[[1]][6])
  #cat(" (ploidy ", est.ploidy, ") ", sep = "")
  
  ## Assessing estimated parental dosage
  dpdq <- as.numeric(strsplit(x = SuperMASSA.output[4], split = "\\(|,")[[1]][c(5,8)])
  
  ## Reading SuperMASSA output
  d<-scan(file=out_file, what="character", nlines=1, quiet = TRUE)
  P<-matrix(as.numeric(gsub("[^0-9]", "", unlist(d))), ncol=2, byrow=TRUE)
  A<-read.table(file=out_file, skip=1)
  file.remove(out_file)
  M<-matrix(NA, nrow = length(all.ind.names), ncol = est.ploidy+1,
            dimnames = list(all.ind.names, c(0:est.ploidy)))
  M.temp<-apply(A[,2:(est.ploidy+2)], 2, parse.geno)
  if(is(M.temp, "numeric")) M.temp <- t(as.matrix(M.temp))
  M.temp<-M.temp[,order(P[,2], decreasing=TRUE)]
  if(is(M.temp,"numeric")) M.temp <- t(as.matrix(M.temp))
  dimnames(M.temp)<-list(as.character(A[,1]), c(0:est.ploidy))
  M[rownames(M.temp),]<-M.temp
  mrk1<-apply(M, 1, function(x) exp(x-matrixStats::logSumExp(x)))
  
  probs <- t(mrk1)
  probs <- data.frame(mrk = supermassa_4parallel[[1]],
                      ind = rownames(probs),
                      probs)
  rownames(probs)<-NULL
  colnames(probs)<-c("mrk", "ind", 0:est.ploidy)
  
  pgeno_temp <- strsplit(SuperMASSA.output[grep("ploidy", SuperMASSA.output)], split = ":")
  pgeno_temp <- unlist(strsplit(pgeno_temp[[1]][length(pgeno_temp[[1]])], split = ","))
  pgeno_temp <- pgeno_temp[c(1,3)]
  
  pgeno <- rep(NA, 2)
  for(t in 0:(est.ploidy)){
    pgeno[grep(t, pgeno_temp)] <- t
  }
  
  ind_geno <- unlist(lapply(strsplit(SuperMASSA.output[grep("\t", SuperMASSA.output)], split = "\t"), "[", 1))
  ind_geno <- gsub('.{1}$', '', ind_geno)
  geno <- unlist(lapply(strsplit(SuperMASSA.output[grep("\t", SuperMASSA.output)], split = "\t"), "[", 2))
  
  ind_geno_tot <- cbind(supermassa_4parallel[[2]], rep(NA, n.ind))
  ind_geno_tot[,2][match(ind_geno,ind_geno_tot[,1])] <- geno
  geno_one <- rep(NA, n.ind)
  
  
  geno_temp <- sapply(strsplit(ind_geno_tot[,2], ","), "[[", 1)
  geno_one <- rep(NA, n.ind)
  for(t in 0:(est.ploidy)){
    geno_one[grep(t, geno_temp)] <- t
  }
  
  names(geno_one) <- all.ind.names
  
  return(list(mk= supermassa_4parallel[[1]],geno_one, probs, pgeno))
}
