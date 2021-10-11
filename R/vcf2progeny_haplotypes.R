#' Convert phased vcf to onemap_progeny_haplotypes object. Turns possible to use
#' the plot.onemap_progeny_haplotypes (By now only for outcrossing population).
#'
#' @param vcfR.object vcfR object
#' @param ind.id vector of characters with progeny individuals to be evaluated
#' @param group_names string with chromosomes or group names to be evaluated
#' @param parent1 character with parent 1 ID
#' @param parent2 character with parent 2 ID
#' @param crosstype character defining crosstype (outcross, f2 intercross, f2 backcross, ril sib, ril self)
#'
#' @export
vcf2progeny_haplotypes <- function(vcfR.object,
                                   ind.id=NULL,
                                   group_names = NULL,
                                   parent1,
                                   parent2,
                                   crosstype){
  if(is.null(ind.id)){
    stop("You should define one individual.")
  }
  n.mk <- dim(vcfR.object@gt)[1]
  n.ind <- dim(vcfR.object@gt)[2]-1
  INDS <- dimnames(vcfR.object@gt)[[2]][-1]
  
  if(length(which(INDS %in% ind.id)) != length(ind.id)){
    stop("At least one of the individuals in ind.id was not found in vcfR object.")
  }
  
  MKS <- vcfR.object@fix[,3]
  if (any(MKS == "." | is.na(MKS))) MKS <- paste0(vcfR.object@fix[,1],"_", vcfR.object@fix[,2])
  
  # Geno matrix
  GT_matrix <- matrix(rep(NA,n.ind*n.mk), ncol=n.ind, nrow=n.mk)
  GT <- which(strsplit(vcfR.object@gt[1,1], split=":")[[1]]=="GT")
  
  for(i in 2:(n.ind+1))
    GT_matrix[,i-1] <- unlist(lapply(strsplit(vcfR.object@gt[,i], split=":"), "[[", GT))
  
  CHROM <- vcfR.object@fix[,1]
  
  if(is.null(group_names)) group_names <- CHROM[1]
  
  if(length(which(unique(CHROM) %in% group_names)) != length(group_names)){
    stop("At least one of the groups in group_names was not found in vcfR object.")
  }
  
  progeny_haplotypes_obj_chr <- data.frame()
  for(chr in 1:length(group_names)){ ### Need optimization
    CHROM.now <- which(CHROM %in% group_names[chr])
    POS <- as.numeric(vcfR.object@fix[,2])[CHROM.now]
    
    colnames(GT_matrix) <- INDS
    rownames(GT_matrix) <- MKS
    
    P1.idx <- grep(parent1, INDS)
    P2.idx <- grep(parent2, INDS)
    
    P1_1 <- sapply(strsplit(GT_matrix[CHROM.now,P1.idx], "[|]"), "[",1)
    P1_2 <- sapply(strsplit(GT_matrix[CHROM.now,P1.idx], "[|]"), "[",2)
    P2_1 <- sapply(strsplit(GT_matrix[CHROM.now,P2.idx], "[|]"), "[",1)
    P2_2 <- sapply(strsplit(GT_matrix[CHROM.now,P2.idx], "[|]"), "[",2)
    
    progeny_haplotypes_obj_ind <- data.frame()
    for(ind in 1:length(ind.id)){
      ind.idx <- grep(ind.id[ind], INDS)
      ind.number <- grep(ind.id[ind], INDS[-c(P1.idx, P2.idx)])
      ind_1 <- sapply(strsplit(GT_matrix[CHROM.now,ind.idx], "[|]"), "[",1)
      ind_2 <- sapply(strsplit(GT_matrix[CHROM.now,ind.idx], "[|]"), "[",2)
      
      p.names <- c("H1", "H2", "H1", "H2")
      Hs <-  rep(list(rep(NA, length(ind_1))),2)
      progeny_haplotypes_obj <- data.frame()
      for(w in 1:2){ # w is the progeny individual haplotype
        ref.frags<-rep(1, length(ind_1))
        comp <- list(cbind(P1_1, P1_2, P2_1,P2_2,H1=ind_1, H2=ind_2))
        idx.cum <- 1
        while(any(is.na(Hs[[w]]))){
          # count how many equal characters are consecutive
          frags <- rep(list(list(0,0,0,0)),length(comp))
          for(z in 1:length(comp)){
            for(j in 1:4){ # j is the parental haplotype
              if(is.vector(comp[[z]]))
                comp[[z]] <- t(as.matrix(comp[[z]]))
              idx.comp <- comp[[z]][,j] == comp[[z]][,4+w]
              frags[[z]][[j]] <- sequence(rle(as.character(idx.comp))$length)
            }
          }
          # Find the higher fragment in ind1
          new.zs <- list()
          inter.cum <-max(ref.frags)
          for(z in 1:length(frags)){
            max.ind1 <- unlist(lapply(frags[[z]], max))
            # I could't adapt to inbred because it found two possible match
            idx.ind1 <- which.max(max.ind1)
            which.max.ind1 <- which.max(frags[[z]][[idx.ind1]])
            frag <- (which.max.ind1 - max.ind1[idx.ind1]+1):which.max.ind1
            Hs[[w]][ref.frags==idx.cum][frag] <- p.names[idx.ind1]
            # The fragment should be removed and the split the remaining in two
            ref.frags.new <- ref.frags
            if(frag[1] == 1)
              frag1 <- NULL else  frag1 <- comp[[z]][1:(frag[1]-1),]
            if(frag[length(frag)] == dim(comp[[z]])[1])
              frag2 <- NULL else  frag2 <- comp[[z]][(frag[length(frag)]+1):dim(comp[[z]])[1],]
            if(!is.null(frag1)){
              new.zs <- c(new.zs,list(frag1))
              inter.cum <- inter.cum + 1
              ref.frags.new[ref.frags==idx.cum][1:(frag[1]-1)] <- inter.cum
            }
            
            if(!is.null(frag2)){
              new.zs <- c(new.zs,list(frag2))
              inter.cum <- inter.cum + 1
              ref.frags.new[ref.frags==idx.cum][(frag[length(frag)]+1):length(ref.frags.new[ref.frags==idx.cum])] <- inter.cum
            }
            ref.frags <- ref.frags.new
            idx.cum <- idx.cum + 1
          }
          comp <- new.zs
        }
        
        num.mk <- length(CHROM.now)
        df.H <- data.frame(ind = rep(ind.id[ind], 2*num.mk),
                           grp = rep(CHROM[CHROM.now], 2),
                           pos = rep(POS, 2),
                           prob = rep(0, 2*num.mk),
                           parents = rep(c("P1", "P2")[w], 2*num.mk),
                           homologs = rep(c("H1","H2"), each=num.mk))
        
        for(i in 1:length(Hs[[w]])){
          df.H$prob[which(df.H$pos == POS[i] & df.H$homologs == Hs[[w]][i])] <- 1
        }
        # bind homologs
        progeny_haplotypes_obj <- rbind(progeny_haplotypes_obj, df.H)
      }
      # bind individuals
      progeny_haplotypes_obj_ind <- rbind(progeny_haplotypes_obj_ind, progeny_haplotypes_obj)
    }
    # bind chromosomes
    progeny_haplotypes_obj_chr <- rbind(progeny_haplotypes_obj_chr,progeny_haplotypes_obj_ind)
  }
  
  crosstype <- switch(crosstype, "outcross" = "outcross", "f2 intercross"="f2",
                      "f2 backcross"="backcross", "ril sib"="rils", "ril self"="rils")
  
  flag <- "most.likely"
  progeny_haplotypes_obj_chr$marker <- rep(MKS,4)
  progeny_haplotypes_obj_chr <- progeny_haplotypes_obj_chr[, c("ind", "marker","grp", "pos", "prob", "parents", "homologs")]
  progeny_haplotypes_obj_chr$allele <- paste0(progeny_haplotypes_obj_chr$parents, "_",progeny_haplotypes_obj_chr$homologs)
  colnames(progeny_haplotypes_obj_chr)[7] <- "parents.homologs"
  
  class(progeny_haplotypes_obj_chr) <- c("onemap_progeny_haplotypes", crosstype, "data.frame", flag)
  return(progeny_haplotypes_obj_chr)
}
