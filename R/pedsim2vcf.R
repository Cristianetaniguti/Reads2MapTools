globalVariables(c("read.table", "rnbinom", "rbinom"))

#' Function to simulate depths and to convert PedigreeSim file to vcf file
#'
#' Given PedigreeSim .dat .map and .chrom files generate vcf file with depths count
#' in AD format field estimated by negative binomial. The function receives genotypes codification in
#' genotypes.dat according with Wu et. al 2002a table I observed bands. The a allele will be considered
#' the reference allele. Null alleles are not supported.
#'
#' @param inputfile file .dat output from PedigreeSim software
#' @param map.file file .map input in PedigreeSim software
#' @param chrom.file file.chrom input in PedigreeSim software
#' @param out.file path to vcf output file
#' @param mean.depth mean of the negative binomial distribution to generate depth counts
#' @param disper.par dispersion parameter for negative binomial distribution
#' @param mean.phred Sequencing error parameter
#' @param chr.mb Chromosome size in mega base.
#' @param method Choose negative binomial ("neg.binom"), poisson ("poisson") distributions or updog ("updog") model to simulate counts values
#' @param miss.perc Percentage of missing data
#' @param pos vector with the position of each marker. Use "cm" if you want to use mapfile position information.
#' @param chr Chromosome where the marker is positioned
#' @param phase if TRUE the genotypes in VCF will be phased
#' @param bias The bias parameter for updog model. Pr(a read after selected) / Pr(A read after selected).
#' @param od The overdispersion parameter for updog model. See the Details of the rho variable in betabinom.
#' @param disper.par Dispertion parameter for negative binomial
#' @param reference.alleles vector defining the reference alleles for each marker
#' @param haplo.ref character indicating the reference haplotype of genotypes.dat file
#' @param use.as.alleles if \code{TRUE} uses codification in genotypes dat to define the reference and alternative
#' alleles fields in VCF, arguments reference.alleles or haplo.ref will be used to define which are the reference alleles
#' @param counts If \code{TRUE} also simulates allele counts using approach defined in \code{method}
#' @param p.mean.depth mean of the negative binomial distribution to generate depth counts for parents
#' @param segregation.distortion.freq numeric between 0 and 1 to define the frequency of distorted markers. For example, if 0.3, the distortion will be applied to 30\% of the markers. This does not consider linkage (see run_pedsim).
#' @param segregation.distortion.mean numeric defining the mean p-value expected in the chi-square test of distorted markers. A normal distribution is used to sample values using the defined mean and standard deviation. This does not consider linkage (see run_pedsim).
#' @param segregation.distortion.sd numeric defining the standard deviation p-value expected in the chi-square test of distorted markers. A normal distribution is used to sample values using the defined mean and standard deviation. This does not consider linkage (see run_pedsim).
#' @param segregation.distortion.seed define seed to set the sample procedures during segregation distortion simulation. This does not consider linkage (see run_pedsim).
#' @param n_selected_loci number of selected loci
#' @param selection_str_mean selection mean intensity 
#' @param selection_str_var selection intensity variance
#' @param pop.size population size
#' 
#' @return vcf file located in out.file defined path
#'
#' @seealso vcf file description <http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/>
#'
#' @author Cristiane Taniguti, \email{chtaniguti@@tamu.edu}
#'
#' @importFrom updog rflexdog
#'
#' @references
#' Voorrips, R. E. , Maliepaard, C. A. (2012) The simulation of meiosis in
#' diploid and tetraploid organisms using various genetic models.
#' \emph{BMC Bioinformatics} 13  - ISSN 1471-2105 - 12 p.
#'
#' Bilton,T. P. , Schofield,  M. R. , Black,  M. A. , Chagné, D. , Wilcox, P. L.
#' and Dodds, K. G. (2018) Accounting for Errors in Low Coverage High-Throughput
#' Sequencing Data When Constructing Genetic Maps Using Biparental Outcrossed Populations.
#' \emph{Genetics}. ISSN 0016-6731.
#' 
#' @export
pedsim2vcf <- function(inputfile=NULL,
                       map.file=NULL,
                       chrom.file=NULL,
                       out.file="out.pedsim2vcf.txt",
                       miss.perc = 0,
                       counts=TRUE,
                       mean.depth=20,
                       p.mean.depth = 20,
                       chr.mb = 10,
                       method = c("updog", "neg.binom"),
                       mean.phred=20,
                       bias=1,
                       od=0.001,
                       disper.par=2,
                       pos=NULL,
                       chr=NULL,
                       phase = FALSE,
                       haplo.ref=NULL,
                       reference.alleles = NULL,
                       use.as.alleles=FALSE,
                       segregation.distortion.freq = NULL,
                       segregation.distortion.mean = NULL,
                       segregation.distortion.sd = NULL,
                       segregation.distortion.seed = 8181,
                       n_selected_loci = NULL, # number of loci under selection
                       selection_str_mean = NULL, # Presence of one allele compared to the other
                       selection_str_var = NULL,
                       pop.size = NULL,
                       selected_mks = NULL,
                       map.size = NULL){
  
  # Do the checks here
  if(is.null(inputfile) | is.null(map.file) | is.null(chrom.file))
    stop("You must define the PedigreeSim output files genotypes, map.file and chrom.file\n")
  
  data <- read.table(paste(inputfile), stringsAsFactors = FALSE, header = TRUE)
  
  if(is.null(chr)){
    chr.info <- read.table(map.file, header = TRUE, stringsAsFactors = FALSE)
    chr <- chr.info$chromosome
  }
  
  if(is.null(pos)){
    chr.info <- read.table(map.file, header = TRUE, stringsAsFactors = FALSE)
    pos.info <- read.table(chrom.file, header = TRUE, stringsAsFactors = FALSE)
    pos <- chr.info$position*((chr.mb*1000)/mean(pos.info$length))
    pos <- round(pos,0)
  } else if(any(pos=="cM")){
    chr.info <- read.table(map.file, header = TRUE, stringsAsFactors = FALSE)
    pos <- round(chr.info$position,2)
  }
  
  if(!is.null(map.size)) {
    chr.info <- read.table(map.file, header = TRUE, stringsAsFactors = FALSE)
    keep.mks <- which(chr.info$position <= map.size)
    chr <- chr[keep.mks]
    pos <- pos[keep.mks]
    data <- data[keep.mks,]
  }
  
  # Infos
  rownames(data) <- data[,1]
  data <- data[,-1]
  n.ind <- dim(data)[2]/2
  n.mk <- dim(data)[1]
  
  # Genotypes matrix
  idx <- rep(1:(n.ind), each = 2)
  gt_matrix <- gt_ref <- gt_alt <-  het_matrix <- matrix(rep(NA,n.ind*n.mk), nrow = n.mk, ncol = n.ind)
  
  data <- as.matrix(data)
  keep.alleles <- data
  
  if(any(data == "o"))
    stop("Null alleles are not supported.")
  
  if(is.null(haplo.ref)){
    h.ref <- data[,1]
  } else{
    h.ref <- data[,which(colnames(data) == haplo.ref)]
  }
  
  alt <- list()
  for(i in 1:n.mk){
    alt <- levels(factor(unlist(data[i,])))[which(levels(factor(unlist(data[i,]))) != h.ref[i])]
    data[i,][which(data[i,]== h.ref[i])] <- "0"
    for(w in 1:length(alt)){
      data[i,][which(data[i,]== alt[w])] <- w
    }
  }
  
  for(i in 1:(length(idx)/2)){
    gt_matrix[,i] <- paste0(data[,which(idx == i)[1]], "|", data[,which(idx == i)[2]])
    het_matrix[,i] <- data[,which(idx == i)[1]] != data[,which(idx == i)[2]]
    gt_ref[,i] <- data[,which(idx == i)[1]]
    gt_alt[,i] <- data[,which(idx == i)[2]]
  }
  
  if(!phase){
    gt_matrix <- gsub("[|]", "/", gt_matrix)
    gt_matrix[which(gt_matrix == "1/0")] <- "0/1"
    gt_matrix[which(gt_matrix == "2/0")] <- "0/2"
    gt_matrix[which(gt_matrix == "3/0")] <- "0/3"
    gt_matrix[which(gt_matrix == "2/1")] <- "1/2"
    gt_matrix[which(gt_matrix == "3/1")] <- "1/3"
    gt_matrix[which(gt_matrix == "3/2")] <- "2/3"
  }
  
  # Add segregation distortion by selection
  if(!is.null(n_selected_loci) & !is.null(selection_str_mean) & !is.null(selection_str_var & !is.null(pop.size))){
    gt_matrix_pro <- gt_matrix[,-c(1:2)]
    if(is.null(selected_mks))
      selected_mks <- sample(nrow(gt_matrix_pro), n_selected_loci)
    
    for(i in 1:length(selected_mks)){
      temp_geno <- gt_matrix_pro[selected_mks[i],]
      n_geno <- table(temp_geno)
      # select individuals to be removed
      geno_one <- which(temp_geno == names(n_geno)[1])
      rm_ind <- round(rnorm(1, selection_str_mean, selection_str_var)*length(geno_one),0)
      gt_matrix_pro <- gt_matrix_pro[,-geno_one[sample(length(geno_one), rm_ind)]]
    }
    
    gt_matrix_pro <- gt_matrix_pro[,sample(ncol(gt_matrix_pro), pop.size)]
    gt_matrix <- cbind(gt_matrix[,1:2], gt_matrix_pro)
    
    n.ind <- dim(gt_matrix)[2]
    n.mk <- dim(gt_matrix)[1]
  }
  
  # Add segregation distortion by genotyping error
  if(!(is.null(segregation.distortion.freq) & is.null(segregation.distortion.mean) & is.null(segregation.distortion.sd))){
    p_values <- abs(rnorm(round(nrow(gt_matrix)*segregation.distortion.freq,0), 
                          mean = segregation.distortion.mean, 
                          sd = segregation.distortion.sd))
    
    # Three possible segregation pattern: only biallelic codominant markers (B3.7, D1.10 and D2.15)
    # Freedom degrees B3.7 - 2; D1.10 - 1; D2.15 - 1
    b3 <- which(gt_matrix[,1] == "0/1" & gt_matrix[,2] == "0/1" | 
                  gt_matrix[,1] == "0|1" & gt_matrix[,2] == "0|1" | 
                  gt_matrix[,1] == "1|0" & gt_matrix[,2] == "0|1" | 
                  gt_matrix[,1] == "0|1" & gt_matrix[,2] == "1|0")
    fd <- rep(1, nrow(gt_matrix))
    fd[b3] <- 3 # phased genotypes, the heterozygous have two categories
    
    set.seed(segregation.distortion.seed)
    mk.segr <- sample(1:nrow(gt_matrix), round(nrow(gt_matrix)*segregation.distortion.freq,0))
    mk.segr.names <- rownames(data)[mk.segr]
    chi <- round(qchisq(p_values, fd[mk.segr], lower.tail = F),0) 
    expected <- ncol(gt_matrix)/(fd[mk.segr]+1)
    n.change <- sqrt(chi*expected) # number of genotypes that need to be changed
    
    # This change the genotypes line by line
    counts.mk <- counts.mk.change <- list()
    gt_matrix_change <- gt_matrix
    for(i in 1:length(mk.segr)){
      mk <- gt_matrix[mk.segr[i],3:ncol(gt_matrix)]
      mk_change <- mk
      
      counts.mk[[i]] <- table(mk)
      genos.mk <- unique(mk)
      if(length(genos.mk) == 1) { # for, non-informative markers
        replace.genos <- c("0|0", "1|1")
        which.geno <- 1
      } else if(length(genos.mk) == 2) {
        which.geno <- sample(1:length(genos.mk), 1) # Sample genotype that will be changed
        replace.genos <- genos.mk[-which.geno]
      } else {
        which.geno <- which(genos.mk == "0|1" | genos.mk == "1|0")
        replace.genos <- c("0|0", "1|1")
      }
      # Sample individuals genotypes and to each genotype will be changed
      mk_change[which(mk %in% genos.mk[which.geno])][sample(1:length(which(mk %in% genos.mk[which.geno])), n.change[i])] <- sample(replace.genos, n.change[i], replace = T)
      gt_matrix_change[mk.segr[i],3:ncol(gt_matrix)] <- mk_change
      counts.mk.change[[i]] <- table(mk_change)
    }
    
    gt_matrix <- gt_matrix_change
  }
  
  ## Simulate counts
  if(counts==TRUE){
    if(method=="neg.binom" ){
      # Negative binomial to estimate the depths (code adaptaded from Gusmap)
      depth <- prob.mat <- matrix(rep(NA, n.ind*n.mk),nrow=n.mk, ncol=n.ind)
      
      prob.mat[which(gt_matrix == "0/0" | gt_matrix == "0|0")] <- 1
      prob.mat[het_matrix] <- 0.5
      prob.mat[is.na(prob.mat)] <- (10^(-mean.phred/10))
      
      depth[which(!is.na(gt_matrix))] <- rnbinom(sum(!is.na(gt_matrix)),mu=mean.depth,size = disper.par)
      
      # Avoiding missing
      idx <- which(depth==0)
      
      while(length(idx) >0){
        depth[which(depth==0)] <- rnbinom(length(idx),mu=mean.depth,size=disper.par)
        idx <- which(depth==0)
      }
      ref_matrix <- matrix(rbinom(n.ind*n.mk, depth, prob.mat), nrow = n.mk)
      
      
      alt_matrix <- depth-ref_matrix
      
      info <- paste0("DP=", apply(depth,1,sum))
      
    } else if(method=="updog"){
      
      up_matrix <- size_matrix <- ref_matrix <- matrix(rep(NA, length(gt_matrix)), nrow = dim(gt_matrix)[1])
      
      idx <- which(het_matrix)
      up_matrix[idx] <- 1
      idx <- which(gt_matrix == "0/0" | gt_matrix == "0|0")
      up_matrix[idx] <- 0
      idx <- which(is.na(up_matrix))
      up_matrix[idx] <- 2
      
      size_matrix <- matrix(rnbinom(length(gt_matrix),mu=mean.depth,size=disper.par), nrow = dim(gt_matrix)[1])
      
      # Parents with other depth
      if(!is.null(p.mean.depth)){
        size_matrix[,1:2] <- rnbinom(length(size_matrix[,1:2]),mu=p.mean.depth,size=disper.par)
      }
      
      mis <- which(size_matrix==0)
      
      while(length(mis) > 0){
        size_matrix[mis] <- rnbinom(length(mis),mu=mean.depth,size=disper.par)
        mis <- which(size_matrix==0)
      }
      
      for(i in 1:dim(up_matrix)[1]){
        ref_matrix[i,] <- rflexdog(sizevec = size_matrix[i,],
                                   geno=up_matrix[i,],
                                   ploidy = 2,
                                   seq=10^(-mean.phred/10),
                                   bias=bias,
                                   od = od)
      }
      
      alt_matrix <- size_matrix-ref_matrix
      
      info <- apply(size_matrix, 1, sum)
    }
    
    # VCF format field
    format <- rep("GT:AD", n.mk)
    
    # Select the number of less frequent allele
    minor_matrix <- matrix(rep(NA, n.mk*n.ind), nrow = n.mk, ncol = n.ind)
    minor_matrix[which(ref_matrix >= alt_matrix)] <- alt_matrix[which(ref_matrix >= alt_matrix)]
    minor_matrix[which(ref_matrix < alt_matrix)] <- ref_matrix[which(ref_matrix < alt_matrix)]
    
    # Probabilities calculation by binomial distribution
    tot_matrix <- ref_matrix + alt_matrix
    tot_matrix[which(tot_matrix==0)] <- NA
    het_matrix <- choose(tot_matrix, minor_matrix)*(0.5^minor_matrix)*(0.5^(tot_matrix-minor_matrix))
    homo_matrix <- choose(tot_matrix, minor_matrix)*((10^(-mean.phred/10))^minor_matrix)*((1-(10^((-mean.phred/10))))^(tot_matrix-minor_matrix))
    homo.oalt <- choose(tot_matrix, minor_matrix)*((10^(-mean.phred/10))^(tot_matrix-minor_matrix))*((1-(10^((-mean.phred/10))))^minor_matrix)
    
    # Reviewed matrix
    check_matrix <- matrix(rep(NA, n.mk*n.ind), nrow = n.mk, ncol = n.ind)
    
    # search in initial file the alleles from heterozygotes
    idx <- which(het_matrix >= homo_matrix | het_matrix == Inf)
    check_matrix[idx][which(gt_ref[idx] != gt_alt[idx])] <- gt_matrix[idx][which(gt_ref[idx] != gt_alt[idx])]
    
    if(phase){
      check_matrix[idx][which(gt_ref[idx] == gt_alt[idx])] <- paste0("0|",sapply(strsplit(gt_matrix[idx][which(gt_ref[idx] == gt_alt[idx])], "[|]"), unique))
    } else {
      check_matrix[idx][which(gt_ref[idx] == gt_alt[idx])] <- paste0("0/",sapply(strsplit(gt_matrix[idx][which(gt_ref[idx] == gt_alt[idx])], "/"), unique))
    }
    
    # search in initial file the alleles from homozigotes
    idx <- which(het_matrix < homo_matrix)
    check_matrix[idx][which(gt_ref[idx] == gt_alt[idx])] <- gt_matrix[idx][which(gt_ref[idx] == gt_alt[idx])]
    
    if(phase){
      allele <- sapply(strsplit(gt_matrix[idx][which(gt_ref[idx] != gt_alt[idx])], "[|]"), "[", 1)
      check_matrix[idx][which(gt_ref[idx] != gt_alt[idx])] <- paste0(allele, "|", allele)
    } else {
      allele <- sapply(strsplit(gt_matrix[idx][which(gt_ref[idx] != gt_alt[idx])], "/"), "[", 1)
      check_matrix[idx][which(gt_ref[idx] != gt_alt[idx])] <- paste0(allele, "/", allele)
    }
    
    chang <- table((gt_matrix == check_matrix))
    
    if(dim(chang) ==2){
      cat("Counts simulation changed", (chang[1]/(chang[1] + chang[2]))*100,
          "% of the given genotypes\n")
    } else if(names(chang) == "FALSE"){
      cat("All genotypes were changed")
    } else{
      cat("None genotypes were changed")
    }
    
    # Reference allele is the most frequent
    # for(i in 1:dim(check_matrix)[[1]]){
    #   z <- table(check_matrix[i,])
    #   w <- strsplit(ad_matrix[i,], split = ",")
    #   if(length(z) == 2){
    #     if(all(names(z) == c("0/1", "1/1"))){
    #       idx1 <- which(check_matrix[i,] == "1/1")
    #       idx2 <- which(check_matrix[i,] == "0/0")
    #       check_matrix[i,][idx1] <- "0/0"
    #       w[idx1] <- lapply(w[idx1], "rev")
    #       check_matrix[i,][idx2] <- "1/1"
    #       w[idx2] <- lapply(w[idx2], "rev")
    #       ad_matrix[i,] <- unlist(lapply(w, function(x) paste0(x, collapse = ",")))
    #     }
    #   } else if(length(z) == 3){
    #     if(all(names(z) == c("0/0", "0/1", "1/1"))){
    #       if(z[[1]] < z[[3]]){
    #         idx1 <- which(check_matrix[i,] == "1/1")
    #         idx2 <- which(check_matrix[i,] == "0/0")
    #         check_matrix[i,][idx1] <- "0/0"
    #         w[idx1] <- lapply(w[idx1], "rev")
    #         check_matrix[i,][idx2] <- "1/1"
    #         w[idx2] <- lapply(w[idx2], "rev")
    #         ad_matrix[i,] <- unlist(lapply(w, function(x) paste0(x, collapse = ",")))
    #       }
    #     }
    #   }
    # }
    
    
    ad_matrix <- matrix(NA, nrow = nrow(check_matrix), ncol = ncol(check_matrix))
    for(j in 1:dim(check_matrix)[1]){
      for(i in c(3,2,1)){
        if(any(grepl(i, check_matrix[j,]))){
          ncol <- i + 1
          break
        }
      }
      ad_matrix_temp <- matrix(0, nrow= dim(check_matrix)[2], ncol = ncol)
      for(i in 0:(ncol-1)){
        idx <- which(gt_ref[j,] == i & gt_alt[j,] == i)
        idx.sub <- which(ref_matrix[j,idx] != 0)
        ad_matrix_temp[idx[idx.sub],i+1] <- ref_matrix[j,idx[idx.sub]]
        idx.sub <- which(alt_matrix[j,idx] != 0)
        ad_matrix_temp[idx[idx.sub],i+1] <- alt_matrix[j,idx[idx.sub]]
        
        idx <- which(gt_ref[j,] == i & gt_alt[j,] != i)
        idx.sub <- gt_ref[j,][idx] == i
        ad_matrix_temp[idx[idx.sub],i+1] <- ref_matrix[j,idx[idx.sub]]
        idx.sub <- gt_alt[j,][idx] == i
        ad_matrix_temp[idx[idx.sub],i+1] <- alt_matrix[j,idx[idx.sub]]
        
        idx <- which(gt_ref[j,] != i & gt_alt[j,] == i)
        idx.sub <- gt_ref[j,][idx] == i
        ad_matrix_temp[idx[idx.sub],i+1] <- ref_matrix[j,idx[idx.sub]]
        idx.sub <- gt_alt[j,][idx] == i
        ad_matrix_temp[idx[idx.sub],i+1] <- alt_matrix[j,idx[idx.sub]]
      }
      ad_matrix[j,] <- apply(ad_matrix_temp, 1, function(x) paste0(x, collapse = ","))
    }
    
    vcf_format <- matrix(paste0(check_matrix, ":", ad_matrix), nrow = n.mk)
    
  } else {
    vcf_format <- check_matrix <- gt_matrix
    format <- rep("GT", n.mk)
    info <- "."
  }
  
  # Adding missing data
  miss <- sample(1:length(vcf_format), length(vcf_format)*(miss.perc/100))
  if(length(miss)>0){
    if(counts){
      vcf_format[miss] <- "./.:0,0"
    } else{ vcf_format[miss] <- "./." }
  }
  
  colnames(vcf_format) <- c("P1", "P2", paste0("F1_", formatC(1:(dim(gt_matrix)[2]-2), width = nchar(dim(gt_matrix)[2]-2), format = "d", flag = "0")))
  
  id <- rownames(data)
  
  # Defining alleles in field REF and ALT
  ## REF
  if(use.as.alleles){
    if(!is.null(reference.alleles)){
      ref <- reference.alleles 
      if(!is.null(map.size)) ref <- ref[keep.mks]
    } else { 
      ref <- h.ref 
    }
  } else{
    ref <- sample(c("A","T", "C", "G"), n.mk, replace = TRUE)
  }
  alt <- rep(NA, length(ref))
  ## ALT
  done <- vector() # vector to store markers already evaluated
  guide <- 1:dim(check_matrix)[1]
  for(j in c(3,2,1)){
    if(length(done) != 0){
      alleles <- apply(check_matrix[-done,], 1, function(x) any(grepl(j, x)))
    } else {
      alleles <- apply(check_matrix, 1, function(x) any(grepl(j, x)))
    }
    if(sum(alleles) != 0){
      alt_temp <- vector()
      if(length(done) != 0)
        ref_temp <- ref[-done][alleles]
      else ref_temp <- ref[alleles]
      
      for(i in 1:length(alleles)){
        if(use.as.alleles){
          if(length(done) == 0)
            temp <- levels(factor(keep.alleles[which(alleles)[i],]))
          else temp <- levels(factor(keep.alleles[-done,][which(alleles)[i],]))
        } else  temp <- sample(c("A","T", "C", "G"), 4, replace = F)
        alt_temp <- rbind(alt_temp, temp[-which(temp == ref_temp[i])])
      }
      rm.colum <- 3-j
      if(rm.colum != 0 & !use.as.alleles){
        alt_temp <- alt_temp[,-c(1:rm.colum)]
      }
      if(length(done) != 0){
        if(is(alt_temp, "matrix"))
          alt[-done][alleles] <- apply(alt_temp, 1, function(x) paste0(x, collapse = ","))
        else alt[-done][alleles] <- alt_temp
        done <- c(done,guide[-done][alleles])
      } else {
        if(is(alt_temp, "matrix"))
          alt[alleles] <- apply(alt_temp, 1, function(x) paste0(x, collapse = ","))
        else alt[alleles] <- alt_temp
        done <- c(done,guide[alleles])
      }
    }
  }
  
  qual <- rep(".", n.mk)
  filter <- "PASS"
  # transformar em vetor
  vcf_file_mks <- data.frame("CHROM"=chr, "POS"=pos, "ID"= id,"REF"=ref, "ALT"=alt,
                             "QUAL"=qual, "FILTER"=filter,"INFO"=info,"FORMAT"=format,vcf_format, stringsAsFactors = FALSE)
  
  vcf_vector <- apply(vcf_file_mks, 1, function(x) paste(x, collapse = "\t"))
  # Remove empty space before position
  vcf_vector <- gsub(pattern = " ", replacement = "",x = vcf_vector)
  
  header1 <- paste0(colnames(vcf_file_mks), collapse = "\t")
  header <- paste0("##fileformat=VCFv4.1", "\n", "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", '\n',
                   "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">",'\n',
                   "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", "\n", "#", header1)
  
  vcf <- c(header,vcf_vector)
  write.table(vcf, file = paste(out.file), quote = FALSE, row.names = FALSE,  col.names = FALSE)
}
