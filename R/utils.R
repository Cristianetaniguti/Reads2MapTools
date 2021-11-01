
#' Write VCF header random format
#'
#'@export
add_head <- function(vcf, outname, type="simu"){

  vcf_vector <- apply(vcf, 1, function(x) paste(x, collapse = "\t"))
  header1 <- paste0(colnames(vcf), collapse = "\t")

  if(type == "radinitio"){
    header <- paste0("##fileformat=VCFv4.2", "\n",
                     "##source=tskit 0.3.4", "\n",
                     "##FILTER=<ID=PASS,Description=\"All filters passed\">", "\n",
                     "##contig=<ID=1,length=1999993>","\n",
                     "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",'\n',
                     "#", header1)
  } else {
    ##source=RADinitio version 1.1.1 - radinitio.merge_vcf()
    header <- paste0("##fileformat=VCFv4.2", "\n",
                     "##source=RADinitio version 1.1.1 - radinitio.merge_vcf()", "\n",
                     "##FILTER=<ID=PASS,Description=\"All filters passed\">", "\n",
                     "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",'\n',
                     "#", header1)
  }

  vcf <- c(header,vcf_vector)
  write.table(vcf, file = outname, quote = FALSE, row.names = FALSE,  col.names = FALSE)
}


#'  
#' @import onemap
#' 
#' @export
phaseToOPGP_OM <- function(x){
  ## code from here taken from the onemap function print.sequence()
  link.phases <- matrix(NA, length(x$seq.num), 2)
  link.phases[1, ] <- rep(1, 2)
  for (i in 1:length(x$seq.phases)) {
    switch(EXPR = x$seq.phases[i],
           link.phases[i + 1, ] <- link.phases[i, ] * c(1, 1),
           link.phases[i +  1, ] <- link.phases[i, ] * c(1, -1),
           link.phases[i + 1, ] <- link.phases[i, ] * c(-1, 1),
           link.phases[i + 1, ] <- link.phases[i, ] * c(-1, -1))
  }
  if (is(x$data.name, "outcross")) {
    link.phases <- apply(link.phases, 1, function(x) paste(as.character(x), collapse = "."))
    parents <- matrix("", length(x$seq.num), 4)
    for (i in 1:length(x$seq.num)) 
      parents[i, ] <- onemap:::return_geno(x$data.name$segr.type[x$seq.num[i]], link.phases[i])
    ## Our code below
    #transpose the parents and set to baseline
    parents[which(parents == 'a')] <-'A'
    parents[which(parents == 'b')] <- 'B'
    
    parents = t(parents)
    if (length(which(apply(parents[1:2,],2,function(x) !(all(x=='A'))))) > 0)
      if(parents[1,which(apply(parents[1:2,],2,function(x) !(all(x=='A'))))[1]] == 'B')
        parents[1:2,] <- parents[2:1,]
    
    if (length(which(apply(parents[3:4,],2,function(x) !(all(x=='A'))))) > 0)
      if(parents[3,which(apply(parents[3:4,],2,function(x) !(all(x=='A'))))[1]] == 'B')
        parents[3:4,] <- parents[4:3,]
    
    phases <- GUSMap:::parHapToOPGP(parents)
    # If there are multiallelic markers, it returns NULL
    multi <- which(sapply(phases, is.null))
    phases[multi] <- NA
    phases[which(phases == 1 | phases == 4)] <- 17 
    phases[which(phases == 2 | phases == 3)] <- 18
    phases[which(phases == 5 | phases == 8)] <- 19
    phases[which(phases == 6 | phases == 7)] <- 20
    phases[which(phases == 9 | phases == 12)] <- 21
    phases[which(phases == 10 | phases == 11)] <- 22
    phases[which(phases == 13 | phases == 16)] <- 23 
    phases[which(phases == 14 | phases == 15)] <- 24 
    
    ## Now from the parental haplotypes, determine the OPGPs
    return(unlist(phases))
  }
}

#' @export
write_report <- function(tab, out_name, max_cores=1) {
  vroom::vroom_write(tab, paste0(out_name, ".tsv.gz"), num_threads= max_cores)
}

