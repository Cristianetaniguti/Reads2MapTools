
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
