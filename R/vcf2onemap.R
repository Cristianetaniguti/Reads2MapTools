#' Wrapper function for Scott Brainard java function converter.jar
#'
#' The java script convert the VCF file to onemap raw file
#'
#' @param vcf_file absolute path of the VCF file (can be gzipped)
#' @param data_type cross type (f2 backcross, f2 intercross, ri self, ri sib, or outcross)
#' @param parent1 ID of the parent 1, as encoded in the VCF file
#' @param parent2 ID of the parent 2, as encoded in the VCF file
#' @param outputfile absolute path to use for the OneMap data file and .log files
#' @param verify_uniform_offspring check all aa x aa and aa x bb segregation types, and record violations of
#' expected progeny genotypes (i.e., "aa" and "ab", respectively), in the .log file
#' @param log_filtered_markers false by default. Even when false, some summary statistics about how many markers
#' were filtered and why. If set to true, separate .log files will be generated for markers that are filtered due to:
#' @param only_phased retain only markers that are already phased in the VCF file (true/false, default=false)
#' @param java.converter.path path to converter.jar script
#'
#' @export
vcf2onemap <- function(vcf_file = NULL,
                       data_type=c("f2 backcross", "f2 intercross", "ri self", "ri sib", "outcross"),
                       parent1=NULL,
                       parent2=NULL,
                       outputfile=NULL,
                       verify_uniform_offspring=TRUE,
                       log_filtered_markers = FALSE,
                       only_phased = FALSE,
                       java.converter.path = NULL){

  if(is.null(vcf_file) | is.null(parent1) | is.null(parent2) | is.null(java.converter.path)) stop("Please, define all required arguments\n")

  if(!any(data_type %in% c("f2 backcross", "f2 intercross", "ri self", "ri sib", "outcross"))) stop("This data type is not implemented \n")
  if(verify_uniform_offspring) verify_uniform_offspring="true" else verify_uniform_offspring = "false"
  if(log_filtered_markers) log_filtered_markers="true" else log_filtered_markers="false"
  if(only_phased) only_phased="true" else only_phased="false"

  system(paste0("java -cp ",
                java.converter.path,
                " org.uwm.vcfconverter.Converter female_parent=",
                parent1, " male_parent=",
                parent2,  " vcf_file=",
                vcf_file, " output_file=",
                outputfile, " data_type=",
                data_type, " verify_uniform_offspring=",
                verify_uniform_offspring, " log_filtered_markers=",
                log_filtered_markers))

}
