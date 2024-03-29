#' Wrapper function for PedigreeSim software
#'
#' With the function arguments users can define the basic characteristics
#' of the maps to be simulated or users can provide a VCF and a reference linkage map files
#' to define the number and type of markers, and genetic distances.
#' The PedigreeSim software offers other options that are not available here.
#'
#' @param chromosome a vector of characters defining the chromosomes names
#' @param n.marker a vector of integers defining the number of markers within chromosome
#' @param tot.size.cm a vector of double numbers defining the total size in centimorgan of each chromosome
#' @param ref.vcf object of class vcfR
#' @param ref.map data.frame with position in base pair and in centimorgan (reference linkage map)
#' @param centromere a vector of double numbers defining the centromere position
#' @param n.ind a integer defining the number of individuals to be simulated in the progeny
#' @param mk.types a vector of characters defining the marker types according with Wu2002a
#' (the types will be distribuited randomly by the chromosomes)
#' @param seed integer defining the simulation seed
#' @param n.types a vector of integers defining the number of markers of each type
#' @param pop a character defining the population type. Supported: "F1", "F2", "BC"
#' @param path.pedsim path to the .jar PedigreeSim file. If you still don't have it,
#' download in https://www.wur.nl/en/show/Software-PedigreeSim.htm
#' @param name.mapfile a character with the mapfile name
#' @param name.founderfile a character with the founderfile name
#' @param name.chromfile a character with the chromfile name
#' @param name.parfile a character with the parfile name
#' @param name.out a character with the PedigreeSim outputs names
#' @param rm.tempfiles logical TRUE/FALSE, if TRUE remove the temporary files generated
#' @param parent1 parent1 ID defined in reference VCF
#' @param parent2 parent2 ID defined in reference VCF
#'
#' @return PedigreeSim outputs including the population genotypes codified as the observed bands genotypes
#' according to Wu et. al 2002 table I if a VCF file is not included. If a reference VCF is included, the genotypes
#' in output are accordind to the reference and alternative allele defined in VCF.
#'
#' @author Cristiane Taniguti, \email{chtaniguti@@usp.br}
#'
#' @importFrom stringr str_pad
#'
#' @references
#'
#' Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002)
#' Simultaneous maximum likelihood estimation of linkage and linkage phases in
#' outcrossing species. \emph{Theoretical Population Biology} 61: 349-363.
#'
#' @examples
#' \dontrun{
#' # For outcross
#' run_pedsim(chromosome = c("Chr1", "Chr2"), n.marker = c(40,32), tot.size.cm = c(100,150), centromere = c(50, 75),
#'            n.ind = 200, mk.types = c("A1", "A2", "A3", "A4", "B1.5", "B2.6", "B3.7",
#'                                      "C.8", "D1.9", "D1.10", "D1.11", "D1.12", "D1.13",
#'                                      "D2.14", "D2.15", "D2.16", "D2.17", "D2.18"),
#'                                      n.types = rep(4,18), pop = "F1", name.mapfile = "mapfile.map", 
#'                                      name.founderfile="founderfile.gen", name.chromfile="sim.chrom", 
#'                                      name.parfile="sim.par", name.out="sim_out")
#' }
#'
#' @export
run_pedsim <- function(chromosome = c("Chr1", "Chr2", "Chr3"),
                       n.marker = c(18,19,20),
                       tot.size.cm = c(18,19,20),
                       ref.vcf = NULL,
                       ref.map = NULL,
                       centromere = c(9,9,10),
                       n.ind = 150,
                       mk.types = c("A1", "D1.10", "D2.15"),
                       n.types = c(7,25,25),
                       pop = c("F1", "F2", "BC"),
                       seed = 8080,
                       name.mapfile = "mapfile.map",
                       name.founderfile="founderfile.gen",
                       name.chromfile="sim.chrom",
                       name.parfile="sim.par",
                       name.out="sim_out",
                       rm.tempfiles=F,
                       parent1 = NULL,
                       parent2 = NULL){

  # Checks
  if(is.null(ref.map) & is.null(ref.vcf)){
    sizes <-sapply(list(chromosome, n.marker, tot.size.cm, centromere), length)
    if(length(unique(sizes))!=1){
      stop("Sizes of vectors chromosome, n.marker, tot.size.cm and centromere must be equal")
    }

    sizes <-sapply(list(mk.types, n.types), length)
    if(length(unique(sizes))!=1){
      stop("Sizes of vectors mk.types and n.types must be equal")
    }

    if(sum(n.marker) != sum(n.types)) stop("Number of marker in vectors n.marker and n.types need to be equal")

    ## Map file
    # Marker names
    marker1 <- "M"
    marker2 <- 1:sum(n.marker)
    marker2 <- str_pad(marker2,3,pad="0",side = "left")
    marker <-paste0(marker1,marker2)

    # Chromossome and position
    pos <- chr <- vector()
    for(i in 1:length(chromosome)){
      int <- tot.size.cm[i]/n.marker[i]
      pos <- c(pos,seq(from=0, to=tot.size.cm[i], by=int)[-1])
      chr <- c(chr,rep(chromosome[i],n.marker[i]))
    }

    if(length(pos) == length(chr) -1)
      pos <- c(pos,tot.size.cm)

    if(length(pos) == length(chr) +1)
      pos <- pos[-1]

    map_file <- data.frame(marker=marker, chromosome=chr, position= pos)
    write.table(map_file, file = name.mapfile, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

    # Founderfile
    ## 18 marker types (Wu2002a)
    A1 <- c("a","b","c","d")
    A2 <- c("a", "b", "a", "c")
    A3 <- c("a","b","c","o")
    A4 <- c("a","o","b","o")

    B1.5 <- c("a", "b", "a","o")
    B2.6 <- c("a", "o", "a","b")
    B3.7 <- c("a", "b", "a", "b")

    C.8 <- c("a","o","a","o")

    D1.9 <- c("a", "b", "c", "c")
    D1.10 <- c("a","b", "a","a")
    D1.11 <- c("a", "b", "o", "o")
    D1.12 <- c("b", "o", "a", "a")
    D1.13 <- c("a", "o", "o", "o")

    D2.14 <- c("c", "c", "a", "b")
    D2.15 <- c("a", "a", "a", "b")
    D2.16 <- c("o", "o", "a", "b")
    D2.17 <- c("a","a", "b", "o")
    D2.18 <- c("o", "o", "a", "o")

    A.H.B <- c("a", "a", "b", "b")
    C.A <- c("a","a","o","o")
    D.B <- c("o","o", "b","b")

    all.types <- rbind(A1, A2, A3, A4, B1.5, B2.6, B3.7,
                       C.8, D1.9, D1.10, D1.11, D1.12, D1.13,
                       D2.14, D2.15, D2.16, D2.17, D2.18,A.H.B,C.A,D.B)

    sele <- which(rownames(all.types) %in% mk.types)

    if(pop=="F2"){
      if(!all(sele %in% c(19:21))) stop(cat("The marker types ", mk.types, " are not expected for inbred parents in F2 populations\n")) # rows for F2 marker types
    }

    if(pop=="BC"){
      if(sele != 19) stop(cat("The marker types ", mk.types, " are not expected for inbred parents in BC populations\n")) # rows for F2 marker types
    }

    p.geno <- vector()
    for(j in 1:length(sele)){
      vec <- all.types[sele[j],]
      for(i in 1:n.types[j]){
        vec.sort <- c(sample(vec[1:2],2,replace = F),sample(vec[3:4],2,replace = F))
        p.geno <- rbind(p.geno,vec.sort)
      }
    }

    p.geno.s <- p.geno[sample(1:dim(p.geno)[1], dim(p.geno)[1], replace = F ),]

    founder_file <- cbind(marker, p.geno.s)
    colnames(founder_file) <- c("marker", "P1_1", "P1_2", "P2_1", "P2_2")

    write.table(founder_file, file = name.founderfile, quote=FALSE, col.names = TRUE, row.names = FALSE, sep = "\t" )

    ## CHROM file
    sim.chrom <- data.frame(chromosome, length=tot.size.cm, centromere, prefPairing=rep("0.0",length(centromere)),
                            quadrivalents=rep("0.0",length(centromere)))

    write.table(sim.chrom, file = name.chromfile, quote=FALSE, col.names = TRUE, row.names = FALSE, sep = "\t" )

    ## Parameters file
    parameter <- paste0("PLOIDY = 2
                      MAPFUNCTION = HALDANE
                      MISSING = NA
                      CHROMFILE = ", name.chromfile,"
                      POPTYPE = ", pop,"
                      SEED =", seed, "
                      POPSIZE = ", n.ind,"
                      MAPFILE = ", name.mapfile,"
                      FOUNDERFILE = ", name.founderfile,"
                      OUTPUT = ", name.out)

    write.table(parameter, file = name.parfile, quote=FALSE, col.names = FALSE, row.names = FALSE, sep = "\t" )
  } else {
    ref_map <- remove_outlier(ref.map, thr=0) # Remove inverted markers

    founderfile <- create_haplo(vcfR.obj = ref.vcf, ref.map = ref_map, seed = seed,
                                P1 = parent1, P2= parent2, filename = name.founderfile)

    mapfile <- create_mapfile(ref.vcf, ref_map, filename = name.mapfile, rm_outlier = F)

    create_parfile(seed, n.ind, filename = name.parfile, name.mapfile = name.mapfile,
                   name.founderfile = name.founderfile, name.chromfile = name.chromfile,
                   name.out = name.out)

    create_chromfile(mapfile[[1]], filename = name.chromfile)
  }

  # Run pedigreesim
  system(paste("java -jar",system.file(package = "Reads2MapTools", "PedigreeSim/PedigreeSim.jar"), name.parfile))

  
  
  
  
  
  
  if(rm.tempfiles){
    file.remove(name.parfile, paste0(name.out, "_alleledose.dat"),
                paste0(name.out, ".hsa"), paste0(name.out, ".hsb"), paste0(name.out, ".ped"))
  }
  if(!is.null(ref.vcf)) return(mapfile[[2]])
}

