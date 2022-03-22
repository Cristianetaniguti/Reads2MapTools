#' @export
create_maps_report_simu <- function(input.seq, 
                                    tot_mks,
                                    gab, 
                                    SNPCall, 
                                    GenoCall, 
                                    fake, 
                                    CountsFrom,
                                    real_phases, seed, depth, max_cores) {
  
  # Check genome position
  pos <- as.numeric(input.seq$data.name$POS[input.seq$seq.num])
  sort.pos <- sort(pos)
  if(any(pos != sort.pos)){
    cat("The markers are not ordered by genome position")
    input.seq <- make_seq(input.seq$twopt, input.seq$seq.num[order(as.numeric(input.seq$data.name$POS[input.seq$seq.num]))])
  }
  
  if(fake == "without-false"){
    true_mks <- input.seq$seq.num[which(input.seq$data.name$POS[input.seq$seq.num] %in% tot_mks[,2])]
    seq_true <- make_seq(input.seq$twopt, true_mks) # only true markers are mapped
  } else if(fake == "with-false"){
    real.mks <- input.seq$data.name$POS[input.seq$seq.num] %in% tot_mks[,2]
    real.mks[which(real.mks == T)] <- "true marker"
    real.mks[which(real.mks == F)] <- "false positive"
    seq_true <- input.seq
  }
  
  if(length(seq_true$seq.num) > 100){
    # batch size close to 60 and the overlap is 3/5 of the size (according with paper)
    #div <- round((length(input.seq$seq.num)/60),0)
    #size = round(length(input.seq$seq.num)/div,0)
    #overlap = round(size*(3/5),0)
    #around = 10
    
    batch_size <- pick_batch_sizes(seq_true,
                                   size = 50,
                                   overlap = 30,
                                   around = 10)
    
    map_df <- map_avoid_unlinked(input.seq = seq_true, size = batch_size, 
                                 phase_cores = max_cores, overlap = 30)
    
  } else if(length(seq_true$seq.num) > 2) {
    map_df <- map_avoid_unlinked(seq_true)
  } else {
    map_df <- NA
  }
  
  if(!is.na(map_df)){
    phases <- phaseToOPGP_OM(x = map_df)
    types <- input.seq$data.name$segr.type[map_df[[1]]]
    pos <- input.seq$data.name$POS[map_df[[1]]]
  } else {
    phases <- types <- pos <- NA
  }
  
  if(fake == "without-false"){
    if(!is.na(map_df)){
      real_type <- gab$segr.type[match(as.numeric(map_df$data.name$POS[map_df[[1]]]), as.numeric(gab$POS))]
      real_type[which(is.na(real_type))] <- "non-informative"
      real_phase <- real_phases[which(real_phases[,1] %in% input.seq$data.name$POS[map_df[[1]]]),2]
      poscM <- tot_mks$pos.map[which(as.numeric(as.character(tot_mks$pos)) %in% as.numeric(as.character(pos)))]
      poscM.norm <- c(0,cumsum(diff(poscM)))
      diff= sqrt((poscM.norm - c(0,cumsum(haldane(map_df$seq.rf))))^2)
      
      map_info <- data.frame(seed,
                             depth,
                             "mk.name"= colnames(input.seq$data.name$geno)[map_df[[1]]],
                             "pos" = input.seq$data.name$POS[map_df[[1]]],
                             "rf" = c(0,cumsum(haldane(map_df[[3]]))),
                             "type"= types,
                             "real.type" = real_type,
                             "est.phases"= unlist(phases),
                             "real.phases"= real_phase,
                             "real.mks" = "true marker",
                             "SNPCall" = SNPCall,
                             "GenoCall" = GenoCall,
                             "CountsFrom" = CountsFrom,
                             "fake" = "without-false",
                             "poscM" = poscM,
                             "poscM.norm" = poscM.norm,
                             "diff" = diff)
    } else {
      map_info <- data.frame(seed,
                             depth,
                             "mk.name"= NA,
                             "pos" = NA,
                             "rf" = NA,
                             "type"= NA,
                             "real.type" = NA,
                             "est.phases"= NA,
                             "real.phases"= NA,
                             "real.mks" = NA,
                             "SNPCall" = NA,
                             "GenoCall" = NA,
                             "CountsFrom" = NA,
                             "fake" = "without-false",
                             "poscM" = NA,
                             "poscM.norm" = NA,
                             "diff" = NA)
    }
  } else { # Including fake markers is not possible to do the comparisions
    # The fake markers can also be multiallelic markers
    real.mks <- real.mks[input.seq$seq.num %in% map_df$seq.num]
    map_info <- data.frame(seed,
                           depth,
                           "mk.name"= colnames(input.seq$data.name$geno)[map_df[[1]]],
                           "pos" = input.seq$data.name$POS[map_df[[1]]],
                           "rf" = c(0,cumsum(haldane(map_df[[3]]))),
                           "type"= types,
                           "real.type" = NA,
                           "est.phases"= unlist(phases),
                           "real.phases"= NA,
                           "real.mks" = real.mks,
                           "SNPCall" = SNPCall,
                           "GenoCall" = GenoCall,
                           "CountsFrom" = CountsFrom,
                           "fake" = "with-false",
                           "poscM" = NA,
                           "poscM.norm" = NA, 
                           "diff" = NA)
  }
  
  return(list(map_df, map_info))
}



#' 
#' @import onemap
#' 
#' @export
create_filters_report_simu <- function(onemap_obj, SNPCall, CountsFrom, GenoCall, seed, depth, threshold = NULL) {
  if(!is.null(threshold)){
    onemap_prob <- filter_prob(onemap_obj, threshold = threshold)
  } else onemap_prob <- onemap_obj
  onemap_mis <- filter_missing(onemap_prob, threshold = 0.25)
  bins <- find_bins(onemap_mis)
  onemap_bins <- create_data_bins(onemap_mis, bins)
  segr <- test_segregation(onemap_bins)
  distorted <- select_segreg(segr, distorted = T)
  no_distorted <- select_segreg(segr, distorted = F, numbers = T)
  twopts <- rf_2pts(onemap_bins) # redundant markers are removed
  seq1 <- make_seq(twopts, no_distorted)
  total_variants <- onemap_obj[[3]]
  filters_tab <- data.frame("higher than 25% missing" = onemap_obj$n.mar - onemap_mis$n.mar,
                            "n_markers"= total_variants,
                            "distorted_markers"=length(distorted), # After red removed
                            "redundant_markers"=total_variants - length(bins[[1]]),
                            "after_filters"= length(seq1$seq.num),
                            "SNPCall" = SNPCall,
                            "GenoCall" = GenoCall,
                            "CountsFrom" = CountsFrom, 
                            seed, 
                            depth)
  
  write_report(filters_tab, paste0(SNPCall, "_",GenoCall, "_",CountsFrom, "_", seed,"_",depth,"_filters_report"))
  return(seq1)
}



#' @export
create_gusmap_report_simu <- function(vcf_file, gab, SNPCall, GenoCall, fake, CountsFrom, tot_mks, real_phases, seed, depth){
  ## Maps with gusmap
  RAfile <- VCFtoRA(vcf_file, makePed = T)
  filelist = list.files(pattern = ".*_ped.csv")
  
  ped.file <- read.csv(filelist)
  ID <- c(1:(dim(ped.file)[1]))
  idx.P1 <- which(as.character(ped.file$SampleID) == "P1")
  idx.P2 <- which(as.character(ped.file$SampleID) == "P2")
  
  mother <- c(rep(idx.P1,dim(ped.file)[1]))
  father <- c(rep(idx.P2,dim(ped.file)[1]))
  fam <- c(rep("F1", dim(ped.file)[1]))
  idx.P1 <- which(as.character(ped.file$SampleID) == "P1")
  idx.P2 <- which(as.character(ped.file$SampleID) == "P2")
  mother[c(idx.P1, idx.P2)] <- ""
  father[c(idx.P1, idx.P2)] <- ""
  fam[c(idx.P1, idx.P2)] <- c("", "")
  ped.file$IndividualID <- ID
  ped.file$Mother <- father # Inverted
  ped.file$Father <- mother
  ped.file$Family <- fam
  
  write.csv(ped.file, file = "ped.file.csv")
  filelist = list.files(pattern = ".*.ra.tab")
  RAdata <- readRA(filelist, pedfile = "ped.file.csv", 
                   filter = list(MAF=0.05, MISS=0.25, BIN=0, DEPTH=0, PVALUE=0.05), sampthres = 0)
  
  mydata <- makeFS(RAobj = RAdata, pedfile = "ped.file.csv", 
                   filter = list(MAF = 0.05, MISS = 0.25,
                                 BIN = 1, DEPTH = 0, PVALUE = 0.05, MAXDEPTH=1000))
  
  # Suggested in vignette
  #mydata$rf_2pt(nClust=1)
  #mydata$createLG()
  #mydata$computeMap(mapped = FALSE, inferOPGP = F) Error
  
  # Alternative way
  if(fake == "without-false"){
    keep.mks <- which(mydata$.__enclos_env__$private$pos %in% tot_mks$pos)
    pos <- mydata$.__enclos_env__$private$pos[keep.mks]
    depth_Ref_m <- mydata$.__enclos_env__$private$ref[[1]][,keep.mks]
    depth_Alt_m <- mydata$.__enclos_env__$private$alt[[1]][,keep.mks]
    config <- mydata$.__enclos_env__$private$config[[1]]
    # If there are inferred config
    infer <- which(is.na(config))
    if(length(infer) > 0)
      config[infer] <- mydata$.__enclos_env__$private$config_infer[[1]][infer]
    config <- config[keep.mks]
  } else {
    pos <- mydata$.__enclos_env__$private$pos
    depth_Ref_m <- mydata$.__enclos_env__$private$ref[[1]]
    depth_Alt_m <- mydata$.__enclos_env__$private$alt[[1]]
    config <- mydata$.__enclos_env__$private$config[[1]]
    # If there are inferred config
    infer <- which(is.na(config))
    if(length(infer) > 0)
      config[infer] <- mydata$.__enclos_env__$private$config_infer[[1]][infer]
    real.mks <- mydata$.__enclos_env__$private$pos %in% tot_mks$pos
    real.mks[which(real.mks == T)] <- "true marker"
    real.mks[which(real.mks == F)] <- "false positive"
  }
  
  # Automatically remove markers with pos = inf
  inf.mk <- 0
  rast.pos <- pos
  while(length(inf.mk) > 0){
    depth_Ref <- list(depth_Ref_m)
    depth_Alt <- list(depth_Alt_m)
    
    phases.gus <- GUSMap:::infer_OPGP_FS(depth_Ref_m, depth_Alt_m, 
                                         config, ep=0.001, reltol=1e-3)
    
    rf_est <- GUSMap:::rf_est_FS(init_r = 0.01, ep = 0.001, 
                                 ref = depth_Ref, 
                                 alt = depth_Alt, 
                                 OPGP=list(as.integer(phases.gus)),
                                 nThreads = 1)
    
    # Remove Inf markers
    inf.mk <- which(is.infinite(haldane(rf_est$rf)))
    if(length(inf.mk) > 0){
      depth_Ref_m <- depth_Ref_m[,-inf.mk]
      depth_Alt_m <- depth_Alt_m[,-inf.mk]
      config <- config[-inf.mk]
      rast.pos <- rast.pos[-inf.mk]
    }
  }
  
  dist.gus <- c(0,cumsum(haldane(rf_est$rf))) # haldane mapping function - 
  # I checked with gusmap example 
  # doing as suggested in vignette
  phases.gus[which(phases.gus == 1 | phases.gus == 4)] <- 17
  phases.gus[which(phases.gus == 2 | phases.gus == 3)] <- 18
  phases.gus[which(phases.gus == 5 | phases.gus == 8)] <- 19
  phases.gus[which(phases.gus == 6 | phases.gus == 7)] <- 20
  phases.gus[which(phases.gus == 9 | phases.gus == 12)] <- 21
  phases.gus[which(phases.gus == 10 | phases.gus == 11)] <- 22
  phases.gus[which(phases.gus == 13 | phases.gus == 16)] <- 23
  phases.gus[which(phases.gus == 14 | phases.gus == 15)] <- 24
  real_phase <- real_phases[which(real_phases$pos%in%rast.pos),][,2]
  
  config[which(config==1)] <- "B3.7"
  config[which(config==2 | config==3)] <- "D1.10"
  config[which(config==4 | config==5)] <- "D2.15"
  
  real_type <- gab$segr.type[match(as.numeric(rast.pos), as.numeric(gab$POS))]
  real_type[which(is.na(real_type))] <- "non-informative"
  poscM <- tot_mks$pos.map[which(as.numeric(as.character(tot_mks$pos)) %in% as.numeric(as.character(rast.pos)))]
  poscM.norm <- poscM-poscM[1]
  diff= sqrt((poscM.norm - dist.gus)^2)
  
  if(fake == "without-false"){
    map_info <- data.frame(seed,
                           depth,
                           "mk.name"= mydata$.__enclos_env__$private$SNP_Names[keep.mks][which(pos%in%rast.pos)],
                           "pos" = rast.pos,
                           "rf" = dist.gus,
                           "type"= config,
                           "real.type" = real_type,
                           "est.phases"= phases.gus,
                           "real.phases"= real_phase,
                           "real.mks" = "true marker",
                           "SNPCall" = SNPCall,
                           "GenoCall" = GenoCall,
                           "CountsFrom" = CountsFrom,
                           "fake" = fake,
                           "poscM" = poscM,
                           "poscM.norm" = poscM.norm,
                           "diff" = diff)
  } else {
    map_info <- data.frame(seed,
                           depth,
                           "mk.name"= mydata$.__enclos_env__$private$SNP_Names[which(pos%in%rast.pos)],
                           "pos" = rast.pos,
                           "rf" = dist.gus,
                           "type"= config,
                           "real.type" = NA,
                           "est.phases"= phases.gus,
                           "real.phases"= NA,
                           "real.mks" = real.mks[which(pos%in%rast.pos)],
                           "SNPCall" = SNPCall,
                           "GenoCall" = GenoCall,
                           "CountsFrom" = CountsFrom,
                           "fake" = fake,
                           "poscM" = NA,
                           "poscM.norm" = NA,
                           "diff" = NA)
  }
  
  map_df <- mydata
  return(list(map_df, map_info))
}

#' @export
update_fake_info <- function(info_fake, simu_onemap_obj, ref_alt_alleles, simulated_phases){
  info_correct <- info_fake
  est.pos <- info_fake[[2]]$pos
  real.type <- simu_onemap_obj$segr.type[match(est.pos, simu_onemap_obj$POS)]
  real.type[which(is.na(real.type))] <- "non-informative"
  poscM <- ref_alt_alleles$pos.map[which(as.numeric(as.character(ref_alt_alleles$pos)) %in% as.numeric(as.character(est.pos)))]
  poscM.norm <- poscM-poscM[1]
  diff <- sqrt((poscM.norm - info_fake[[2]]$rf)^2)
  real.phase <- simulated_phases[which(simulated_phases$pos%in%est.pos),][,2]
  
  info_correct[[2]]$real.type <- real.type
  info_correct[[2]]$real.phases <- real.phase
  info_correct[[2]]$fake <- "without-false"
  info_correct[[2]]$poscM <- poscM
  info_correct[[2]]$poscM.norm <- poscM.norm
  info_correct[[2]]$diff <- diff
  
  return(info_correct)
}

