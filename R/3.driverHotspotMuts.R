# Author: J. Alberto Nakauma-González
# Date: 03-08-2023
# e-mail: j.nakaumagonzalez@erasmusmc.nl 
# Script used to find driver genes of hotspot mutations associated to APOBEC activity
#

# Clean everything and load packages and reference genome ---------------------
rm(list=ls())
dev.off()
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load packages
pacman::p_load('plyr', 'dplyr', 'tidyr', 'ggplot2', 'GenomicRanges', 'IRanges', 'ggpubr')


# Load data
load("/APOBECproject/RData/results_hairpinLoopsStability.RData")
ChIPSeqData <- read.csv(file = "data/ChIPSeqData_bladderC.csv")


# Count hotspot mutations
results_hairpinLoopsStability <- results_hairpinLoopsStability %>%
  dplyr::add_count(chr, pos, name = "totalHotspotPos") %>%
  dplyr::add_count(chr, pos, ref, alt, name = "totalSpecificHotspotMut")


# get chromatin regions for each mutation (defined by ChIPSeq experiments of normal tissue; available for bladder and breast cancer in this repository)
results_hairpinLoopsStability <- results_hairpinLoopsStability %>% dplyr::mutate(bin = ceiling(pos / 1E6)) %>%
  dplyr::left_join(ChIPSeqData, by = c('chr','bin')) %>%
  dplyr::filter(!is.na(DNA_Accessibility))

# Define Tri-nucleotide context
results_hairpinLoopsStability <- results_hairpinLoopsStability %>%
  dplyr::mutate(triNucleotideContext = paste0(upstreamNuc, gsub('>.*', '', mutChange), downstreamNuc))


# Only keep unique positions (with the most dominant specific mutations in case of multiple nucleotide changes in the same position)
results_hairpinLoopsStability <- results_hairpinLoopsStability %>%
  dplyr::arrange(-totalSpecificHotspotMut) %>% dplyr::distinct(chr, pos, .keep_all = TRUE)

# get all mutations for later
results_hairpinLoopsStability_all <- results_hairpinLoopsStability



########## Analyse APOBEC mutations ###########

# Keep only TpC hotspot mutations and hairpin loops (also NpC in hairpin loops)
results_hairpinLoopsStability <- results_hairpinLoopsStability %>% dplyr::filter(isTpC == TRUE | formLoop == 1) %>%
  dplyr::mutate(isApC = base::grepl('A\\[C', mutContext)) %>%
  dplyr::mutate(isCpC = base::grepl('C\\[C', mutContext)) %>%
  dplyr::mutate(isGpC = base::grepl('G\\[C', mutContext)) %>%
  dplyr::filter(isTpC == TRUE | isApC == TRUE | isCpC == TRUE | isGpC == TRUE) %>%
  dplyr::mutate(isTpC = ifelse(isTpC == TRUE, "TpC", ifelse(isApC == TRUE, "ApC", ifelse(isCpC == TRUE, "CpC", "GpC")))) %>%
  dplyr::mutate(isTpC = factor(isTpC, levels = c("TpC", "ApC", "CpC", "GpC")))


# get hotspot mutations (totalHotspotPos > 1) outside or within hairpin-loops
results_hairpinLoopsStability_notHairpin <- dplyr::filter(results_hairpinLoopsStability, formLoop == 0, totalHotspotPos > 1)
results_hairpinLoopsStability_Hairpin <- results_hairpinLoopsStability %>% dplyr::filter(formLoop == 1, totalHotspotPos > 1) %>%
  dplyr::mutate(triNucleotideContext = ifelse(isTpC == 'TpC', triNucleotideContext, as.character(isTpC))) %>%
  dplyr::mutate(triNucleotideContext = factor(triNucleotideContext, levels = c("TCA", "TCC", "TCG", "TCT", "ApC", "CpC", "GpC")))
  

########## driver hotspot mutations in non-hairpin loops ########################################################

# scale number of hotspots to zero, for proper application of the poisson distribution 
results_hairpinLoopsStability_notHairpin <- results_hairpinLoopsStability_notHairpin %>%
  dplyr::mutate(totalHotspotPos_2minus = totalHotspotPos - 2)

# Run in parallel modelling for poisson process
summResults_gmaPoissModel_noHairpin <- pbapply::pblapply(base::split(results_hairpinLoopsStability_notHairpin, results_hairpinLoopsStability_notHairpin$triNucleotideContext), function(x){
  
  # initialize a column to fill in p-values
  x$pValue <- 1
  
  # fit a model excluding the focal hotspot mutation
  for(iHotspotMut in c(1:nrow(x))) {

    # get number of mutations of focal hotspot mut and termite program if nMutsFocal == 0
    nMutsFocal <- x[iHotspotMut,]$totalHotspotPos_2minus
    if (nMutsFocal == 0) {next}
    
    # exclude all hotspot muts with higuer frequency than nMutsFocal
    data <- x[-iHotspotMut,]
    data <- data %>% dplyr::filter(totalHotspotPos_2minus <= nMutsFocal)
    
    # fit data to poisson process
    gml_PoissonModel <- stats::glm(formula = totalHotspotPos_2minus ~ DNA_Accessibility_bins, data= data, family=poisson(link=log))
    
    # get data from focal hotspot mutation and calculate expected number of mutations (lambda) and get observed number of mutations
    focalHPdata = data.frame(DNA_Accessibility_bins = x[iHotspotMut,]$DNA_Accessibility_bins)
    lambda_par <- stats::predict(gml_PoissonModel, focalHPdata, type="response")
    n_muts <- x[iHotspotMut,]$totalHotspotPos_2minus
    
    x$pValue[iHotspotMut] <- poisson.test(n_muts, r=lambda_par, alternative = 'greater')$p.value
  }
  
  
  return(x)
  
  # define number of cores
}, cl = 2)

# get all results in data frame
summResults_gmaPoissModel_noHairpin <- BiocGenerics::do.call(BiocGenerics::rbind, summResults_gmaPoissModel_noHairpin)

# Adjust pValues for multiple testing
summResults_gmaPoissModel_noHairpin <- summResults_gmaPoissModel_noHairpin %>% dplyr::mutate(pAdj = p.adjust(pValue, method = "BH"))

# Plot Pvalue
summResults_gmaPoissModel_noHairpin <- summResults_gmaPoissModel_noHairpin %>%
  dplyr::mutate(labelGene = ifelse(pAdj <= 0.05, geneName, ""))

# add hotspot id
summResults_gmaPoissModel_noHairpin <- summResults_gmaPoissModel_noHairpin %>%
  dplyr::mutate(hs_id = paste0(chr, ":", pos, " ", geneName))
  



########## driver hotspot mutations in hairpin loops ########################################################

# transform strings of ATCG to 0011 (G, C == 1; A, T == 0)
results_hairpinLoopsStability_Hairpin <- results_hairpinLoopsStability_Hairpin %>%
  dplyr::mutate(loopSeq_binary = gsub('A|T','0', loopSeq)) %>%
  dplyr::mutate(loopSeq_binary = gsub('G|C','1', loopSeq_binary))


# Is the hotspot happening in 1001 or 101 context in the loop?
results_hairpinLoopsStability_Hairpin <- results_hairpinLoopsStability_Hairpin %>%
  dplyr::mutate(isSeq1001or101_up = ifelse(base::substr(loopSeq_binary, start = loopPos-2, stop = loopPos) == '101' |
                                             base::substr(loopSeq_binary, start = loopPos-3, stop = loopPos) == '1001', 1, 0),
                isSeq1001or101_down = ifelse(base::substr(loopSeq_binary, start = loopPos, stop = loopPos+2) == '101' |
                                               base::substr(loopSeq_binary, start = loopPos, stop = loopPos+3) == '1001', 1, 0)) %>%
  dplyr::mutate(has_1001or101 = ifelse(isSeq1001or101_up+isSeq1001or101_down == 0, "No", "Yes"))


# scale number of hotspots to zero, for proper application of the poisson distirbution 
results_hairpinLoopsStability_Hairpin <- results_hairpinLoopsStability_Hairpin %>%
  dplyr::mutate(totalHotspotPos_2minus = totalHotspotPos - 2)


# Run in parallel modelling for poisson process (~1m30s)
summResults_gmaPoissModel_hairpin <- pbapply::pblapply(base::split(results_hairpinLoopsStability_Hairpin, results_hairpinLoopsStability_Hairpin$triNucleotideContext), function(x){

  # initialize a column to fill in p-values
  x$pValue <- 1

  # fit a model excluding the focal hotspot mutation and hotspots with more mutations than the focal hotspot mutation
  for(iHotspotMut in c(1:nrow(x))) {
    
    # get number of mutations of focal hotspot mut and termite program if nMutsFocal == 0
    nMutsFocal <- x[iHotspotMut,]$totalHotspotPos_2minus
    if (nMutsFocal == 0) {next}
    
    # exclude all hotspot muts with higuer frequency than nMutsFocal
    data <- x[-iHotspotMut,]
    data <- data %>% dplyr::filter(totalHotspotPos_2minus <= nMutsFocal)
    
    # fit data to poisson process
    gml_PoissonModel <- stats::glm(formula = totalHotspotPos_2minus ~ hairpinStab + DNA_Accessibility_bins + has_1001or101, data= data, family=poisson(link=log))
    
    # get data from focal hotspot mutation and calculate expected number of mutations (lambda) and get observed number of mutations
    focalHPdata = data.frame(hairpinStab = x[iHotspotMut,]$hairpinStab, DNA_Accessibility_bins=x[iHotspotMut,]$DNA_Accessibility_bins, has_1001or101= x[iHotspotMut,]$has_1001or101)
    lambda_par <- stats::predict(gml_PoissonModel, focalHPdata, type="response")
    n_muts <- x[iHotspotMut,]$totalHotspotPos_2minus

    x$pValue[iHotspotMut] <- poisson.test(n_muts, r=lambda_par, alternative = 'greater')$p.value
  }
  
  
  return(x)
  
  # define number of cores
}, cl = 2)

# get all results in data frame
summResults_gmaPoissModel_hairpin <- BiocGenerics::do.call(BiocGenerics::rbind, summResults_gmaPoissModel_hairpin)

# Adjust pValues for multiple testing
summResults_gmaPoissModel_hairpin <- summResults_gmaPoissModel_hairpin %>% dplyr::mutate(pAdj = p.adjust(pValue, method = "BH"))

# add labels for significant hotspot mutations
summResults_gmaPoissModel_hairpin <- summResults_gmaPoissModel_hairpin %>%
  dplyr::mutate(labelGene = ifelse(pAdj <= 0.05, geneName, ""))

# check drivers
summResults_gmaPoissModel_hairpin  %>% dplyr::select(totalHotspotPos, pAdj, geneName) %>% dplyr::arrange(pAdj)

# add hotspot id
summResults_gmaPoissModel_hairpin <- summResults_gmaPoissModel_hairpin %>%
  dplyr::mutate(hs_id = paste0(chr, ":", pos, " ", geneName))


# combine results of all driver APOBEC-derived hotspot mutations (hairpin and no hairpin loops)
summResults_gmaPoissModel_All <- rbind(dplyr::select(summResults_gmaPoissModel_hairpin, chr, pos, totalHotspotPos, formLoop, pAdj, labelGene, consequenceAnnotation, overlappingGenes),
                                       dplyr::select(summResults_gmaPoissModel_noHairpin, chr, pos, totalHotspotPos, formLoop, pAdj, labelGene, consequenceAnnotation, overlappingGenes)) %>%
  dplyr::arrange(pAdj)



# export results
write.csv(summResults_gmaPoissModel_All, file = "results/summResults_gmaPoissModel_All.csv", row.names = FALSE)







