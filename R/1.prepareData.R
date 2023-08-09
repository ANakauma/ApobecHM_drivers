# Author:  J. Alberto Nakauma-González
# Date:   03-08-2023
# e-mail: j.nakaumagonzalez@erasmusmc.nl 
# Function: Prepare the data set for identifying drivers
#           a) Identifies the trinucleotide context of single nucleotide variants
#           b) Calculates APOBEC enrichment score in tri- and tetra-nucleotide context
#           c) Identifies APOBEC positive samples
#           d) Kataegis estimation


# Clean everything before starting a new project and set working directory (default R script path)---------------------------
rm(list=ls())
dev.off()
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load libraries ----------------------------------------------------------
pacman::p_load('plyr', 'dplyr', 'tidyr', 'BiocParallel', 'GenomicRanges', 'IRanges')

# load reference human genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
#ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
require(ref_genome, character.only = TRUE)


# Auxiliary functions ----------------------------------------------------
#function to change to 6 base substitutions types
convertMutationalContext <- function(muts){
  
  muts = base::gsub('[G>T]', '[C>A]', muts, fixed = T)
  muts = base::gsub('[G>C]', '[C>G]', muts, fixed = T)
  muts = base::gsub('[G>A]', '[C>T]', muts,fixed = T)
  muts = base::gsub('[A>T]', '[T>A]', muts, fixed = T)
  muts = base::gsub('[A>G]', '[T>C]', muts, fixed = T)
  muts = base::gsub('[A>C]', '[T>G]', muts, fixed = T)
  
  return(muts)
}

#return complement base pair when mutated base pair has changed
convertToComplementary <- function(bp){
  new_bp = ifelse(bp == "A", "T",
                  ifelse(bp == "T", "A",
                         ifelse(bp == "G", "C", "G")))
  return(new_bp)
}


# Import mutations --------------------------------------------------------
## Import data set, only SNVs are needed with at least sampleId, chr, pos, ref and alt columns 
SNVs_data <- read.table(file = "data/SNV_data.txt", header = TRUE)

#Check that alt and ref are 1L length
SNVs_data <- SNVs_data %>% dplyr::filter(nchar(ref) == 1 & nchar(alt) == 1)

#Delete chr M and Y in both data sets
SNVs_data <- SNVs_data %>% subset(chr!="chrY") %>% subset(chr!="chrM")
SNVs_data$chr <- factor(SNVs_data$chr)
SNVs_data$chr <- factor(SNVs_data$chr, levels = unique(SNVs_data$chr))



#----------------------------------------------------------------------------------------------------------------
#################################### Find trinucleotide context for SNVs ###############################
#----------------------------------------------------------------------------------------------------------------

# Find the upstream and downstream nucleotides for each SNV
SNVs_data$upstreamNuc <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, SNVs_data$chr, SNVs_data$pos - 1, SNVs_data$pos - 1))
SNVs_data$downstreamNuc <- as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, SNVs_data$chr, SNVs_data$pos + 1, SNVs_data$pos + 1))

#Add mutContext
SNVs_data$mutContext = sprintf('%s[%s>%s]%s', SNVs_data$upstreamNuc, SNVs_data$ref, SNVs_data$alt, SNVs_data$downstreamNuc)
  
# Convert to the six types of mutational contexts
SNVs_data$mutContext_6 <- SNVs_data$mutContext
SNVs_data$mutContext_6 <- convertMutationalContext(SNVs_data$mutContext_6)
  
#Change flanking bases to their inverse complementary for cases where mutations have been converted
SNVs_data <- SNVs_data %>% dplyr::mutate(downstreamNuc_6 = ifelse(mutContext == mutContext_6,
                                                                  as.character(downstreamNuc),
                                                                  convertToComplementary(upstreamNuc)))
SNVs_data <- SNVs_data %>% dplyr::mutate(upstreamNuc_6 = ifelse(mutContext == mutContext_6,
                                                                as.character(upstreamNuc),
                                                                convertToComplementary(downstreamNuc)))

# update the final mutational context  
SNVs_data$mutContext_6Updated = sprintf('%s[%s>%s]%s', SNVs_data$upstreamNuc_6, SNVs_data$ref, SNVs_data$alt, SNVs_data$downstreamNuc_6)
SNVs_data$mutContext_6Updated <- convertMutationalContext(SNVs_data$mutContext_6Updated)
SNVs_data <- SNVs_data %>% dplyr::mutate(upstreamNuc = upstreamNuc_6,
                                           downstreamNuc = downstreamNuc_6,
                                           mutContext = mutContext_6Updated)
SNVs_data$upstreamNuc_6 <- NULL
SNVs_data$downstreamNuc_6 <- NULL
SNVs_data$mutContext_6Updated <- NULL
SNVs_data$mutContext_6 <- NULL
  

# Find APOBEC mutations (T[C->T]A, T[C->T]T, T[C->G]A, T[C->G]T)
SNVs_data$APOBEC <- ifelse(SNVs_data$mutContext == 'T[C>T]A' |
                             SNVs_data$mutContext == 'T[C>T]T' |
                             SNVs_data$mutContext == 'T[C>G]A' |
                             SNVs_data$mutContext == 'T[C>G]T', 1, 0)
  
# Specific Mutational Changes (6 combinations)
SNVs_data <-  SNVs_data %>% dplyr::mutate(mutChange = gsub('].*', '', gsub('.*\\[', '', mutContext)) )
  
# save data 
save(SNVs_data, file = "/APOBECproject/RData/SNVs_data.RData")



#----------------------------------------------------------------------------------------------------------------
#################################### Calculate APOBEC Enrichment ###############################
#----------------------------------------------------------------------------------------------------------------


# Number of C+G and TCW (WGA) bases in the whole genome
total_Context_CorG <- 1100439365
total_Context_TCW <- 223735123
# For hg38 use:
#total_Context_CorG <- 1128757468
#total_Context_TCW <- 231224774

# Count total Number of C>T and C>G mutations and APOBEC (TCW>TTW and TCW>TGW mutations) muts (Exclude sex chromosomes)
APOBECClassSample <- SNVs_data %>% dplyr::mutate(C_to_T = ifelse(mutChange == 'C>T', 1, 0), C_to_G = ifelse(mutChange == 'C>G', 1, 0), ) %>%
  dplyr::filter(!(chr %in% c('chrX', 'chrY'))) %>% dplyr::group_by(sampleId) %>%
  dplyr::add_count(name = 'totalSNVs') %>%
  dplyr::mutate(APOBEC_Freq = sum(APOBEC), C_T_Freq = sum(C_to_T), C_G_Freq = sum(C_to_G)) %>% dplyr::ungroup() %>%
  dplyr::distinct(sampleId, totalSNVs, APOBEC_Freq, C_T_Freq, C_G_Freq) %>%
  dplyr::mutate(total_CtoT_CtoG = C_G_Freq + C_T_Freq)

# Calculate fold enrichment for APOBEC mutations
APOBECClassSample <- APOBECClassSample %>%
  dplyr::mutate(foldEnrichment = (APOBEC_Freq * total_Context_CorG) / (total_CtoT_CtoG * total_Context_TCW))

# Initialize pValues for test
APOBECClassSample$E_pValue <- 1

# Calculate enrichment pValues
for (iPatient in c(1:nrow(APOBECClassSample))) {
  challenge.df = matrix(c(APOBECClassSample$APOBEC_Freq[iPatient], total_Context_TCW,
                          APOBECClassSample$total_CtoT_CtoG[iPatient]-APOBECClassSample$APOBEC_Freq[iPatient], total_Context_CorG - total_Context_TCW), nrow = 2)
  
  APOBECClassSample$E_pValue[iPatient] <- fisher.test(challenge.df, alternative = "greater")$p.value
}

#Adjust p-values with Benjamini-Hochberg method
APOBECClassSample$E_pAdj <- p.adjust(APOBECClassSample$E_pValue, method = "BH")

# Tumor is APOBEC enriched if p_Adj < 0.05
APOBECClassSample <- APOBECClassSample %>% dplyr::mutate(APOBEC_enrich = ifelse(E_pAdj < 0.05, "Yes", "No"))
APOBECClassSample$APOBEC_enrich <- factor(APOBECClassSample$APOBEC_enrich, levels = c("No", "Yes"))

# classify tumors into high, medium, low or no APOBEC mutagenesis 
APOBECClassSample <- APOBECClassSample %>%
  dplyr::mutate(APOBEC_mutagenesis = ifelse(APOBEC_enrich == "Yes",
                                            ifelse(foldEnrichment > 3, "High", ifelse(foldEnrichment > 2, 'Medium', "Low")), "No"))
APOBECClassSample$APOBEC_mutagenesis <- factor(APOBECClassSample$APOBEC_mutagenesis, levels = c("No", "Low", "Medium", "High"))





#----------------------------------------------------------------------------------------------------------------
########################### estimate tetra nucleotide context for APOBEC mutations ###############################
#----------------------------------------------------------------------------------------------------------------

# Count total SNVs and keep mutations with NTCA context (Exclude sex chr)
tetraContext_APOBEC <- SNVs_data %>% dplyr::filter(!(chr %in% c("chrX", "chrY"))) %>%
  dplyr::left_join(dplyr::select(APOBECClassSample, sampleId, totalSNVs, total_CtoT_CtoG), by = "sampleId") %>%
  dplyr::filter(APOBEC == 1) %>%
  dplyr::filter(downstreamNuc == "A")

# Get info for mut context tetra nucleotide
tetraContext_APOBEC <- tetraContext_APOBEC %>%
  dplyr::mutate(mutContextTetra_nuc = as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, chr, pos - 2, pos + 2))) %>%
  dplyr::mutate(upstreamNuc2_tmp = as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, chr, pos - 2, pos - 2))) %>%
  dplyr::mutate(downstreamNuc2_tmp = as.character(Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, chr, pos + 2, pos + 2)))


# convert to 6 types of mutations (C and T, change A>T and G>C)
tetraContext_APOBEC <- tetraContext_APOBEC %>%
  dplyr::mutate(upstreamNuc2_tmp = ifelse(ref == "G" | ref == "A", convertToComplementary(upstreamNuc2_tmp), upstreamNuc2_tmp),
                downstreamNuc2_tmp = ifelse(ref == "G" | ref == "A", convertToComplementary(downstreamNuc2_tmp), downstreamNuc2_tmp))

# get mutational tetra context for all APOBEC mutations
tetraContext_APOBEC <- tetraContext_APOBEC %>%
  dplyr::mutate(upstreamNuc2 = ifelse(ref == "G" | ref == "A", downstreamNuc2_tmp, upstreamNuc2_tmp)) %>%
  dplyr::mutate(mutContextTetra = paste0(upstreamNuc2, mutContext))

# Define YTCA (A3A; Y = T or C) or RTCA (A3B; R = G or A)
tetraContext_APOBEC <- tetraContext_APOBEC %>%
  dplyr::mutate(APOBECmut_tetraContext = ifelse(upstreamNuc2 == "T" | upstreamNuc2 == "C", "YTCA", "RTCA"))

# Count total mutations in tetra context per sample
tetraContext_APOBEC <- tetraContext_APOBEC %>%
  dplyr::add_count(sampleId, name = "totalNTCAmuts") %>%
  dplyr::mutate(isYTCA = ifelse(APOBECmut_tetraContext == "YTCA", 1, 0),
                isRTCA = ifelse(APOBECmut_tetraContext == "RTCA", 1, 0)) %>%
  dplyr::group_by(sampleId) %>%
  dplyr::mutate(totalYTCA = sum(isYTCA), totalRTCA = sum(isRTCA)) %>% dplyr::ungroup() %>%
  dplyr::mutate(YTCW_RTCW_ratio = totalYTCA/totalRTCA) %>%
  dplyr::distinct(sampleId, .keep_all = TRUE) %>%
  dplyr::select(sampleId, totalSNVs, totalYTCA, totalRTCA, total_CtoT_CtoG)


# Enrichment for YTCA and for RTCA
########### Add APOBEC mut_load Enrichment calculations ###########
# Number of context bases in the whole genome  excluding sex and uncharacterized chromosomes
total_Context_CorG <- 1100439365
total_Context_YTCA <- 63373254
total_Context_RTCA <- 41712804
# For hg38 use:
#total_Context_CorG <- 1128757468
#total_Context_YTCA <- 65769909
#total_Context_RTCA <- 42574131

# Calculate fold enrichment for YTCA mutations
tetraContext_APOBEC <- tetraContext_APOBEC %>%
  dplyr::mutate(foldEnrichment_YTCA = (totalYTCA * total_Context_CorG) / (total_CtoT_CtoG * total_Context_YTCA),
                foldEnrichment_RTCA = (totalRTCA * total_Context_CorG) / (total_CtoT_CtoG * total_Context_RTCA))
tetraContext_APOBEC <- tetraContext_APOBEC %>%
  dplyr::select(sampleId, totalYTCA, totalRTCA, foldEnrichment_YTCA, foldEnrichment_RTCA)

# Get the complete APOBEC class data
APOBECClassSample <-  APOBECClassSample %>% dplyr::full_join(tetraContext_APOBEC, by = 'sampleId')

# Save APOBEC class data 
save(APOBECClassSample, file = "/APOBECproject/RData/APOBECClassSample.RData")



#------------------------------------ Estimate Kataegis -------------------------------------------------
# install R package
BiocManager::install("katdetectr")
library(katdetectr)

# group SNVs per patient, including the read counts for the reference (AD_Ref), the alteration (AD_alt) and total depth (DP) 
snvRangesList <- base::split(SNVs_data %>% dplyr::select(chr, pos, ref, alt, AD_Ref, AD_alt, DP, sample), SNVs_data$sampleId)

# Estimate kataegis per sample (cl = number of cores)
cohortKataegis <- pbapply::pblapply(snvRangesList, function(x){
  vR <- VariantAnnotation::VRanges(seqnames =  x$chr, IRanges(x$pos, x$pos),
                                   ref = x$ref, alt = x$alt,
                                   refDepth = x$AD_Ref, altDepth = x$AD_alt, totalDepth = x$DP,
                                   sampleNames = x$sample)
  
  katDetect_result <- katdetectr::detectKataegis(vR)
  return(katDetect_result)
}, cl = 8)

# save result
save(cohortKataegis, file = "/APOBECproject/RData/cohortKataegis.RData")







