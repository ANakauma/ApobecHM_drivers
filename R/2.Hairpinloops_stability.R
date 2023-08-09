# Author: Minouk Noordsij + J. Alberto Nakauma-González
# Date: 03-08-2023
# e-mail: j.nakaumagonzalez@erasmusmc.nl 
# This script predicts the hairpin structures around mutation sites (or any genomic position):
#     a) Estimates hairpin-loops and the thermodynamic stability
#     b) Identifies Twin mutations and Didymi

# Clean everything and load packages and reference genome ---------------------
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Load packages
pacman::p_load('plyr', 'dplyr', 'tidyr', 'BiocParallel', 'readODS', 'stringr')

# load reference genomes (both can be uploaded or only the required one)
require("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
require("BSgenome.Hsapiens.UCSC.hg38", character.only = TRUE)




# ------------------------------------------------------------------------------------------------------
# -------------------------------------  Auxiliary functions  ------------------------------------------
# ------------------------------------------------------------------------------------------------------
# Function to change to 6 base substitutions types
convertMutationalContext <- function(muts){
  
  muts = base::gsub('[G>T]', '[C>A]', muts, fixed = T)
  muts = base::gsub('[G>C]', '[C>G]', muts, fixed = T)
  muts = base::gsub('[G>A]', '[C>T]', muts, fixed = T)
  muts = base::gsub('[A>T]', '[T>A]', muts, fixed = T)
  muts = base::gsub('[A>G]', '[T>C]', muts, fixed = T)
  muts = base::gsub('[A>C]', '[T>G]', muts, fixed = T)
  
  return(muts)
}

# Function that will return complement base pair for a single base
convertToComplementary <- function(bp){
  new_bp = ifelse(bp == "A", "T",
                  ifelse(bp == "T", "A",
                         ifelse(bp == "G", "C", "G")))
  return(new_bp)
}

# Function that will return complement sequence for multiple bases as input
# Input can be either a single string containing all bases or 
# a vector containing individual bases
convertToComplementarySeq <- function(bps){
  new_bps = ""
  splitstring <- TRUE
  if (nchar(bps)[1]!=1){
    bps <- strsplit(bps, "")[[1]]
    splitstring <- FALSE}
  for (bp in bps){
    new_bp <- ifelse(bp == "A", "T",
                     ifelse(bp == "T", "A",
                            ifelse(bp == "G", "C", "G")))
    
    new_bps <- paste(new_bps,new_bp,sep="")}
  if (splitstring==TRUE){
    new_bps <- strsplit(new_bps, "")[[1]]}
  return(new_bps)
}

#Function that will return reverse complement sequence
# Input can be either a single string containing all bases or 
# a vector containing individual bases
convertToReverseComplementarySeq <- function(bps){ 
  new_bps = ""
  splitstring <- TRUE
  if (nchar(bps)[1]!=1){
    bps<- strsplit(bps, "")[[1]]
    splitstring <- FALSE}
  for (bp in bps){
    new_bp <- ifelse(bp == "A", "T",
                     ifelse(bp == "T", "A",
                            ifelse(bp == "G", "C", "G")))
    
    new_bps <- paste(new_bp,new_bps,sep="")}
  if (splitstring==TRUE){
    new_bps <- strsplit(new_bps, "")[[1]]}
  return(new_bps)
}

#Function that will return reverse sequence
# Input can be either a single string containing all bases or 
# a vector containing individual bases
convertToReverseSeq <- function(bps){ 
  splitstring <- TRUE
  if (nchar(bps)[1]!=1){
    bps <- strsplit(bps, "")[[1]]
    splitstring <- FALSE
  }
  new_bps <- bps[length(bps):1]
  if (splitstring==FALSE){
    new_bps <- paste(new_bps,sep="",collapse="")
  }
  return(new_bps)
}


# ------------------------------------------------------------------------------------------------------
# ------------------------------------  Parameter files for functions ----------------------------------
# ------------------------------------------------------------------------------------------------------
# Used for nearest neighbour stability (Themodynamics)
terminal_mismatches = readODS::read_ods("data/TerminalMismatches.ods")
triloops = readODS::read_ods("data/Triloops.ods")
tetraloops = readODS::read_ods("data/Tetraloops.ods")
loopsize = readODS::read_ods("data/Loopsize.ods")
NN_parameters = readODS::read_ods("data/NNparameters.ods")

# Used for substrate optimality as estimated by Buisson et al., Science 2019 
substrate_optimality = readODS::read_ods("data/BuissonSubstrateOptimality.ods")
buisson_triloops = readODS::read_ods("data/BuissonTriLoops.ods")
buisson_triloops$context=stringr::str_remove_all(buisson_triloops$context, "[()]")
buisson_tetraloops = readODS::read_ods("data/BuissonTetraLoops.ods")
buisson_tetraloops$context=stringr::str_remove_all(buisson_tetraloops$context, "[()]")




# ------------------------------------------------------------------------------------------------------
# ------------------------------------  Function for hairpin stability  --------------------------------
# ------------------------------------------------------------------------------------------------------

# Function substrateOptimality ------------------------------------------------
# Function that returns the predicted substrate optimality of a hairpin with
# stem_seq, the DNA sequence of matching part of the the hairpin stem (leaving 
# out any mismatches),
# loop_pos, the position (integer) of the mutation within the loop,
# loop_seq, the loop sequence including the first matching basepair 5' to 3' and
# mismatch_count, the number (integer) of mismatches
substrateOptimality <- function(stem_seq, loop_pos, loop_seq, mismatch_count){
  loop_size <- nchar(loop_seq)-2
  stability_score <- 0
  # sum of 1*#AT basepair + 3*#CG basepair in the stem
  for (nt in stem_seq){
    stability_score <- stability_score + ifelse(nt=="A" || nt=="T", 1, 3)
  }
  # mismatch penalty of -2 
  stability_score <- stability_score - (2*mismatch_count)
  
  if(loop_size==3&&loop_pos==3){
    if(stability_score<5){
      substrate_optimality_value<-(subset(buisson_triloops, context==loop_seq))[2]
    }
    else if(stability_score<7){
      substrate_optimality_value<-(subset(buisson_triloops, context==loop_seq))[3]
    }
    else if(stability_score<9){
      substrate_optimality_value<-(subset(buisson_triloops, context==loop_seq))[4]
    }
    else if(stability_score<11){
      substrate_optimality_value<-(subset(buisson_triloops, context==loop_seq))[5]
    }
    else if(stability_score<13){
      substrate_optimality_value<-(subset(buisson_triloops, context==loop_seq))[6]
    }
    else if(stability_score<15){
      substrate_optimality_value<-(subset(buisson_triloops, context==loop_seq))[7]
    }
    else if(stability_score<17){
      substrate_optimality_value<-(subset(buisson_triloops, context==loop_seq))[8]
    }
    else if(stability_score<19){
      substrate_optimality_value<-(subset(buisson_triloops, context==loop_seq))[9]
    }
    else{
      substrate_optimality_value<-(subset(buisson_triloops, context==loop_seq))[10]
    }
  }
  
  else if(loop_size==4&&(loop_pos==3||loop_pos==4)){
    if(stability_score<5){
      substrate_optimality_value<-(buisson_tetraloops %>% subset(looppos==loop_pos) %>% subset(context==loop_seq))[4]
    }
    else if(stability_score<7){
      substrate_optimality_value<-(buisson_tetraloops %>% subset(looppos==loop_pos) %>% subset(context==loop_seq))[5]
    }
    else if(stability_score<9){
      substrate_optimality_value<-(buisson_tetraloops %>% subset(looppos==loop_pos) %>% subset(context==loop_seq))[6]
    }
    else if(stability_score<11){
      substrate_optimality_value<-(buisson_tetraloops %>% subset(looppos==loop_pos) %>% subset(context==loop_seq))[7]
    }
    else if(stability_score<13){
      substrate_optimality_value<-(buisson_tetraloops %>% subset(looppos==loop_pos) %>% subset(context==loop_seq))[8]
    }
    else{
      substrate_optimality_value<-(buisson_tetraloops %>% subset(looppos==loop_pos) %>% subset(context==loop_seq))[9]
    }
  }
  
  else if(loop_size>2&&loop_size<12){
    if(stability_score<7){
      substrate_optimality_value<-(substrate_optimality %>% subset(looppos==loop_pos) %>% subset(looplen==loop_size))[3]
    }
    else if(stability_score<9){
      substrate_optimality_value<-(substrate_optimality %>% subset(looppos==loop_pos) %>% subset(looplen==loop_size))[4]
    }
    else if(stability_score<11){
      substrate_optimality_value<-(substrate_optimality %>% subset(looppos==loop_pos) %>% subset(looplen==loop_size))[5]
    }
    else if(stability_score<13){
      substrate_optimality_value<-(substrate_optimality %>% subset(looppos==loop_pos) %>% subset(looplen==loop_size))[6]
    }
    else if(stability_score<15){
      substrate_optimality_value<-(substrate_optimality %>% subset(looppos==loop_pos) %>% subset(looplen==loop_size))[7]
    }
    else if(stability_score<17){
      substrate_optimality_value<-(substrate_optimality %>% subset(looppos==loop_pos) %>% subset(looplen==loop_size))[8]
    }
    else if(stability_score<19){
      substrate_optimality_value<-(substrate_optimality %>% subset(looppos==loop_pos) %>% subset(looplen==loop_size))[9]
    }
    else{
      substrate_optimality_value<-(substrate_optimality %>% subset(looppos==loop_pos) %>% subset(looplen==loop_size))[10]
    }
  }
  
  else{
    print("Incorrect loop size")
    return(-1)
  }
  
  return(substrate_optimality_value[[1]])
}



# Function nearestNeighbourStability --------------------------------------------
# This function computes the stability of a hairpin loop using nearest 
# neighbour parameters with
# stem_seq, the downstream stem sequence from 5' to 3',
# comp_stem_seq, the upstream stem sequence from 3' to 5',
# loop_seq, the loop sequence including the closing matching basepair 5' to 3',
# terminal_mismatch, the 5' mismatch + closing basepair + 3' mismatch at the stembase,
# bulge_loops_stem, position of bulge loops (if present) in the downstream stem sequence and
# bulge_loops_comp_stem, position of bulge loops (if present) in the upstream stem sequence
nearestNeighbourStability <- function(stem_seq, comp_stem_seq, loop_seq, terminal_mismatch, bulge_loops_stem, bulge_loops_comp_stem){
  loop_size <- nchar(loop_seq)-2
  if(loop_size<3 || loop_size>20){
    return(0)
  }
  
  stem_length_up <- length(comp_stem_seq)-length(bulge_loops_comp_stem)
  
  stem_length_down <- length(stem_seq)-length(bulge_loops_stem)
  
  if(stem_length_up!=stem_length_down){
    return("Error with stem length")
  }
  
  stability_score <- terminal_mismatches[terminal_mismatches$Seq==terminal_mismatch,2] + 
    loopsize[loopsize$LoopSize==loop_size,2] +
    4.0*length(bulge_loops_comp_stem) + 4.0*length(bulge_loops_stem)
  
  if (loop_size==3){
    if(nrow(triloops[triloops$Seq==loop_seq,])!=0){
      stability_score <- stability_score + triloops[triloops$Seq==loop_seq,2]
    }
  }
  else if (loop_size==4){
    if(nrow(tetraloops[tetraloops$Seq==loop_seq,])!=0){
      stability_score <- stability_score + tetraloops[tetraloops$Seq==loop_seq,2]
    }
  }
  
  prev_nt <- "X"
  skip <- 0
  skip_comp <- 0
  for (i in 1:stem_length_down){
    
    if (i %in% bulge_loops_stem){ 
      skip <- skip + 1}
    nt <- stem_seq[i+skip]
    if (i %in% bulge_loops_comp_stem){ 
      skip_comp <- skip_comp + 1}
    comp_nt <- comp_stem_seq[i+skip_comp]
    
    if (prev_nt!="X"){
      prop_seq <- paste(prev_nt,nt,sep="")
      
      if (convertToComplementary(prev_nt)==prev_comp_nt){
        Y <- ifelse(comp_nt=="A", 2, 
                    ifelse(comp_nt=="C", 3,
                           ifelse(comp_nt=="G", 4, 5)))
        stability_score <- stability_score + NN_parameters[NN_parameters$Seq==prop_seq,Y]}
      else if (convertToComplementary(nt)==comp_nt){
        prop_seq <- paste(comp_nt,prev_comp_nt,sep="")
        Y <- ifelse(prev_nt=="A", 2, 
                    ifelse(prev_nt=="C", 3,
                           ifelse(prev_nt=="G", 4, 5)))
        stability_score <- stability_score + NN_parameters[NN_parameters$Seq==prop_seq,Y]
      }
      else{
        return("Impossible to have adjacent mismatches")
      }
    }
    
    prev_nt <- nt
    prev_comp_nt <- comp_nt
  }
  
  return(stability_score)
}


# Function findHairpin ----------------------------------------------------
# This function looks for a hairpin loop around a cytosine mutation site and
# it returns the hairpin with the best stability
# downstream_seq: the sequence downstream from the mutation 5' to 3'
# upstream_seq: the sequence upstream from the mutation 3' to 5'
# max_loopsize: maximum number of nucleotides allowed in the loop
# min_stemsize: minimum number of matching basepairs at the beginning of the stem
# mismatches_allowed: if TRUE single bp mismatches and single nt bulge loops are allowed
findHairpin <- function(downstream_seq, upstream_seq, max_loopsize, min_stemsize, mismatches_allowed){
  hairpin <- c(0,0,0,0,0,1,1)
  new_hairpin <- c(NA,NA,NA,NA,NA,NA,NA)
  '%notin%' <- Negate('%in%')
  
  for (upstream_shift in 1:max_loopsize){ #Shift downstream 
    seq_up <- upstream_seq[upstream_shift:(upstream_shift+(min_stemsize-1))]
    com_seq_up <- convertToComplementarySeq(seq_up)
    for (downstream_shift in 1:(max_loopsize-upstream_shift)){ #Shift upstream
      loop_size <- upstream_shift+downstream_shift-1
      if (downstream_shift == 0 || loop_size<3){  
        next #minimum loop size of 3
      }
      if (((upstream_shift == 1 && downstream_seq[downstream_shift-1]=="G")
           || (downstream_shift == 1 && upstream_seq[upstream_shift-1]=="G"))
          && loop_size>4){
        next #makes sure the first 2 bases of the loop are not a matching GC pair if the loop is bigger than 4
      }
      seq_down = downstream_seq[downstream_shift:(downstream_shift+(min_stemsize-1))]
      if(all(com_seq_up == seq_down)){ #Check if sequences match
        all_hairpins_found <- FALSE
        seq_loop <- paste(append(append(upstream_seq[(upstream_shift):1], "C"), 
                                 downstream_seq[1:(downstream_shift)]),collapse="")
        mismatch_iter <- 1
        iter <- 0
        bulge_found <- c()
        while(all_hairpins_found==FALSE){
          iter <- iter + 1
          seq_stem <- seq_down
          seq_stem_no_mm <- seq_down
          seq_comp_stem <- seq_up
          
          bulge_down <- c()
          bulge_up <- c()
          
          match <- TRUE
          mismatch <- FALSE
          match_count <- 1
          mismatch_count <- 0
          down_count <- 0 
          up_count <- 0
          
          while (match==TRUE){
            down_count <- down_count + 1
            up_count <- up_count + 1
            
            # finish is tail is very long (Max number of nucleotides has been reached)
            if(is.na(downstream_seq[downstream_shift+(min_stemsize-1)+down_count+1]) |
               is.na(upstream_seq[upstream_shift+(min_stemsize-1)+up_count+1])){break}
            
            # check if next 2 nucleotides match
            if(downstream_seq[downstream_shift+(min_stemsize-1)+down_count] == 
               convertToComplementary(upstream_seq[upstream_shift+(min_stemsize-1)+up_count])){
              match_count <- match_count + 1
              if (mismatch==TRUE){
                seq_stem <- append(seq_stem, downstream_seq[downstream_shift+(min_stemsize-1)+down_count-1])
                seq_comp_stem <- append(seq_comp_stem, upstream_seq[upstream_shift+(min_stemsize-1)+up_count-1])
                match_count <- match_count + 1
                mismatch <- FALSE
              }
              
              seq_stem <- append(seq_stem, downstream_seq[downstream_shift+(min_stemsize-1)+down_count])
              seq_stem_no_mm <- append(seq_stem_no_mm, downstream_seq[downstream_shift+(min_stemsize-1)+down_count])
              seq_comp_stem <- append(seq_comp_stem, upstream_seq[upstream_shift+(min_stemsize-1)+up_count])
            } 
            else{ # mismatch
              # calculate the stability of current sequence
              terminal_mismatch <- paste0(upstream_seq[upstream_shift+(min_stemsize-1)+up_count],
                                          upstream_seq[upstream_shift+(min_stemsize-1)+up_count-1],
                                          downstream_seq[downstream_shift+(min_stemsize-1)+down_count-1],
                                          downstream_seq[downstream_shift+(min_stemsize-1)+down_count],collapse="")
              stability <- nearestNeighbourStability(seq_stem,seq_comp_stem,seq_loop,terminal_mismatch,bulge_down,bulge_up)
              
              # when substrate optimality cannot be estimated, then it should be == 1 
              substrate_optimality <- substrateOptimality(seq_stem_no_mm,upstream_shift,seq_loop,mismatch_count)
              if(length(substrate_optimality) == 0) substrate_optimality = 1

              new_hairpin[1] <- min_stemsize + down_count - 1  #downstream stemlength
              new_hairpin[2] <- min_stemsize + up_count - 1    #upstream stemlength
              new_hairpin[3] <- loop_size                      #loopsize
              new_hairpin[4] <- upstream_shift                 #loop pos of mutation
              new_hairpin[5] <- stability                      #nearest neighbour stability
              new_hairpin[6] <- substrate_optimality           #substrate optimality
              hairpin <- rbind(hairpin, new_hairpin)
              
              
              if (mismatches_allowed==FALSE){
                match <- FALSE
                all_hairpins_found <- TRUE
              }
              else{
                mismatch_bulge_down <- FALSE
                mismatch_bulge_up <- FALSE
                if (upstream_seq[upstream_shift+(min_stemsize-1)+up_count] ==
                    convertToComplementary(downstream_seq[downstream_shift+(min_stemsize-1)+down_count+1]) 
                    && length(bulge_down)<1 && length(bulge_up)<1 
                    && (length(seq_stem)+1) %notin% bulge_found){
                  
                  mismatch_bulge_down <- TRUE #BULGE LOOP DOWNSTREAM
                }
                if (upstream_seq[upstream_shift+(min_stemsize-1)+up_count+1] ==
                    convertToComplementary(downstream_seq[downstream_shift+(min_stemsize-1)+down_count])
                    && length(bulge_down)<1 && length(bulge_up)<1 
                    && -(length(seq_comp_stem)+1) %notin% bulge_found){
                  mismatch_bulge_up <- TRUE #BULGE LOOP UPSTREAM
                }
                if (upstream_seq[upstream_shift+(min_stemsize-1)+up_count+1] ==
                    convertToComplementary(downstream_seq[downstream_shift+(min_stemsize-1)+down_count+1])){
                  mismatch <- TRUE  #MISMATCH
                }
                if (mismatch == FALSE && mismatch_bulge_down == FALSE && mismatch_bulge_up == FALSE){
                  match <- FALSE # no possibilities for continuing the sequence
                }
                else{
                  mismatch_count <- mismatch_count + 1
                  if (mismatch_bulge_down == TRUE){
                    seq_stem <- append(seq_stem, downstream_seq[downstream_shift+(min_stemsize-1)+down_count])
                    up_count <- up_count - 1
                    bulge_down <- append(bulge_down, length(seq_stem))
                    bulge_found <- append(bulge_found, length(seq_stem))
                    
                    if (mismatch_bulge_up == TRUE){
                      # bulge up and bulge down
                      mismatch_iter <- mismatch_iter + 1
                      all_hairpins_found <- FALSE
                    }
                    if (mismatch == TRUE){ 
                      # bulge down and mismatch
                      mismatch_iter <- mismatch_iter + 1
                      all_hairpins_found <- FALSE 
                      mismatch <- FALSE
                    }
                    
                  }
                  else if (mismatch_bulge_up == TRUE){
                    seq_comp_stem <- append(seq_comp_stem, upstream_seq[upstream_shift+(min_stemsize-1)+up_count])
                    down_count <- down_count - 1
                    bulge_up <- append(bulge_up, length(seq_comp_stem))
                    bulge_found <- append(bulge_found, -(length(seq_comp_stem)))
                    
                    
                    if (mismatch == TRUE){
                      # bulge up and mismatch
                      mismatch_iter <- mismatch_iter + 1
                      all_hairpins_found <- FALSE
                      mismatch <- FALSE
                    }
                    
                  }
                }
              } 
            } # end of else (mismatch)
          } # end of while(match==TRUE)
          if (mismatch_iter<=iter){
            all_hairpins_found <- TRUE
          }
        } # end of while(all_hairpins_found==FALSE)
        
      } 
    }
  }
  
  if (!is.null(dim(hairpin))){
    max_optimality <- max(hairpin[,6])
    hairpin <- hairpin[!duplicated(hairpin[,]),]
    if(min(hairpin[,5])>=0){
      hairpin = subset(hairpin,hairpin[,6] == max(hairpin[,6]))
      if ((dim(hairpin))[1]>1){
        hairpin <- subset(hairpin,hairpin[,5] == min(hairpin[,5]))
      }
    }
    else{
      hairpin <- subset(hairpin,hairpin[,5] == min(hairpin[,5]))
      if ((dim(hairpin))[1]>1){
        hairpin <- subset(hairpin,hairpin[,6] == max(hairpin[,6]))
      }
    }
    if ((dim(hairpin))[1]>1){
      hairpin <- subset(hairpin,hairpin[,3] == min(hairpin[,3]))
    }
    hairpin[7] <- max_optimality
  }
  return(hairpin) #Contains upstream start pos, downstream start pos, stemsize up, stemsize down, stability
}


# Function runFindHairpin in paralell -----------------------------------------------------
runFindHairpin <- function(data, hg_version = 'hg19', nCores = 6){

  # Reference genomes hg19 and hg38
  if(hg_version == 'hg19'){
    data$full_seq <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, data$chr, start=data$pos-50, end=data$pos+50, as.character=TRUE)
  }else if(hg_version == 'hg38'){
    data$full_seq = Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, data$chr, start=data$pos-50, end=data$pos+50, as.character=TRUE)
  }else{
    print('Please specify hg19 or hg38')
    return()
  }
  

  # split the data into list for parallel computing
  data_splitList <- split(data, 1:nrow(data))
  
  # Read all (main/summary) results
  hairpinStabilityValues <- pbapply::pblapply(data_splitList, function(iMutPos){

    full_seq <- strsplit(iMutPos$full_seq,"")[[1]]
    up_seq <- full_seq[52:101]
    down_seq <- convertToReverseSeq(full_seq[1:50])
    
    hairpin <- findHairpin(up_seq, down_seq, max_loopsize = 10, min_stemsize = 2, mismatches_allowed = TRUE)
    iMutPos$hairpinPosDown <-iMutPos$pos + hairpin[3]-hairpin[4] +1
    iMutPos$hairpinPosUp <- iMutPos$pos - hairpin[4]
    iMutPos$hairpinStemDown <- hairpin[1]
    iMutPos$hairpinStemUp <- hairpin[2]
    iMutPos$loopSize <- hairpin[3]
    iMutPos$loopPos <- hairpin[4]
    iMutPos$hairpinStab <- ifelse(hairpin[5]<0,hairpin[5],0)
    iMutPos$substrateOptimality <- hairpin[6]
    iMutPos$maxSubstrateOptimality <- hairpin[7]
  
  
  return(iMutPos)
    
    # define number of cores
  }, cl = nCores)
  
  # get all results in data frame
  hairpinStabilityValues <- BiocGenerics::do.call(BiocGenerics::rbind, hairpinStabilityValues)
  
  return(hairpinStabilityValues)
}




# ------------------------------------------------------------------------------------------------------
# ----------------------- Calculate hairpin loop stability for all mutations  --------------------------
# ------------------------------------------------------------------------------------------------------
# import all SNVs, ideally with trinucleotide context from 1.prepareData.R
load("/APOBECproject/RData/SNVs_data.RData")

# Only keep one variant per position (to not repeat calculation for hotspot mutations)
SNV_data_uniquePos <- SNVs_data %>% dplyr::distinct(chr, pos)


# Calculate hairpin loop stability for all hotspot mutations
# You can run it all at once, but we recommend to dot it in batches of 10,000 mutations in case of not having a powerful machine
iBatch <- c(1, seq(from = 10000, to = nrow(SNV_data_uniquePos), by = 10000), nrow(SNV_data_uniquePos))
iBatch <- unique(iBatch)

results_hairpinLoopsStability <- NULL
for (i in c(2:length(iBatch))) {
  print(i)
  if (is.null(results_hairpinLoopsStability)) {
    results_hairpinLoopsStability <- runFindHairpin_parallel(SNV_data_uniquePos[iBatch[1]:iBatch[2],], hg_version = 'hg19', nCores = 4)
  }else{
    results_hairpinLoopsStability <- rbind(results_hairpinLoopsStability, runFindHairpin_parallel(SNV_data_uniquePos[(iBatch[i-1]+1):iBatch[i],], hg_version = 'hg19', nCores = 4))
  }
}


# Get information per hotspot site (exclude failed cases = usually around TATATA repeated regions)
results_hairpinLoopsStability <- results_hairpinLoopsStability %>%
  dplyr::filter(full_seq != "Error in new_hairpin[5] <- stability : replacement has length zero\n" ) %>%
  dplyr::mutate(pos = as.numeric(pos),
                hairpinPosDown = as.numeric(hairpinPosDown),
                hairpinPosUp = as.numeric(hairpinPosUp),
                hairpinStemDown = as.numeric(hairpinStemDown),
                hairpinStemUp = as.numeric(hairpinStemUp),
                loopSize = as.numeric(loopSize),
                loopPos = as.numeric(loopPos),
                hairpinStab = as.numeric(hairpinStab),
                substrateOptimality = as.numeric(substrateOptimality),
                maxSubstrateOptimality = as.numeric(maxSubstrateOptimality)) %>%
  dplyr::full_join(SNVs_data, by = c('chr', 'pos'))


# Identify mutations that form hairpin loops and those that do not and TpC mutations
results_hairpinLoopsStability <- results_hairpinLoopsStability %>%
  dplyr::mutate(formLoop = ifelse(loopSize == 0, 0, 1),
                isTpC = base::grepl('T\\[C', mutContext))

# Add loop seq
results_hairpinLoopsStability <- results_hairpinLoopsStability %>%
  dplyr::mutate(loopSeq = ifelse(formLoop == 0, NA,
                                 Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, results_hairpinLoopsStability$chr,
                                                    start=results_hairpinLoopsStability$hairpinPosUp + 1,
                                                    end=results_hairpinLoopsStability$hairpinPosDown - 1, as.character=TRUE)))

results_hairpinLoopsStability$full_seq <- NULL


# save results
save(results_hairpinLoopsStability, file = "/APOBECproject/RData/results_hairpinLoopsStability.RData")





#-----------------------------------------------------------------------------------------------------------
# ----------------------------------- Identify Twin mutations and didymi 
#-----------------------------------------------------------------------------------------------------------

results_TwinMutations <- dplyr::filter(results_hairpinLoopsStability, formLoop == 1) %>%
  dplyr::add_count(chr, pos, name = "totalHotspotPos") %>%
  dplyr::distinct(chr, pos, .keep_all = TRUE) %>%
  dplyr::mutate(loopGenPos = paste0(chr, ":", hairpinPosUp, '-', hairpinPosDown)) %>%
  dplyr::add_count(loopGenPos, name = 'nHotspotPosLoop') %>%
  dplyr::filter(nHotspotPosLoop == 2) %>%
  dplyr::group_by(loopGenPos) %>%
  dplyr::mutate(nHotspotMutsInLoop = sum(totalHotspotPos)) %>%
  dplyr::mutate(isTpC = ifelse(isTpC == TRUE, 1, 0)) %>%
  dplyr::mutate(isAnyTpC = sum(isTpC)) %>%
  dplyr::ungroup()


# Identify Didimy

# Keep the same loop to do extra analysis
# TGAACA (AGTTCT); 010010, 010010
# Many loops are very specific
# GAAC (GTTC) are the most commonly mutated hotspots
# G, C == 1
# A, T == 0
# TGAACA (AGTTCT) == 010010 (010010)
# GAAC (GTTC)  == 1001 (1001)

# transform strings of ATCG to 0011
results_didimy <- results_TwinMutations %>%
  dplyr::mutate(loopSeq_binary = gsub('A|T','0', loopSeq)) %>%
  dplyr::mutate(loopSeq_binary = gsub('G|C','1', loopSeq_binary))


# Is the hotspot happening in 1001 or 101 context in the loop?
results_didimy <- results_didimy %>%
  dplyr::mutate(isSeq1001or101_up = ifelse(base::substr(loopSeq_binary, start = loopPos-2, stop = loopPos) == '101' |
                                             base::substr(loopSeq_binary, start = loopPos-3, stop = loopPos) == '1001', 1, 0),
                isSeq1001or101_down = ifelse(base::substr(loopSeq_binary, start = loopPos, stop = loopPos+2) == '101' |
                                               base::substr(loopSeq_binary, start = loopPos, stop = loopPos+3) == '1001', 1, 0)) %>%
  dplyr::mutate(has_1001or101 = ifelse(isSeq1001or101_up+isSeq1001or101_down == 0, "No", "Yes")) 


results_didimy <- results_didimy %>%
  dplyr::filter(base::grepl('[ACGT]\\[C>T', mutContext), nHotspotMutsInLoop >= 2, has_1001or101 == "Yes", isAnyTpC > 0) %>%
  dplyr::add_count(loopGenPos, name = "nMutPositionsInLoop") %>%
  dplyr::filter(nMutPositionsInLoop == 2) %>%
  dplyr::mutate(isDidimy = "Yes")


# export didymi data 
write.csv(results_didimy, file = "results/results_didimy.csv", row.names = FALSE)


