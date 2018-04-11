library(rhdf5)
library(ggplot2)

##########################################################################################
# Extract data from *.fast5 -- functions
##########################################################################################

#' Get raw data from one *.fast5 file
#' 
#' @param path2file The path to the nanopore *.fast5 file to read in
#' @return Dataframe containing raw data and time
extractRaw <- function(path2file){
  tmpPath <- h5ls(path2file)[(which(h5ls(path2file) == "/Raw/Reads") + 1) ,1]
  data.raw <- h5read(path2file, tmpPath)$Signal
  H5close()
  
  time <- c(1:length(data.raw))
  #time <- time/4000 # Time in seconds
  df <- data.frame(data.raw, time)
  colnames(df) <- c("Raw", "Time")
  return(df)
}

#' This function will extract and return the actual dwellTime per base
#' 
#' @param path2fil The path to the nanopore *.fast5 file to extract
#' @param dwellTimePolyAEnd The dwellTime where the polyA ends in the sample
#' @return Numeric value that is representing the dwellTime for one base
extractDwellTimePerBP <- function(path2file, dwellTimePolyAEnd){
  tmpPath <- h5ls(path2file)[which(h5ls(path2file) == "/Analyses/Basecall_1D_000/BaseCalled_template")[1], 1]
  data.event <- h5read(path2file, tmpPath)$Events
  H5close()
  
  baseIdx <- which(data.event$start >= dwellTimePolyAEnd)[1]
  len <- length(data.event$move)
  return((tail(data.event$start, n=1)+tail(data.event$length, n=1)-data.event$start[baseIdx])/sum(data.event$move[baseIdx:len]))
}

extractFastqLength <- function(path2file, dwellTimePolyAEnd){
  tmpPath <- h5ls(path2file)[which(h5ls(path2file) == "/Analyses/Basecall_1D_000/BaseCalled_template")[1], 1]
  data.event <- h5read(path2file, tmpPath)$Fastq
  H5close()
  fastq <- strsplit(data.event, "\n")[[1]][2]
  
  # Should I cut the fastq such that it matches witch moves?
  #baseIdx <- which(data.event$Fastq >= dwellTimePolyAEnd)[1]
  #m <- sum(data.event$move[1:baseIdx])
  #return(nchar(fastq) - m)
  
  return(nchar(fastq))
}

##########################################################################################
# Plot functions
##########################################################################################

#' Plot raw data as graph. Plots also PolyA if the rawData contains PolyA column.
#' 
#' @param rawData Dataframe containing raw data and time you want to plot, additional 
#' PolyA column
#' @param showPolyA True if you want to mark PolyA in the plot, otherwise false 
#' (column PolyA must exist)
#' @return The ggplot2 plot that is plotting the data
plotRaw <- function(rawData, showPolyA = FALSE){
  if(showPolyA == TRUE && !is.null(rawData$PolyA))
    return(ggplot(rawData, aes(Time, Raw, colour = PolyA)) + geom_line())
  else
    return(ggplot(rawData, aes(Time, Raw)) + geom_line())
}

#' This function will read in all given raw data, find polyA, mark polyA and plot the 
#' data. The plot will be saved to a subfolder called "Plots". 
#' NB! Speed is ca. 60 plots/saved-images per minute on my intel i3 mashine. 
#' 
#' @param path2File Charactervector containing all path to raw *.fast5 files
#' @param path_plot The path to the folder where you want to save the plots/images
plotAndSaveAll <- function(path2File, path_plot){
  idx <-  1;
  while(idx <= length(path2File)){
    raw <- extractRaw(path2File[idx])
    polyA <- findPolyA(raw)
    raw <- markPolyA(raw, polyA)
    plot <- plotRaw(raw, T)
    ggsave(paste0(idx, ".png"), plot = plot, path = path_plot, width = 17.3, height = 7.06, dpi=75)
    idx <- idx + 1;
  }
}

# Other plot
# ggplot(raw.plot, aes(x=Time)) + geom_line(aes(y=Raw, colour="Raw")) + geom_line(aes(y=Mean, colour="Mean")) + scale_color_manual(values=c("blue", "#FE2E2E")

##########################################################################################
# Math functions & Find PolyA functions
##########################################################################################

#' Run over data with a scope and calculate mean for every scope
#' 
#' @param scopeSize The size of the scope
#' @param rawData Dataframe containing the raw data you want to process
#' @return Dataframe containing the calculated mean stats
calcMeanOfRaw <- function(rawData, scopeSize = 10){
  raw <- rawData$Raw
  meanOfScopes <- vector(mode = "double", length = 0)
  
  idx = 1
  while(!FALSE){
    idx2 = idx + scopeSize
    if(idx2 >= length(raw)) break()
    meanOfScopes[length(meanOfScopes)+1] <- mean(raw[idx:idx2])
    idx <- idx2
  }
  
  timeOfScopes <- 1:length(meanOfScopes)
  timeOfScopes <- timeOfScopes * scopeSize
  
  stats <- data.frame(meanOfScopes, timeOfScopes)
  colnames(stats) <- c("Mean", "Time")
  return(stats)
}

#' Find polyA tail (seems working correctly in most cases, needs fine-tuning)
#' 
#' @param rawData Dataframe containing the rawData you want to find a PolyA on.
#' @return Dataframe containing start, end and length of PolyA in DwellTime. 
findPolyA <- function(rawData){
  scopeSize <- 10; sPolyA <- 0; ePolyA <- 0;
  
  # Make sure data is normalized and compute stats
  rawData$Raw <- scale(rawData$Raw)
  stat <- calcMeanOfRaw(rawData, scopeSize)
  
  # I assume a polyA at specific height (maybe with some threshold).
  # In normalized data e.g. between 0.2 and 1,5
  hTable <- rle(as.numeric((as.numeric(stat$Mean > 0.2)) & (as.numeric(stat$Mean < 1.5))))
  hTable$lengths <- hTable$lengths * scopeSize
  len <- length(hTable$lengths)
  
  # Looping throw data, looking for PolyA
  idx <- 1
  while((idx + 5) <= len){
    
    # Before polyA, I assume that the raw seq is lower than 0 over a long time threshold
    # e.g. about 950 seems to work great. 
    if(hTable$values[idx] == 0 && hTable$lengths[idx] > 950){
      possiblePolyA <- TRUE
      sIdx <- idx
      idx <- idx + 1
        
      # Check for noisy beginning of PolyA with a small pattern-threshold of e.g. 100 and
      # with stepsize/range of 2 (idx+1 and idx+3). I assume a minimum length of PolyA 
      # with e.g. about 200 steps in data
      if(hTable$lengths[idx] <= 200){
        possiblePolyA <- FALSE
        if(max(as.numeric(hTable$lengths[c((idx+2), (idx+4))] > 200)) == 1){
          tmpIdx <-  2 * which(hTable$lengths[c((idx+2), (idx+4))] > 200)[1] + idx
          tmpC <- c((idx+1), (idx+3))
          if(tmpIdx < (idx+3))
            tmpC <- tmpC[1]
            
          if(max(as.numeric(hTable$lengths[tmpC] > 100)) == 0){
            possiblePolyA <- TRUE
            idx <- tmpIdx
          }
        }
      }
        
      # Find end of PolyA
      if(possiblePolyA){
        sPolyA <- sum(hTable$lengths[1:sIdx])
          
        # PolyA may have drops, detecting these by a threshold e.g. 150
        while((idx+2) < len){
          if(hTable$lengths[idx+1] < 150 && hTable$lengths[idx+2] > 150)
            idx <- idx + 2
          else
            break()
        }
          
        ePolyA <- sum(hTable$lengths[1:idx])
        break()
      }
    }
    idx <- idx + 1
  }
  
  # Prepare result and return it
  result <- data.frame(sPolyA, ePolyA, ePolyA-sPolyA)
  colnames(result) <- c("PolyA_Start", "PolyA_End", "PolyA_Length")
  return(result)
}

#' This functions marks the polyA in the raw data
#' 
#' @param rawData The raw data you want to mark the PolyA on
#' @param polyA The dataframe returned from findPolyA() or one row from statsPolyA() 
#' conatining information about start and end of the PolyA
#' @return Raw data with additional polyA column, marking the polyA
markPolyA <- function(rawData, polyA){
  rawData$PolyA <- 0
  rawData$PolyA[polyA$PolyA_Start:polyA$PolyA_End] <- 1
  return(rawData)
}

#' This function reads in the filelist, finds the polyA for every raw and converts
#' dwellTime to BasePair
#' NB! Runtime is slow, 5 minutes for 2909 samples on my intel i3 mashine
#' 
#' @param f5List dataframe from find_fast5_filenames() containing name and path to fast5
#' files you want to process
#' @return dataframe containing information about start, end and length of PolyA 
statsPolyA <- function(f5List){
  polyA_data <- data.frame(
    PolyA_Start = integer(),
    PolyA_End = integer(),
    PolyA_Length = integer(),
    DwellTimePerBP = numeric(),
    Fastq_Length = numeric()
  )
  
  percent <- 0;
  for(f5 in f5List$fast5_path){
    raw <- extractRaw(f5)
    lraw <- length(raw$Raw)
    polyA <- findPolyA(raw)
    dwellPerBP <- extractDwellTimePerBP(f5, polyA$PolyA_End)
    fastq <- extractFastqLength(f5, polyA$PolyA_End)
    
    # Shifting start and end position and computing length in bp
    len <- polyA$PolyA_Length/dwellPerBP
    st <- (lraw - polyA$PolyA_End)/dwellPerBP
    en <- (lraw - polyA$PolyA_Start)/dwellPerBP 
    
    df <- data.frame(st, en, len, dwellPerBP, polyA$PolyA_Length, fastq)
    colnames(df) <- c("PolyA_Start", "PolyA_End", "PolyA_Length", "DwellTimePerBP", 
                      "PolyA_Length_DwellTime", "Fastq_Length")
    
    polyA_data <- rbind(polyA_data, df)
    
    if(percent %% 10 == 0)
      cat(c(percent, " samples computed! \n"))
    
    percent <- percent + 1
  }
  
  df <- data.frame(f5List[, c(3,1,2)], polyA_data, stringsAsFactors = FALSE)
  return(df)
}
