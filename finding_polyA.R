library(rhdf5)
library(ggplot2)

path_0 <- "../Data/0/"
f5List <- paste0(path_0, f5FileList)

##########################################################################################
# Extract raw from *.fast5 and plot the data
##########################################################################################

#' Get raw data from one *.fast5 file
#' 
#' @param path2file The path to the nanopore *.fast5 file to read in
#' @return Dataframe containing raw data and time
extractRaw <- function(path2file){
  tmpPath <- h5ls(path2file)[(which(h5ls(path2file) == "/Raw/Reads") + 1) ,1]
  data.raw <- h5read(path2file, tmpPath)$Signal
  
  time <- c(1:length(data.raw))
  #time <- time/4000 # Time in seconds
  df <- data.frame(data.raw, time)
  colnames(df) <- c("Raw", "Time")
  return(df)
}

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

#' #' Plot all raw data from given *.fast5 filelist. The raw data will be plotted on 
#' #' top of each other.
#' #' (The plot you get is just chaos)
#' #'
#' #' @param f5FileList List over *.fast5 that you want to plot the raw from.
#' #' @return The ggplot2 plot that shows the data on top of each other
#' plotAllRaw <- function(f5FileList){
#'   plot <- ggplot()
#'   #colors <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", 
#'   #"#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", 
#'   #"#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", 
#'   #"#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")
#' 
#'   for(f5 in f5FileList){
#'     rawData <- extractRaw(f5)
#'     rawData.stats <- calcStatOfRaw(rawData)
#'     if(length(rawData$Raw) > 20000)
#'       plot <- plot + geom_line(data = rawData.stats[1:(20000/10),], aes(x=Time, y=Mean))
#'     else
#'       plot <- plot + geom_line(data = rawData.stats, aes(x=Time, y=Mean))
#'   }
#'   return(plot)
#' }

# Other plot
# ggplot(raw.plot, aes(x=Time)) + geom_line(aes(y=Raw, colour="Raw")) + geom_line(aes(y=Mean, colour="Mean")) + scale_color_manual(values=c("blue", "#FE2E2E"))


##########################################################################################
# Math function & Find PolyA
##########################################################################################

#' Calculate kind of "running" mean with a alpha value to simplyfy the raw data. 
#' 
#' @param rawData Dataframe containing raw data you want to process
#' @param alphaV The alpha value let you control how much you want to smooth out the mean
#' @return Dataframe containing mean, reverse mean and time on the data.
calcMeanOfRaw <- function(rawData, alphaV = 0.05){
  raw <- rawData$Raw
  revRaw <- rev(raw)
  
  meanA <- vector(mode = "double", length = 0)
  meanB <- vector(mode = "double", length = 0)
  
  meanA[1] <- raw[1]
  meanA[2] <- mean(raw[1:2])
  meanB[1] <- revRaw[1]
  meanB[2] <- mean(revRaw[1:2])
  
  for(value in c(2:(length(raw)-1))){
    meanA[length(meanA)+1] <- meanA[length(meanA)] * (1-alphaV) + 
      mean(raw[value:(value+1)]) * alphaV
    meanB[length(meanB)+1] <- meanB[length(meanB)] * (1-alphaV) + 
      mean(revRaw[value:(value+1)]) * alphaV
  }
  df <- data.frame(meanA, rev(meanB), rawData$Time)
  colnames(df) <- c("Mean", "MeanB", "Time")
  return(df)
}


#' Run over data with a scope and calculate some stats for every scope 
#' (mean, max, min, max-min)
#' 
#' @param scopeSize The size of the scope
#' @param rawData Dataframe containing the raw data you want to process
#' @return Dataframe containing the calculated stats
calcStatOfRaw <- function(rawData, scopeSize = 10){
  raw <- rawData$Raw
  
  meanOfScopes <- vector(mode = "double", length = 0)
  maxOfScopes <- vector(mode = "double", length = 0)
  minOfScopes <- vector(mode = "double", length = 0)
  varOfScopes <- vector(mode = "double", length = 0)
  
  idx = 1
  while(!FALSE){
    idx2 = idx + scopeSize
    if(idx2 >= length(raw)) break()
    
    meanOfScopes[length(meanOfScopes)+1] <- mean(raw[idx:idx2])
    maxOfScopes[length(maxOfScopes)+1] <- max(raw[idx:idx2])
    minOfScopes[length(minOfScopes)+1] <- min(raw[idx:idx2])
    varOfScopes[length(varOfScopes)+1] <- maxOfScopes[length(maxOfScopes)] - 
      minOfScopes[length(minOfScopes)]
    
    idx <- idx2
  }
  
  timeOfScopes <- 1:length(varOfScopes)
  timeOfScopes <- timeOfScopes * scopeSize
  
  stats <- data.frame(meanOfScopes, varOfScopes, maxOfScopes, minOfScopes, timeOfScopes)
  colnames(stats) <- c("Mean", "Var", "Max", "Min", "Time")
  return(stats)
}

#' Find polyA tail (seems working correctly in most cases, but I'm still working on it)
#' 
#' @param rawData Dataframe containing the rawData you want to find a PolyA on. 
#' @return Dataframe containing rawData with additional PolyA column (binary, 0 or 1)
findPolyA <- function(rawData){
  rawData$PolyA <- 0
  scopeSize <- 10
  
  # Make sure data is normalized and compute stats
  rawData$Raw <- scale(rawData$Raw)
  stat <- calcStatOfRaw(rawData, scopeSize)
  
  # I assume a polyA at specific height (maybe with some threshold).
  # In normalized data e.g. between 0 and 1,5
  hTable <- rle(as.numeric(( as.numeric(stat$Mean > 0)) & (as.numeric(stat$Mean < 1.5))))
  hTable$lengths <- hTable$lengths * scopeSize
  len <- length(hTable$lengths)
  
  # Looping throw data, looking for PolyA
  idx <- 1
  while((idx + 5) <= len){
    possiblePolyA <- FALSE
    
    # Before polyA, I assume that the raw seq is lower than 0 over a long time threshold
    # e.g. about 950 seems to work great. 
    if(hTable$values[idx] == 0 && hTable$lengths[idx] > 950)
      possiblePolyA <- TRUE
    
    # Check for polyA
    if(possiblePolyA){
      sIdx <- idx
      idx <- idx + 1
      
      # Check for noisy beginning of PolyA with a small pattern-threshold of e.g. 30 and
      # with stepsize/range of 2 (idx+1 and idx+3). I assume a minimum length of PolyA 
      # with e.g. about 200 steps in data
      if(hTable$lengths[idx] <= 200){
        possiblePolyA <- FALSE
        if(max(as.numeric(hTable$lengths[c((idx+2), (idx+4))] > 200)) == 1){
          tmpIdx <-  2 * which(hTable$lengths[c((idx+2), (idx+4))] > 200)[1] + idx
          tmpC <- c((idx+1), (idx+3))
          if(tmpIdx < (idx+3))
            tmpC <- tmpC[1]
          
          if(max(as.numeric(hTable$lengths[tmpC] > 30)) == 0){
            possiblePolyA <- TRUE
            idx <- tmpIdx #add
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
        
        # Mark PolyA in data and terminate loop
        rawData$PolyA[sPolyA:ePolyA] <- 1
        break()
      }
    }
    idx <- idx + 1
  }
  
  # TODO: I assume that polyA begins in the lower region growing up to over 0  so I want 
  # to find the exact beginning by going backwards from polyA in data, checking height of 
  # datapoints.
  
  # TODO: At the end of polyA I assume that max-min for stat scopes changes.
  # Mean of data also changes as the sample will goes up and down more crazy
  # Better end of polyA can be detected by checking if a data points are near the mean
  # of the polyA with some threshold
  
  View(data.frame(hTable$values, hTable$lengths))
  
  return(rawData)
}
