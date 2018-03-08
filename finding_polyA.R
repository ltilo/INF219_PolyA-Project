library(rhdf5)
library(ggplot2)

path_0 <- "../Data/0/"
f5List <- paste0(path_0, f5FileList)

########################################################################################
# Extract raw from *.fast5 and plot the data
########################################################################################

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
#' @param rawData Dataframe containing raw data and time you want to plot, additional PolyA column
#' @param showPolyA True if you want to mark PolyA in the plot, otherwise false (column PolyA must exist)
#' @return The ggplot2 plot that is plotting the data
plotRaw <- function(rawData, showPolyA = FALSE){
  if(showPolyA == TRUE && !is.null(rawData$PolyA))
    return(ggplot(rawData, aes(Time, Raw, colour = PolyA)) + geom_line())
  else
    return(ggplot(rawData, aes(Time, Raw)) + geom_line())
}

#' Plot all raw data from given *.fast5 filelist. The raw data will be plotted on top of each other.
#' (The plot you get is just chaos, TODO: normalization of time, is this even possible?)
#' 
#' @param f5FileList List over *.fast5 that you want to plot the raw from.
#' @return The ggplot2 plot that shows the data on top of each other
plotAllRaw <- function(f5FileList){
  plot <- ggplot()
  colors <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
  "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
  "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
  "#8A7C64", "#599861")
  
  for(f5 in f5FileList){
    rawData <- extractRaw(f5)
    rawData.stats <- calcMeanOfRaw(rawData)
    if(length(rawData$Raw) > 20000)
      plot <- plot + geom_line(data = rawData.stats[1:20000,], aes(x=Time, y=MeanB), color=sample(colors, 1))
    else
      plot <- plot + geom_line(data = rawData.stats, aes(x=Time, y=MeanB), color=sample(colors, 1))
  }
  return(plot)
}

# Other plot
# ggplot(raw.plot, aes(x=Time)) + geom_line(aes(y=Raw, colour="Raw")) + geom_line(aes(y=Mean, colour="Mean")) + scale_color_manual(values=c("blue", "#FE2E2E"))


########################################################################################
# Math function & Find PolyA
########################################################################################

#' Calculate kind of "running" mean with a alpha value to simplyfy the raw data. 
#' 
#' @param rawData Dataframe containing raw data you want to process
#' @param alphaV The alpha value that let you control how much you want to smooth out the mean
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
    meanA[length(meanA)+1] <- meanA[length(meanA)] * (1-alphaV) + mean(raw[value:(value+1)]) * alphaV
    meanB[length(meanB)+1] <- meanB[length(meanB)] * (1-alphaV) + mean(revRaw[value:(value+1)]) * alphaV
  }
  df <- data.frame(meanA, rev(meanB), rawData$Time)
  colnames(df) <- c("Mean", "MeanB", "Time")
  return(df)
}


#' Run over data with a scope and calculate some stats for every scope (mean, max, min, max-min)
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
    varOfScopes[length(varOfScopes)+1] <- maxOfScopes[length(maxOfScopes)] - minOfScopes[length(minOfScopes)]
    
    idx <- idx2
  }
  
  timeOfScopes <- 1:length(varOfScopes)
  timeOfScopes <- timeOfScopes * scopeSize
  
  stats <- data.frame(meanOfScopes, varOfScopes, maxOfScopes, minOfScopes, timeOfScopes)
  colnames(stats) <- c("Mean", "Var", "Max", "Min", "Time")
  return(stats)
}

#' Find polyA tail (may working correctly in some case but still working)
#' 
#' @param rawData Dataframe containing the rawData you want to find a PolyA on. 
#' @return Dataframe containing rawData with additional PolyA column (binary, 0 or 1)
findPolyA <- function(rawData){
  rawData$PolyA <- 0
  scopeSize <- 10
  stat <- calcStatOfRaw(rawData, scopeSize)
  
  # I except a polyA at specific height, between e.g. 700 and 900 (maybe with some threshold)
  hTable <- rle(as.numeric(( as.numeric(stat$Mean > 700)) & (as.numeric(stat$Mean < 900))))
  
  # I except not a polyA right at beginning, so there is some time threshold, e.g. 2000 in raw data
  idx <- 1
  len <- length(hTable$lengths)
  while(idx <= len){
    if(sum(hTable$lengths[1:idx]) > (2000/scopeSize)){
      idx <- idx + 1
      break()
    }
    idx <- idx + 1
  }
  
  # Check for polyA, given a minimum length e.g 200 steps in raw data
  startPolyA <- -1
  endPolyA <- -1
  while(idx <= len){
    if(hTable$values[idx] == 1){
      if(hTable$lengths[idx] > (200/scopeSize)){
        startPolyA <- sum(hTable$lengths[1:(idx-1)]) * scopeSize
        endPolyA <- startPolyA + (hTable$lengths[idx] * scopeSize)
        rawData$PolyA[startPolyA:endPolyA] <- 1
      }
    }
    idx <- idx + 1
  }
  
  # Mark polyA if existing
  #if(startPolyA > 0 && endPolyA > 0)
  #  rawData$PolyA[startPolyA:endPolyA] <- 1
  
  return(rawData)
}
