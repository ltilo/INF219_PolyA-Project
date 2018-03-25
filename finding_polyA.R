library(rhdf5)
library(ggplot2)
library(gridExtra)

path_0 <- "../Data/0/"
f5List <- paste0(path_0, f5FileList)

##########################################################################################
# Extract raw from *.fast5 and plot data functions
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

#' This function will read in all given raw data, find polyA, mark polyA and plot the 
#' data. The plot will be saved to a subfolder called "Plots". 
#' Speed is ca. 60 plots/saved-images per minute on my intel i3 mashine. 
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
#'     rawData$Raw <- scale(rawData$Raw)
#'     rawData.stats <- calcMeanOfRaw(rawData, 10)
#'     rawData.stats$Time <- rawData.stats$Time/(10*length(rawData.stats$Time))
#'     
#'     plot <- plot + geom_line(data = rawData.stats, aes(x=Time, y=Mean))
#'   }
#'   return(plot)
#' }

# Other plot
# ggplot(raw.plot, aes(x=Time)) + geom_line(aes(y=Raw, colour="Raw")) + geom_line(aes(y=Mean, colour="Mean")) + scale_color_manual(values=c("blue", "#FE2E2E"))


##########################################################################################
# Math function & Find PolyA
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
#' @return Dataframe containing start, end and length of PolyA. The function
#' marks the PolyA with an additional column in raw data.
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

#' This function reads the given raw data, finds the polyA for every raw and return it.
#' 
#' @param f5List Charactervector containing path to all *.fast5 you want to process
#' @return dataframe containing information about start, end and length of PolyA 
statsPolyA <- function(f5List){
  polyA_data <- data.frame(
    PolyA_Start = integer(),
    PolyA_End = integer(),
    PolyA_Length = integer()
  )
  
  for(f5 in f5List){
    raw <- extractRaw(f5)
    polyA <- findPolyA(raw)
    polyA_data <- rbind(polyA_data, polyA)
  }
  
  return(polyA_data)
}