library(rhdf5)
library(GenomicAlignments)
library(GenomicFeatures)
library(ggplot2)

#'########################################################################################
# Read functions && Find overlaps                                                     ####
#'########################################################################################

read_and_findPolyA <- function(f5list, chunkNumber = 50){
  
  f5PathList <- f5list$fast5_path
  
  # Split data in chunks
  rows <- length(f5PathList)
  rowsPerChunk <- floor(rows/chunkNumber)
  extraLastChunk <- (rows %% chunkNumber)
  counter = 1
  
  result <- list()
  cores <- detectCores()
  
  for(chunk in c(1:chunkNumber)){
    
    if(chunk == chunkNumber)
      r <- rowsPerChunk + extraLastChunk
    else
      r <- rowsPerChunk
    
    f5Subset <- f5PathList[counter:(counter+r-1)]
    counter <- counter + r
    
    tmp <- mclapply(f5Subset, function(path2file){
        #Read event data
        tmpPath <- h5ls(path2file)[which(h5ls(path2file) == "/Analyses/Basecall_1D_000/BaseCalled_template")[1], 1]
        data.event <- h5read(path2file, tmpPath)$Events
        H5close()
        
        #Merge interesting information about the sample
        data <- list(
          eventMean = scale(data.event$mean),
          eventStart = data.event$start,
          eventMove = data.event$move,
          eventLength = data.event$length[1]
        )
        return(data)
    }, mc.cores = cores)
    
    result[[chunk]] <- findPolyA(tmp)
    
    cat("[INFO]  Chunk", chunk, "of", chunkNumber, "computed! \n")
    
  }
  
  # Unlist and make dataframe
  finalResult <- data.frame()
  for(l in result){
    df <- data.frame(l, stringsAsFactors = F)
    finalResult <- rbind(finalResult, df)
  }
  finalResult <- data.frame(f5list, finalResult)
  
  return(finalResult)
}

makeTranscriptIdTable <- function(hits, data){
  
  query <- from(hits)
  subjects <- to(hits)

  cores <- detectCores()
  result <- mclapply(unique(query), function(q){
    
    # Transcript ID
    result.idxId <- names(transcripts[q])
    
    # Get samples from data
    smpls <- subjects[query == q]
    result.subj <- paste(which(data$subjectHits %in% smpls), collapse = ",")
    smpls <- data[data$subjectHits %in% smpls,]
    
    # Number of Tags
    result.tags <- nrow(smpls)
    
    # Calculate mean
    result.mean <- mean(smpls$lengthPolyA_BasePair)
    
    # Calculate median
    result.median <- median(smpls$lengthPolyA_BasePair)
    
    # Result
    df <- data.frame(result.idxId, result.tags, result.mean, result.median, result.subj, stringsAsFactors = F)
    colnames(df) <- c("transcript_id", "tags", "mean_length", "median_length", "samples_comma_seperated")
    return(df)
    
  }, mc.cores = cores)
  
  # Unlist and make dataframe
  finalResult <- data.frame()
  for(l in result){
    df <- data.frame(l, stringsAsFactors = F)
    finalResult <- rbind(finalResult, df)
  }
  
  return(finalResult)
}

computeCorrelation <- function(data1, data2, mergeVec){
  tmp <- merge(data1, data2, by.x=mergeVec[1], by.y=mergeVec[2])
  colnames(tmp) <- c("transcript_id", "n", "b")
  tmp <- tmp[!is.nan(tmp$n),]
  tmp <- tmp[!is.na(tmp$b),]
  
  covar <- cor(tmp$n, tmp$b)
  return(covar)
}


#'########################################################################################
# Plot -- functions                                                                   ####
#'########################################################################################

plotHistogram <- function(data, column, labelVec, savePath = "Result/"){
  tmp <- as.data.frame(round(data[,column]))
  colnames(tmp) <- "l"
  
  plot <- ggplot(tmp, aes(x=l)) + geom_histogram(stat="bin", binwidth=1, fill="#ff5b5b") + 
    theme_minimal() + labs(title=labelVec[1], x=labelVec[2])
  ggsave(paste0(labelVec[1], ".png"), plot=plot, path=savePath, width=17.3, height=7.06, dpi=125)
  return(plot)
}

plotScatter <- function(data1, data2, mergeVec, labelVec, savePath = "Result/"){
  tmp <- merge(data1, data2, by.x=mergeVec[1], by.y=mergeVec[2])
  colnames(tmp) <- c("transcript_id", "n", "b")
  tmp <- tmp[!is.nan(tmp$n),]
  tmp <- tmp[!is.na(tmp$b),]
  
  plot <- ggplot(tmp, aes(x=n, y=b)) + 
    geom_line(data = data.frame(x=c(0:300)), mapping = aes(x=x, y=x), color="lightgrey", alpha = 0.4) +
    geom_point(color="steelblue", alpha = 0.4) + theme_minimal() +
    labs(title=labelVec[1], x=labelVec[2], y=labelVec[3])
  ggsave(paste0(labelVec[1], ".png"), plot=plot, path=savePath, width=17.3, height=7.06, dpi=125)
  return(plot)
}

plotDensity <- function(data1, data2, mergeVec, labelVec, classVec, savePath = "Result/"){
  tmp <- merge(data1, data2, by.x=mergeVec[1], by.y=mergeVec[2])
  colnames(tmp) <- c("transcript_id", "n", "b")
  tmp <- tmp[!is.nan(tmp$n),]
  tmp <- tmp[!is.na(tmp$b),]
  tmp$n <- round(tmp$n); tmp$b <- round(tmp$b)
  tmp <- rbind(data.frame(l=tmp$n, class=classVec[1]), data.frame(l=tmp$b, class=classVec[2]))
  
  plot <- ggplot(tmp, aes(x=l, fill=class)) + 
    geom_histogram(aes(y=..density..), alpha=0.5, position="identity", binwidth = 1) + 
    geom_density(alpha=0.4) + 
    labs(title=labelVec[1], x=labelVec[2])
  ggsave(paste0(labelVec[1], ".png"), plot=plot, path=savePath, width=17.3, height=7.06, dpi=125)
  return(plot)
}

# plotRaw <- function(x, start, end){
#   tmp <- data.frame(data[[x]]$rawSignal, c(1:(length(data[[x]]$rawSignal))))
#   colnames(tmp) <- c("raw", "idx")
#   tmp$PolyA <- 0
#   tmp$PolyA[start:end] <- 1
#   return(ggplot(tmp, aes(x=idx, y=raw, color=PolyA)) + geom_line() + theme_minimal())
# }

# plotEventMean <- function(x){
#   tmp <- data.frame(data[[x]]$eventStart, data[[x]]$eventMean)
#   colnames(tmp) <- c("start", "mean")
#   return(ggplot(tmp, aes(x=start, y=mean)) + geom_point() + geom_line(color="red"))
# }

# plotEventMove <- function(x){
#   tmp <- data.frame(data[[x]]$eventStart, as.numeric(data[[x]]$eventMove > 0))
#   colnames(tmp) <- c("start", "move")
#   return(ggplot(tmp, aes(x=start, y=move)) + geom_point())
# }

# plotEventReducedMove <- function(x){
#   tmp <- data.frame(data[[x]]$eventStart, as.numeric(data[[x]]$eventMove > 0))
#   colnames(tmp) <- c("start", "move")
#   tmp <- tmp[data[[x]]$eventP_mp_state > 0.8,]
#   tmp <- tmp[tmp$move > 0,]
#   return(ggplot(tmp, aes(x=start, y=move)) + geom_point())
# }

# plotEventP_mp_state <- function(x){
#   tmp <- data.frame(data[[x]]$eventStart, data[[x]]$eventP_mp_state)
#   colnames(tmp) <- c("start", "p_mp_state")
#   return(ggplot(tmp, aes(x=start, y=p_mp_state)) + geom_line( color = "darkgreen") + theme_minimal())
# }

# plotEventReducedP_mp_state <- function(x){
#   tmp <- data.frame(data[[x]]$eventStart, data[[x]]$eventP_mp_state, data[[x]]$eventMove)
#   colnames(tmp) <- c("start", "p_mp_state", "move")
#   tmp <- tmp[tmp$move > 0,]
#   return(ggplot(tmp, aes(x=start, y=p_mp_state)) + geom_line( color = "darkblue") + theme_minimal())
# }

# plotEventCombo <- function(x){
#   l <- data[[x]]
#   tmp <- data.frame(l$eventStart, l$eventMean, l$eventMove, l$eventP_mp_state, l$eventWeights)
#   colnames(tmp) <- c("dwelltime", "mean", "move", "p", "w")
#   return(ggplot(tmp, aes(x=dwelltime, y=mean)) + geom_line() + theme_minimal() +
#     geom_line(color="darkgreen", aes(y=move), alpha=0.6) +
#     #geom_line(color="darkblue", aes(y=w)) +
#     geom_line(color="darkblue", aes(y=p))
#   )
# }

# plotAndSaveAll <- function(x, path){
#   for(i in x){
#     polyA <- findPolyA(data, i)
#     polyA <- polyA*15
#     if(length(polyA) == 1) next()
#     if(polyA[1] < 0 || polyA[2] < 0 || polyA[3] < 0) next()
#     #plot <- plotEventCombo(i)
#     plot <- plotRaw(i, polyA[1], polyA[2])
#     ggsave(paste0(i, ".png"), plot=plot, path=path, width=17.3, height=7.06, dpi=75)
#   }
# }
