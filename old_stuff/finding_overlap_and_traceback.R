library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)

#########################################################################################
# Import *.bam file, *.gtf file and then find overlaps
#########################################################################################

# Reading stuff
#path_bam <- "../Data/shield_pass_0.bam"
#path_gtf <- "../Data/Bartel/Danio_rerio.Zv9.79.gtf"
#path_fastq <- "../Data/shield_pass_0.fastq"

bam <- readGAlignments(path_bam, use.names = T)
gtf = makeTxDbFromGFF(path_gtf, format="gtf")

# Making reference transcripts
transcripts <- exonsBy(gtf, by = "tx", use.names=T)

# Filter only coding genes
transcripts <- transcripts[names(cdsBy(gtf, by = "tx", use.names = T))]

# Finding overlaps
hits <- findOverlaps(transcripts, bam) # hits object containing all hits

# Seperate all overlaps by name, all seq from bam that have hits with transcripts
traceNames <- unique(names(bam[to(hits)]))


#########################################################################################
# Traceback to *.fast5 files that overlaps 
#########################################################################################

#' This function finds the names of the corresponding *.fast5 files based on the fastQ
#' file. 
#' 
#' @param path_fastq String, Path to the fastq file containing all filenames
#' @param traceNames Character vector containing names that overlaps to reference
#' @return Dataframe containing *.fast5 filenames that matches with traceNames
find_fast5_filenames <- function(path_fastq, traceNames) {
  t <- traceNames
  f <- file(path_fastq, "r")
  fast5_filenames <- data.frame(
    fast5_filename = character(),
    traceName = character()
  )
  
  while(!FALSE) {
    line = readLines(f, n = 1)
    if ( length(line) == 0 ) break
    if(substring(line, 1, 1) == "@"){
      # compute match between fastq file and trace list
      line <- substring(line, 2)
      splitLine <- unlist(strsplit(line, " "))
      pos <- match(splitLine[1], t)
      # Check if no match
      if(is.na(pos)) next
      # If match, add filename to file-list
      #fast5_filenames[length(fast5_filenames)+1] <- paste0(splitLine[2], ".fast5")
      
      df <- data.frame(paste0(splitLine[2], ".fast5"), splitLine[1])
      colnames(df) <- c("fast5_filename", "traceName")
      
      fast5_filenames <- rbind(fast5_filenames, df)
    }
  }
  close(f)
  return(fast5_filenames)
}

# Path to first 4000 data
#path_0 <- "../Data/0/"
f5List <- find_fast5_filenames(path_fastq, traceNames)
f5List$fast5_path <- paste0(path_0, f5List$fast5_filename)