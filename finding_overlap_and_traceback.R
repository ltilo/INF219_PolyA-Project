library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)

#########################################################################################
# Import *.bam file, *.gtf file and then find overlaps
#########################################################################################

# Reading stuff
path_bam <- "../Data/shield_pass_0.bam"
path_gtf <- "../Data/Bartel/Danio_rerio.Zv9.79.gtf"
path_fastq <- "../Data/shield_pass_0.fastq"

bam <- readGAlignments(path_bam, use.names = T)
gtf = makeTxDbFromGFF(path_gtf, format="gtf")

# Making reference transcripts
transcripts <- exonsBy(gtf, by = "tx")

# TODO Filter only coding genes


# Finding overlaps
hits = findOverlaps(transcripts, bam) # hits object containing all hits
# hits = countOverlaps(transcripts, foo) # Just an Integer array, need hits-object

# Seperate all overlaps by name, all seq from bam that have hits with transcripts
traceNames <- names(bam[to(hits)])

# Other stuff that may be useful?
#transcripts[from(hits)] # all transcripts that have hits with bam

#########################################################################################
# Traceback to *.fast5 files that overlaps 
#########################################################################################

#' This function finds the names of the corresponding *.fast5 files
#' 
#' @param path_fastq String, Path to the fastq file containing all filenames
#' @param traceNames Character vector containing traceNames that overlaps to reference
#' @return Character vector containing *.fast5 filenames that matches traceNames
find_fast5_filenames <- function(path_fastq, traceNames) {
  t <- traceNames
  f <- file(path_fastq, "r")
  fast5_filenames <- vector(mode = "character", length = 0)
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
      fast5_filenames[length(fast5_filenames)+1] <- paste0(splitLine[2], ".fast5")
      # Delete all entries for this match at the trace list for better performance
      t <- t[! t %in% splitLine[1]]
    }
  }
  close(f)
  return(fast5_filenames)
}

f5FileList <- find_fast5_filename(path_fastq, traceNames)

