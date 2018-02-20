library(GenomicAlignments)
library(GenomicRanges)
library(GenomicFeatures)

#########################################################################################
# Import *.bam file, *.gtf file and then find overlaps
#########################################################################################

# Reading stuff
path_bam <- "path to bam"
path_gtf <- "path to gtf"
bam <- readGAlignments(path_bam, use.names = T)
gtf = makeTxDbFromGFF(path_gtf, format="gtf")

# Making reference transcripts
transcripts <- transcriptsBy(gtf, "gene")
# transcripts <- exonsBy(gtf, by="tx", use.names = T) # Making transcripts exon ???
# transcripts <- transcriptsBy(gtf, "cds")
# transcripts <- transcriptsBy(gtf, "exon")
# transcripts <- extractTranscriptSeqs(gtf, use.names = T) # ??? is not working, no genome, how to import *.fa as BSgenome?

# Finding overlaps
hits = findOverlaps(transcripts, bam) # hits object containing all hits
# hits = countOverlaps(transcripts, foo) # Just an Integer array, need hits-object

# Seperate all overlaps by name
traceNames <- bam[to(hits)]@NAMES # all seq from bam that have hits with transcripts

# Other stuff that may be useful?
transcripts[from(hits)] # all transcripts that have hits with bam


#########################################################################################
# Traceback to *.fast5 files that overlaps 
#########################################################################################

find_fast5_filename <- function(filepath_fastq, traceNames) {
  t <- traceNames
  f <- file(filepath_fastq, "r")
  fast5_filenames <- vector(mode = "character", length = 0)
  while ( TRUE ) {
    line = readLines(f, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if(substring(line, 1, 1) == "@"){
      # compute match between fastq file and trace list
      line <- substring(line, 2)
      splitLine <- unlist(strsplit(line, " "))
      pos <- match(splitLine[1], t)
      # Check if no match
      if(is.na(pos)) next
      # If match, add filename to file-list
      fast5_filenames[length(fast5_filenames)+1] <- splitLine[2]
      # Delete all entries for this match at the trace list for better performans
      t <- t[! t %in% splitLine[1]]
    }
  }
  close(f)
  return(fast5_filenames)
}
