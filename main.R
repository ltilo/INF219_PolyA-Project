#########################################################################################
# Main R-Script -- NB! Takes ca. 8 minutes
#########################################################################################

# Set all path to first 4k Reads of nanopore data
path_bam <- "../Data/shield_pass_0.bam"
path_fastq <- "../Data/shield_pass_0.fastq"
path_0 <- "../Data/0/"

# Set all path to Bartel data (Shield - 6hpf)
path_gtf <- "../Data/Bartel/Danio_rerio.Zv9.79.gtf"
path_b.dmv <- "../Data/Bartel/PAL_seq/dist_median_values_subset_of_standards_zv9/Dre_mock_6hpf_PAL_length_dist_only_10_50_100_lengths.csv"
path_b.dav <- "../Data/Bartel/PAL_seq/dist_all_values_all_standards_zv9/Dre_mock_6hpf_PAL_length_dist_all_lengths.csv"
path_b.dfp <- "../Data/Bartel/PAL_seq/dist_from_publication/GSE52809_Dre_mock_6hpf.txt"

# Find overlap and traceback to fast5 files (NB! ca. 1 minutes on my intel i3)
print("[INFO] Finding overlap and traceback to fast5 files")
source('finding_overlap_and_traceback.R')

# Read nanopore data and find polyA (NB! ca. 5 minutes on my intel i3)
print("[INFO] Finding polyA")
source('finding_polyA.R')
polyA <- statsPolyA(f5List)

# Read in Bartel data
b.dmv <- read.csv(path_b.dmv, header = T, sep = '\t', stringsAsFactors = F)
b.dav <- read.csv(path_b.dav, header = T, sep = '\t', stringsAsFactors = F)
b.dfp <- read.csv(path_b.dfp, header = T, sep = '\t', stringsAsFactors = F)

# Compute new hits list with new idx (NB! ca. 1 minute on my intel i3)
subject <- vector(mode = "integer", length = 0)
for(sampl in to(hits)){
  newIdx <- which(polyA$traceName == names(bam[sampl]))
  subject[length(subject)+1] <- newIdx 
}
newHits <- data.frame(from(hits), subject)
colnames(newHits) <- c("queryHits", "subjectHits")
rm(newIdx, sampl, subject)

# Make the table of tables for nanopore data (NB! ca. 1 minute on my intel i3)
nanopore <- data.frame(
  transcript_id = character(),
  fast5_files.comma_seperated = character(),
  tail_length.comma_seperated = character(),
  tag_count = integer(),
  mean_length = numeric(),
  median_length = numeric(),
  mean_dwell_length = numeric()
)
for(sampl in unique(newHits$queryHits)){
  txId <- names(transcripts[sampl])
  allSubj <- newHits$subjectHits[which(newHits$queryHits == sampl)]
  allLength <- polyA$PolyA_Length[allSubj]
  allSubj <- polyA$fast5_filename[allSubj]
  median_length <- median(allLength)
  mean_length <- mean(allLength)
  mean_dwell_length <- mean(polyA$PolyA_Length_DwellTime[allSubj])
  
  df <- data.frame(
    txId, paste(allSubj, collapse = ","), paste(allLength, collapse = ","),
    length(allLength), mean(allLength), median(allLength), mean_dwell_length, 
    stringsAsFactors = F
  )
  colnames(df) <- colnames(nanopore)
  nanopore <- rbind(nanopore, df)
}
rm(df, sampl, txId, allSubj, allLength, median_length, mean_length, mean_dwell_length)

## Plot the requested plots. 
# Nanopore length distr.
tmp <- as.data.frame(round(nanopore$mean_length))
colnames(tmp) <- "l"
ggplot(tmp, aes(x=l)) + geom_histogram(stat="bin", binwidth=1, fill="#ff5b5b") + theme_minimal() +
  labs(title="Nanopore polyA length", x="polyA length in bp")
rm(tmp)

# Nanopore dwelltime per dp
tmp <- as.data.frame(round(polyA$DwellTimePerBP))
colnames(tmp) <- "d"
ggplot(tmp, aes(x=d)) + geom_histogram(stat="bin", binwidth=1, fill="#ff5b5b") + theme_minimal() +
  labs(title="Nanopore - Dwelltime for one basepair", x="DwellTime for one basepair")
rm(tmp)

# Nanopore vs bartel (dist_median_values_subset_of_standards_zv9)
tmp <- merge(nanopore[,c(1,5)], b.dmv[,c(1,4)], by.x="transcript_id", by.y="X.transcript_id")
# tmp <- merge(nanopore[,c(1,7)], b.dmv[,c(1,4)], by.x="transcript_id", by.y="X.transcript_id")
colnames(tmp) <- c("transcript_id", "n", "b")
ggplot(tmp, aes(x=n, y=b)) + geom_line(data = data.frame(x=c(0:150)), mapping = aes(x=x, y=x), color="grey") +
 geom_point(color="steelblue") + theme_minimal() +
 labs(title="PolyA length (nanopore vs dist_median_values_subset_of_standards_zv9)",
      x="Nanopore polyA length", y="dist_median_values_subset_of_standards_zv9 polyA length")

# ggplot(tmp, aes(x=n, y=b)) + geom_point(color="steelblue") + theme_minimal() +
#   labs(title="PolyA length (nanopore vs dist_median_values_subset_of_standards_zv9)",
#        x="Nanopore polyA length in dwelltime", y="dist_median_values_subset_of_standards_zv9 polyA length")

n_vs_b.dmv_cor <- cor(tmp$n, tmp$b)
n_vs_b.dmv_N <- length(tmp$n)

tmp$n <- round(tmp$n); tmp$b <- round(tmp$b)
tmp <- rbind(data.frame(l=tmp$n, class="nanopore"), data.frame(l=tmp$b, class="dmv"))
ggplot(tmp, aes(x=l, fill=class)) + geom_histogram(aes(y=..density..), alpha=0.5, 
     position="identity", binwidth = 1) + geom_density(alpha=0.4) + 
  labs(title="PolyA length density (dist_median_values_subset_of_standards_zv9 polyA length)",
       x="polyA length in bp")
rm(tmp)

# Nanopore vs bartel (dist_all_values_all_standards_zv9)
tmp <- merge(nanopore[,c(1,5)], b.dav[,c(1,4)], by.x="transcript_id", by.y="X.transcript_id")
# tmp <- merge(nanopore[,c(1,7)], b.dav[,c(1,4)], by.x="transcript_id", by.y="X.transcript_id")
colnames(tmp) <- c("transcript_id", "n", "b")
ggplot(tmp, aes(x=n, y=b)) + geom_line(data = data.frame(x=c(0:150)), mapping = aes(x=x, y=x), color="grey") +
  geom_point(color="steelblue") + theme_minimal() +
  labs(title="PolyA length (nanopore vs dist_all_values_all_standards_zv9)",
       x="Nanopore polyA length", y="dist_all_values_all_standards_zv9 polyA length")

# ggplot(tmp, aes(x=n, y=b)) + geom_point(color="steelblue") + theme_minimal() +
#   labs(title="PolyA length (nanopore vs dist_all_values_all_standards_zv9)",
#        x="Nanopore polyA length in dwelltime", y="dist_all_values_all_standards_zv9 polyA length")

n_vs_b.dav_cor <- cor(tmp$n, tmp$b)
n_vs_b.dav_N <- length(tmp$n)

tmp$n <- round(tmp$n); tmp$b <- round(tmp$b)
tmp <- rbind(data.frame(l=tmp$n, class="nanopore"), data.frame(l=tmp$b, class="dav"))
ggplot(tmp, aes(x=l, fill=class)) + geom_histogram(aes(y=..density..), alpha=0.5, 
    position="identity", binwidth = 1) + geom_density(alpha=0.4) + 
  labs(title="PolyA length density (nanopore vs dist_all_values_all_standards_zv9)",
       x="polyA length in bp")
rm(tmp)

# Nanopore vs bartel (dist_from_publication)
tmp <- merge(nanopore[,c(1,5)], b.dfp[,c(1,4)], by.x="transcript_id", by.y="Transcript.ID")
#tmp <- merge(nanopore[,c(1,7)], b.dfp[,c(1,4)], by.x="transcript_id", by.y="Transcript.ID")
colnames(tmp) <- c("transcript_id", "n", "b")
tmp$b[is.na(tmp$b)] <- 0
ggplot(tmp, aes(x=n, y=b)) + geom_line(data = data.frame(x=c(0:150)), mapping = aes(x=x, y=x), color="grey") +
  geom_point(color="steelblue") + theme_minimal() +
  labs(title="PolyA length (nanopore vs dist_from_publication)",
       x="Nanopore polyA length", y="dist_from_publication polyA length")

# ggplot(tmp, aes(x=n, y=b)) + geom_point(color="steelblue") + theme_minimal() +
#   labs(title="PolyA length (nanopore vs dist_from_publication)",
#        x="Nanopore polyA length in dwelltime", y="dist_from_publication polyA length")

n_vs_b.dfp_cor <- cor(tmp$n, tmp$b)
n_vs_b.dfp_N <- length(tmp$n)

tmp$n <- round(tmp$n); tmp$b <- round(tmp$b)
tmp <- rbind(data.frame(l=tmp$n, class="nanopore"), data.frame(l=tmp$b, class="dav"))
ggplot(tmp, aes(x=l, fill=class)) + geom_histogram(aes(y=..density..), alpha=0.5, 
    position="identity", binwidth = 1) + geom_density(alpha=0.4) + 
  labs(title="PolyA length density (nanopore vs dist_from_publication)",
       x="polyA length in bp")
rm(tmp)

# Dwelltime per one basepair vs. fastq
ggplot(polyA, aes(x=Fastq_Length, y=DwellTimePerBP)) + geom_point(color="orange") + 
  theme_minimal() + labs(title="Fastq length vs. Dwelltime Per One Basepair")
