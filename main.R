#'########################################################################################
# Set path to files                                                                   ####
#'########################################################################################

# Set all path to nanopore shield data
path_bam <- "../Data/shield_pass.bam"
path_fastq <- "../Data/shield_pass.fastq"

path_data <- "../Data/"

# Set all path to Bartel data (Shield - 6hpf)
path_gtf <- "../Data/Bartel/Danio_rerio.Zv9.79.gtf"
path_b.dmv <- "../Data/Bartel/PAL_seq/dist_median_values_subset_of_standards_zv9/Dre_mock_6hpf_PAL_length_dist_only_10_50_100_lengths.csv"
path_b.dav <- "../Data/Bartel/PAL_seq/dist_all_values_all_standards_zv9/Dre_mock_6hpf_PAL_length_dist_all_lengths.csv"
path_b.dfp <- "../Data/Bartel/PAL_seq/dist_from_publication/GSE52809_Dre_mock_6hpf.txt"

# Source necessary files
source('read_plot_functions.R')
Rcpp::sourceCpp('findPolyA.cpp')


#'########################################################################################
# Read *.bam and *.gtf file -- Find overlaps                                          ####
#'########################################################################################

# Read *bam and *gtf file
bam <- readGAlignments(path_bam, use.names = T)
gtf <- makeTxDbFromGFF(path_gtf, format="gtf")
  
# Making reference transcripts and filter only coding genes
transcripts <- exonsBy(gtf, by = "tx", use.names=T)
transcripts <- transcripts[names(cdsBy(gtf, by = "tx", use.names = T))]
  
# Finding overlaps
hits <- findOverlaps(transcripts, bam)

# Analyze fastQ
fast5_register = analyzeFastQ(normalizePath(path_fastq))
fast5_register$fast5_path <- paste0(path_data, fast5_register$fast5_path)
fast5_register <- data.frame(fast5_register, stringsAsFactors = F)

# Compute .*fast5 traceback list
subjectHits <- unique(to(hits))
f5list <- data.frame(subjectHits, names(bam[subjectHits]), stringsAsFactors = F)
rm(subjectHits)
colnames(f5list) <- c("subjectHits", "traceName")
f5list <- merge(fast5_register, f5list, by = "traceName")


#'########################################################################################
# Read *.fast5 files -- Find PolyA tail                                               ####
#'########################################################################################

# Read .*fast5 files and find polyA for every sample
n.shield.data <- read_and_findPolyA(f5list)

# Compute additional StartTime_DevidedBy_SampleFreq_groupedByHour
n.shield.data$StartTime_DevidedBy_SampleFreq_groupedByHour <- n.shield.data$Attribute_StartTime/n.shield.data$Attribute_SampleFrequency
n.shield.data$StartTime_DevidedBy_SampleFreq_groupedByHour <- n.shield.data$StartTime_DevidedBy_SampleFreq_groupedByHour / 3600

# Save data
write.csv(n.shield.data, "nanopore_shield_allOverlappingToBartel_polyA.csv")

# Round and factor data
n.shield.data$StartTime_DevidedBy_SampleFreq_groupedByHour <- round(n.shield.data$StartTime_DevidedBy_SampleFreq_groupedByHour)
n.shield.data$StartTime_DevidedBy_SampleFreq_groupedByHour <- as.factor(n.shield.data$StartTime_DevidedBy_SampleFreq_groupedByHour)
n.shield.data$Attribute_StartMux <- as.factor(n.shield.data$Attribute_StartMux)
n.shield.data$Attribute_ChannelNumber <- as.factor(n.shield.data$Attribute_ChannelNumber)

# Filter data
n.shield.dataOriginal <- n.shield.data
n.shield.data <- n.shield.data[n.shield.data$startPolyA_dwelltime > 0, ]
n.shield.data <- n.shield.data[n.shield.data$DwelltimePerBasepair < 255, ]
n.shield.data <- n.shield.data[n.shield.data$lengthPolyA_BasePair < 500, ]

# Compute final polyA table for nanopore-shield transcripts
n.shield.transcriptTable <- makeTranscriptIdTable(hits, n.shield.data)
n.shield.transcriptTable <- n.shield.transcriptTable[n.shield.transcriptTable$tags > 0,]

# Save data
write.csv(n.shield.transcriptTable, "nanopore_shield_transcriptId_table.csv")


#'########################################################################################
# Read Bartel data                                                                    ####
#'########################################################################################

# Read in Bartel data
b.dmv <- read.csv(path_b.dmv, header = T, sep = '\t', stringsAsFactors = F)
b.dav <- read.csv(path_b.dav, header = T, sep = '\t', stringsAsFactors = F)
b.dfp <- read.csv(path_b.dfp, header = T, sep = '\t', stringsAsFactors = F)


#'########################################################################################
# Plotting                                                                            ####
#'########################################################################################

# Plot polyA length histogram
plotHistogram(n.shield.data, 7, c("Nanopore Shield PolyA Length", "PolyA Length In Basepair"))

# Plot Dwelltime Per BasePair histogram
plotHistogram(n.shield.data, 8, c("Nanopore Shield Dwelltime Per Basepair", "Dwelltime per Basepair"))


### Plots against Bartel data - dist_median_values_subset_of_standards_zv9 (b.dmv)
plotScatter(
  n.shield.transcriptTable[,c(1,3)], 
  b.dmv[,c(1,4)], 
  c("transcript_id", "X.transcript_id"),
  c("Nanopore Shield Vs. dist_median_values_subset_of_standards_zv9 -- PolyA Length",
    "Nanopore Shield PolyA Length", "dist_median_values_subset_of_standards_zv9 PolyA Length"
  )
)
plotDensity(
  n.shield.transcriptTable[,c(1,3)], 
  b.dmv[,c(1,4)], 
  c("transcript_id", "X.transcript_id"),
  c("Nanopore Shield Vs. dist_median_values_subset_of_standards_zv9 -- Density PolyA Length",
    "PolyA Length In Basepair"
  ),
  c("Nanopore Shield", "dist_median_values_subset_of_standards_zv9")
)


### Plots against Bartel data - dist_all_values_all_standards_zv9 (b.dav)
plotScatter(
  n.shield.transcriptTable[,c(1,3)], 
  b.dav[,c(1,4)], 
  c("transcript_id", "X.transcript_id"),
  c("Nanopore Shield Vs. dist_all_values_all_standards_zv9 -- PolyA Length",
    "Nanopore Shield PolyA Length", "dist_all_values_all_standards_zv9 PolyA Length"
  )
)
plotDensity(
  n.shield.transcriptTable[,c(1,3)], 
  b.dav[,c(1,4)], 
  c("transcript_id", "X.transcript_id"),
  c("Nanopore Shield Vs. dist_all_values_all_standards_zv9 -- Density PolyA Length",
    "PolyA Length In Basepair"
  ),
  c("Nanopore Shield", "dist_all_values_all_standards_zv9")
)


### Plots against Bartel data - dist_from_publication (b.dfp)
plotScatter(
  n.shield.transcriptTable[,c(1,3)], 
  b.dfp[,c(1,4)], 
  c("transcript_id", "Transcript.ID"),
  c("Nanopore Shield Vs. dist_from_publication -- PolyA Length",
    "Nanopore Shield PolyA Length", "dist_from_publication PolyA Length"
  )
)
plotDensity(
  n.shield.transcriptTable[,c(1,3)], 
  b.dfp[,c(1,4)], 
  c("transcript_id", "Transcript.ID"),
  c("Nanopore Shield Vs. dist_from_publication -- Density PolyA Length",
    "PolyA Length In Basepair"
  ),
  c("Nanopore Shield", "dist_from_publication")
)

# Plot Dwelltime Vs. DwellTimePerBp (x = Max500000, y = Max250)
plot <- ggplot(n.shield.data, aes(x=DwellTime, y=DwelltimePerBasepair)) + geom_point(alpha = 0.005) +
  theme_minimal() + labs(title="Nanopore Shield - Dwelltime for whole sample Vs. DwelltimePerBp") +
  coord_cartesian(xlim=c(0, 500000), ylim=c(0, 250))
plot

# Boxplot for StartTime
plot <- ggplot(n.shield.dataOriginal, aes(x=StartTime_DevidedBy_SampleFreq_groupedByHour, 
          y=DwelltimePerBasepair, color=StartTime_DevidedBy_SampleFreq_groupedByHour)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore Shield - Boxplot StartTime/SampleFreq grouped by hours Vs. DwellTimePerBp")
plot
plot <- ggplot(n.shield.data, aes(x=StartTime_DevidedBy_SampleFreq_groupedByHour, 
                                          y=DwelltimePerBasepair, color=StartTime_DevidedBy_SampleFreq_groupedByHour)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore Shield - Boxplot StartTime/SampleFreq grouped by hours Vs. DwellTimePerBp")
plot
rm(plot)

# Boxplot for mux groups
plot <- ggplot(n.shield.dataOriginal, aes(x=Attribute_StartMux, 
          y=DwelltimePerBasepair, color=Attribute_StartMux)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore Shield - Boxplot Attribute_StartMux Vs. DwellTimePerBp")
plot
plot <- ggplot(n.shield.data, aes(x=Attribute_StartMux, 
          y=DwelltimePerBasepair, color=Attribute_StartMux)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore Shield - Boxplot Attribute_StartMux Vs. DwellTimePerBp")
plot
rm(plot)

# Boxplot for ChannelNumber
plot <- ggplot(n.shield.dataOriginal, aes(x=Attribute_ChannelNumber, 
          y=DwelltimePerBasepair)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore Shield - Boxplot Attribute_ChannelNumber Vs. DwellTimePerBp")
plot
plot <- ggplot(n.shield.data, aes(x=Attribute_ChannelNumber, 
          y=DwelltimePerBasepair)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore Shield - Boxplot Attribute_ChannelNumber Vs. DwellTimePerBp")
plot
rm(plot)

#'########################################################################################
# Compute Covariance                                                                  ####
#'########################################################################################

cov1 <- computeCorrelation(
  n.shield.transcriptTable[,c(1,3)], 
  b.dmv[,c(1,4)], 
  c("transcript_id", "X.transcript_id")
)

cov2 <- computeCorrelation(
  n.shield.transcriptTable[,c(1,3)], 
  b.dav[,c(1,4)], 
  c("transcript_id", "X.transcript_id")
)

cov3 <- computeCorrelation(
  n.shield.transcriptTable[,c(1,3)], 
  b.dfp[,c(1,4)], 
  c("transcript_id", "Transcript.ID")
)


