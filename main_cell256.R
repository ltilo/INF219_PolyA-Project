#'########################################################################################
# Set path to files                                                                   ####
#'########################################################################################

# Set all path to nanopore 256cell data
path_bam <- "../Data/256cell_pass.bam"
path_fastq <- "../Data/256cell_pass.fastq"

path_data <- "../Data/"

# Set all path to Bartel data (2hpf)
path_gtf <- "../Data/Bartel/Danio_rerio.Zv9.79.gtf"
path_b.dmv <- "../Data/Bartel/PAL_seq/dist_median_values_subset_of_standards_zv9/Dre_mock_2hpf_PAL_length_dist_only_10_50_100_lengths.csv"
path_b.dav <- "../Data/Bartel/PAL_seq/dist_all_values_all_standards_zv9/Dre_mock_2hpf_PAL_length_dist_all_lengths.csv"
path_b.dfp <- "../Data/Bartel/PAL_seq/dist_from_publication/GSE52809_Dre_mock_2hpf.txt"

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
n.256cell.data <- read_and_findPolyA(f5list)

# Compute additional StartTime_DevidedBy_SampleFreq_groupedByHour
n.256cell.data$StartTime_DevidedBy_SampleFreq_groupedByHour <- n.256cell.data$Attribute_StartTime/n.256cell.data$Attribute_SampleFrequency
n.256cell.data$StartTime_DevidedBy_SampleFreq_groupedByHour <- n.256cell.data$StartTime_DevidedBy_SampleFreq_groupedByHour / 3600

# Save data
write.csv(n.256cell.data, "nanopore_256cell_allOverlappingToBartel_polyA.csv")

# Round and factor data
n.256cell.data$StartTime_DevidedBy_SampleFreq_groupedByHour <- round(n.256cell.data$StartTime_DevidedBy_SampleFreq_groupedByHour)
n.256cell.data$StartTime_DevidedBy_SampleFreq_groupedByHour <- as.factor(n.256cell.data$StartTime_DevidedBy_SampleFreq_groupedByHour)
n.256cell.data$Attribute_StartMux <- as.factor(n.256cell.data$Attribute_StartMux)
n.256cell.data$Attribute_ChannelNumber <- as.factor(n.256cell.data$Attribute_ChannelNumber)

# Filter data
n.256cell.dataOriginal <- n.256cell.data
n.256cell.data <- n.256cell.data[n.256cell.data$startPolyA_dwelltime > 0, ]
n.256cell.data <- n.256cell.data[n.256cell.data$DwelltimePerBasepair < 255, ]
n.256cell.data <- n.256cell.data[n.256cell.data$lengthPolyA_BasePair < 500, ]

# Compute final polyA table for nanopore-256cell transcripts
n.256cell.transcriptTable <- makeTranscriptIdTable(hits, n.256cell.data)
n.256cell.transcriptTable <- n.256cell.transcriptTable[n.256cell.transcriptTable$tags > 0,]

# Save data
write.csv(n.256cell.transcriptTable, "nanopore_256cell_transcriptId_table.csv")


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
plotHistogram(n.256cell.data, 7, c("Nanopore 256cell PolyA Length", "PolyA Length In Basepair"))

# Plot Dwelltime Per BasePair histogram
plotHistogram(n.256cell.data, 8, c("Nanopore 256cell Dwelltime Per Basepair", "Dwelltime per Basepair"))


### Plots against Bartel data - dist_median_values_subset_of_standards_zv9 (b.dmv)
plotScatter(
  n.256cell.transcriptTable[,c(1,3)], 
  b.dmv[,c(1,4)], 
  c("transcript_id", "X.transcript_id"),
  c("Nanopore 256cell Vs. dist_median_values_subset_of_standards_zv9 -- PolyA Length",
    "Nanopore 256cell PolyA Length", "dist_median_values_subset_of_standards_zv9 PolyA Length"
  )
)
plotDensity(
  n.256cell.transcriptTable[,c(1,3)], 
  b.dmv[,c(1,4)], 
  c("transcript_id", "X.transcript_id"),
  c("Nanopore 256cell Vs. dist_median_values_subset_of_standards_zv9 -- Density PolyA Length",
    "PolyA Length In Basepair"
  ),
  c("Nanopore 256cell", "dist_median_values_subset_of_standards_zv9")
)


### Plots against Bartel data - dist_all_values_all_standards_zv9 (b.dav)
plotScatter(
  n.256cell.transcriptTable[,c(1,3)], 
  b.dav[,c(1,4)], 
  c("transcript_id", "X.transcript_id"),
  c("Nanopore 256cell Vs. dist_all_values_all_standards_zv9 -- PolyA Length",
    "Nanopore 256cell PolyA Length", "dist_all_values_all_standards_zv9 PolyA Length"
  )
)
plotDensity(
  n.256cell.transcriptTable[,c(1,3)], 
  b.dav[,c(1,4)], 
  c("transcript_id", "X.transcript_id"),
  c("Nanopore 256cell Vs. dist_all_values_all_standards_zv9 -- Density PolyA Length",
    "PolyA Length In Basepair"
  ),
  c("Nanopore 256cell", "dist_all_values_all_standards_zv9")
)


### Plots against Bartel data - dist_from_publication (b.dfp)
plotScatter(
  n.256cell.transcriptTable[,c(1,3)], 
  b.dfp[,c(1,4)], 
  c("transcript_id", "Transcript.ID"),
  c("Nanopore 256cell Vs. dist_from_publication -- PolyA Length",
    "Nanopore 256cell PolyA Length", "dist_from_publication PolyA Length"
  )
)
plotDensity(
  n.256cell.transcriptTable[,c(1,3)], 
  b.dfp[,c(1,4)], 
  c("transcript_id", "Transcript.ID"),
  c("Nanopore 256cell Vs. dist_from_publication -- Density PolyA Length",
    "PolyA Length In Basepair"
  ),
  c("Nanopore 256cell", "dist_from_publication")
)

# Plot Dwelltime Vs. DwellTimePerBp (x = Max500000, y = Max250)
plot <- ggplot(n.256cell.data, aes(x=DwellTime, y=DwelltimePerBasepair)) + geom_point(alpha = 0.1) +
  theme_minimal() + labs(title="Nanopore 256cell - Dwelltime for whole sample Vs. DwelltimePerBp") +
  coord_cartesian(xlim=c(0, 500000), ylim=c(0, 250))
ggsave(paste0(plot$labels$title, ".png"), plot=plot, path="Result/", width=17.3, height=7.06, dpi=125)

# Boxplot for StartTime
plot <- ggplot(n.256cell.dataOriginal, aes(x=StartTime_DevidedBy_SampleFreq_groupedByHour, 
        y=DwelltimePerBasepair, color=StartTime_DevidedBy_SampleFreq_groupedByHour)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore 256cell - Boxplot StartTime/SampleFreq grouped by hours Vs. DwellTimePerBp")
ggsave(paste0("Nanopore 256cell_Boxplot_StartTimeVsDwellTimePerBp", ".png"), plot=plot, path="Result/", width=17.3, height=7.06, dpi=125)

plot <- ggplot(n.256cell.data, aes(x=StartTime_DevidedBy_SampleFreq_groupedByHour, 
       y=DwelltimePerBasepair, color=StartTime_DevidedBy_SampleFreq_groupedByHour)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore 256cell - Boxplot StartTime/SampleFreq grouped by hours Vs. DwellTimePerBp")
ggsave(paste0("Nanopore 256cell_Boxplot_StartTimeVsDwellTimePerBp_zoomed", ".png"), plot=plot, path="Result/", width=17.3, height=7.06, dpi=125)


# Boxplot for mux groups
plot <- ggplot(n.256cell.dataOriginal, aes(x=Attribute_StartMux, 
         y=DwelltimePerBasepair, color=Attribute_StartMux)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore 256cell - Boxplot Attribute_StartMux Vs. DwellTimePerBp")
ggsave(paste0(plot$labels$title, ".png"), plot=plot, path="Result/", width=17.3, height=7.06, dpi=125)

plot <- ggplot(n.256cell.data, aes(x=Attribute_StartMux, 
         y=DwelltimePerBasepair, color=Attribute_StartMux)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore 256cell - Boxplot Attribute_StartMux Vs. DwellTimePerBp")
ggsave(paste0(plot$labels$title, "_zoomed", ".png"), plot=plot, path="Result/", width=17.3, height=7.06, dpi=125)


# Boxplot for ChannelNumber
plot <- ggplot(n.256cell.dataOriginal, aes(x=Attribute_ChannelNumber, 
             y=DwelltimePerBasepair)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore 256cell - Boxplot Attribute_ChannelNumber Vs. DwellTimePerBp")
ggsave(paste0(plot$labels$title, ".png"), plot=plot, path="Result/", width=17.3, height=7.06, dpi=125)

plot <- ggplot(n.256cell.data, aes(x=Attribute_ChannelNumber, 
           y=DwelltimePerBasepair)) + 
  geom_boxplot() + theme_minimal() + 
  labs(title="Nanopore 256cell - Boxplot Attribute_ChannelNumber Vs. DwellTimePerBp")
ggsave(paste0(plot$labels$title, "_zoomed", ".png"), plot=plot, path="Result/", width=17.3, height=7.06, dpi=125)
rm(plot)

#'########################################################################################
# Compute Covariance                                                                  ####
#'########################################################################################

cov1_256 <- computeCorrelation(
  n.256cell.transcriptTable[,c(1,3)], 
  b.dmv[,c(1,4)], 
  c("transcript_id", "X.transcript_id")
)

cov2_256 <- computeCorrelation(
  n.256cell.transcriptTable[,c(1,3)], 
  b.dav[,c(1,4)], 
  c("transcript_id", "X.transcript_id")
)

cov3_256 <- computeCorrelation(
  n.256cell.transcriptTable[,c(1,3)], 
  b.dfp[,c(1,4)], 
  c("transcript_id", "Transcript.ID")
)


