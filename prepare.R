#################################################################################################
##  library(rhdf5)
#################################################################################################

# Finding path and path to files
getwd(); setwd(); # Get/Set path
path <- "/home/florian/Schreibtisch/Uni/Inf219"
files <- list.files(path = path, pattern = ".fast5$", full.names = TRUE)

# Show information about the first hdf5 files
file <- files[1]
h5ls(file)

# Grab some raw data
tmpPath <- h5ls(file)[19,1]
data.raw <- h5read(file, tmpPath)
data.raw <- data.raw$Signal

# Grab some BaseCalled_template data
tmpPath <- h5ls(file)[4,1]
data.events <- h5read(file, tmpPath)$Events

# Plotting the raw data with ggplot2
library(ggplot2)
time <- c(1:length(data.raw))
time <- time/4000
df <- data.frame(data.raw, time)
colnames(df) <- c("Raw", "Time")

ggplot(df, aes(Time, Raw)) + geom_line() #plotting data as line
ggplot(df, aes(Time, Raw)) + geom_point() #plotting data as points
ggplot(df, aes(Time, Raw)) + geom_line() + geom_point() #Plotting both lines & points

#################################################################################################
##  library(IONiseR)
#################################################################################################

# Show some information about data
mData <- readFast5Summary(files)
print(mData)
readInfo(mData)
eventData(mData)
baseCalled(mData)
fastq(m)

# Get FastQ as ShortRead object
library(ShortRead)
seq <- sread(fastq(mData))