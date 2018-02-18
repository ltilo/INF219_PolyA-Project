#!/bin/bash

###############################################################################
# General uncompressing and looking what I have
###############################################################################

# Uncompressing the data into a single *.tar (NB! ca. 30 min on my mashine)
# I noticed quiet slow access to the data before uncompressing.
gunzip shield.tar.gz

# Show tar structure
tar --list -f  shield.tar | grep '/$'

# Extract pass/0 (first 4000 reads)
tar xf shield.tar /pass/0

# Take out the 0 directory in pass and delete pass
mv pass/* ./
rmdir pass

# Show first 10 files of 0 diroctory
cd 0
ls -d $(ls | head -10)
cd ../

# Looking at the poretools stats on all data (NB! ca. 40 min on my mashine)
poretools stats shield.tar

###############################################################################
# Running Poretools, MiniMap2 aligner & Samtools on the first 4000 reads
###############################################################################

# Running poretools stats and obtain fastq for 0 directory
poretools stats 0/
poretools fastq 0/ > shield_pass_0.fastq

# Show head of the fastq file and check if the file seems to be fine
head shield_pass_0.fastq

# Using minimap2 aligner
~/Application/minimap2/minimap2 -ax splice -l14 -t 8 -uf /Bartel/Danio_rerio.Zv9.fa shield_pass_0.fastq > shield_pass_0.sam 

# Using samtools to convert sam file to bam file
samtools view -Sb shield_pass_0.sam > shield_pass_0.bam

###############################################################################
# Running Poretools, MiniMap2 aligner & Samtools on all data
# (NB: This may take some time)
# (Not done yet)
###############################################################################

# Running poretools to obtain fastq file
poretools fastq shield.tar > shield_pass.fastq

# Looking into the fastq file, just for fun
head shield_pass.fastq

# Using minimap2 aligner
~/Application/minimap2/minimap2 -ax splice -l14 -t 8 -uf /Bartel/Danio_rerio.Zv9.fa shield_pass.fastq > shield_pass.sam

# Using samtools to convert sam file to bam file
samtools view -Sb shield_pass.sam > shield_pass.bam
