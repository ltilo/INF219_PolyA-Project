#!/bin/bash

# Unpack the data into a single *.tar (NB! ca. 30 min on my mashine)
gunzip shield.tar.gz

# Run poretools stats on the data (NB! ca. 40 min on my mashine)
poretools stats shield.tar

# Show tar structure
tar --list -f  shield.tar | grep '/$'

# Extract pass/0 (first 4000 reads)
tar xf shield.tar /pass/0

# Take out the 0 directory in pass and delete pass
mv pass/* ./
rmdir pass

# Running poretools stats and fastq for 0 directory
poretools stats 0/
poretools fastq 0/ > shield_pass_0.fastq

# Show head of the fastq file and check if the file seems to be fine
head shield_pass_0.fastq

# Using minimap2 aligner
~/Application/minimap2/minimap2 -ax splice -l14 -t 8 -uf /Bartel/Danio_rerio.Zv9.fa shield_pass_0.fastq > shield_pass_0.sam 

# Using samtools to convert sam file to bam file
samtools view -Sb shield_pass_0.sam > shield_pass_0.bam
