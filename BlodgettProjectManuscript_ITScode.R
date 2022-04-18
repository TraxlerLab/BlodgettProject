library(data.table)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(vegan)
library(phyloseq)
#
####
### AMPtk pre-processing --> ITS DADA2 workflow: 
### https://github.com/nextgenusfs/amptk/blob/master/amptk/dada2_pipeline_nofilt.R 
### https://benjjneb.github.io/dada2/ITS_workflow.html 
####
#
###############################################################
#### Package installation and Overview Summary of Pipeline ####
###############################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ShortRead")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

## AMPtk installation:
# https://amptk.readthedocs.io/en/latest/index.html


### OVERVIEW
# Step 1: run the AMPtk pre-processing step in terminal on local machine
# Step 2: explore the AMPtk preprocessing output in R, check quality, etc. on local machine
# Step 3: run DADA2 (test run on local machine with 1-2 samples, then as a SBATCH on Savio)
# Step 4: run decontam from R on local machine to account for taxa present in negative controls
# Step 5: stats!

##############################################################
#### STEP 1: run the AMPtk pre-processing step in terminal ####
###############################################################
#
# $ cd ~/BlodgettAmpliconSequencingRawData/
# $ conda activate amptk_v1.5.1
# $ amptk illumina -i allITS -o allITS_AMPTKpreprocessed -f AACTTTYRRCAAYGGATCWCT -r AGCCTCCGCTTATTGATATGCTTAART --rescue_forward off --primer_mismatch 6
#### Runtime on my MacbookPro + 4TB external hard drive = 2hrs 20min
#### total number of samples = 664

##################################################################################
#### STEP 2: explore the AMPtk preprocessing output in R, check quality, etc. ####
##################################################################################
library(dada2)
library(ShortRead)
library(Biostrings)

# Check for any primers remaining in sequences after AMPtk processing (from the DADA2 ITS tutorial)
# Create objects that are your primer sequences
FWD <- "AACTTTYRRCAAYGGATCWCT"  #5.8sfun
REV <- "AGCCTCCGCTTATTGATATGCTTAART"  #ITS4fun

# Verify the correct sequence and orientation of those primer sequences!
# Create all orientations of the input sequence
allOrients <- function(primer) { 
  require(Biostrings)
  dna <- DNAString(primer)  #Biostrings package works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients # "Forward" sequence matches what I assigned above as FWD
REV.orients # "Forward" sequence matched what I assigned above as REV

# Count the number of times the primers appear in the forward and reverse read, 
# while considering all possible primer orientations. 
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Designate File Path to the AMPtk output folder
path <- "~/BlodgettAmpliconSequencingRawData/allITS_AMPTKpreprocessed"
# check that the path is correct
list.files(path)
#each sample has five different file types:
# (1) _R1.fq      ## forward reads, left primer removed, low quality reads removed
# (2) _R2.fq      ## reverse reads, left primer removed, low quality reads removed
# (3) .demux.fq   ## merged reads + additional quality filtering/trimming!
# (4) .merged.fq  ## merged reads
# (5) .stats      ## statistics...
## This .stats file coontains a string of 6 different numbers:
## (1) Initial total number of input reads 
## (2) number of reads where the F primer was found                          
## (3) number of reads where the R primer was found         
## (4) number of reads where the multiple primers were found (discarded due to no single primer being found)                                   
## (5) ?number of reads discarded for being too short (<100bp)        
## (6) total number of reads output at end?     
#
# Jon Palmer says he uses the DADA2 pipeline as if working with single reads instead of paired reads
# because quality filtering and trimming is part of the merging process in AMPtk
#
# make a table out of the .stats files:
path <- "~/BlodgettAmpliconSequencingRawData/allITS_AMPTKpreprocessed"
filelist <- list.files(path, pattern=".stats") #list the filenames = sampleIDs
#create a list of all the file paths to the .stats files:
fullpathfilelist <- list.files(path, pattern=".stats", full.names=TRUE) 
#apply the fread function to fullpathfilelist list of file paths:
allfiles <- lapply(fullpathfilelist, fread)
length(allfiles) #check that it's the correct number of files
#turn the allfiles object into a data.table with rbindlist()!
AMPtkStats <- rbindlist(allfiles)
#add the filesnames/sampleIDs onto the list as a new column:
AMPtkStats[ ,FileNames:=filelist]
#remove the .stat from the end of each filename:
AMPtkStats[ ,c("SampleID", "trash"):=tstrsplit(FileNames, ".st", fixed=FALSE)]
#remove excess columns:
AMPtkStats[ ,c("trash", "FileNames"):=NULL] 
ColnamesList <- c("Input", "ReadsWithFprimer", "ReadsWithRprimer", "ReadsWithMultiplePrimers", "TooShortReads", "Output", "SampleID")
names(AMPtkStats) <- ColnamesList
write.csv(AMPtkStats, "/Users/monikafischer/Desktop/AMPtkStats.csv")

# Create a list of the file names for forward and reverse files:
#raw:
path <- "/Users/monikafischer/Desktop/AmpliconSeqStuff/supersmall/supersmallsubset_rawdata"
fnRawFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRawRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))
plotQualityProfile(fnRawFs[2])
plotQualityProfile(fnRawRs[2])

#AMPtk output #1, low quality reads removed
path <- "/Users/monikafischer/Desktop/AmpliconSeqStuff/supersmall/supersmallsubset_preprocessed/"
fnFs <- sort(list.files(path, pattern = "_R1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fq", full.names = TRUE))
plotQualityProfile(fnFs[2]) #plot the quality profile for a single fastq file
plotQualityProfile(fnRs[2]) #plot the quality profile for a single fastq file
plotQualityProfile(fnFs, aggregate = TRUE) #aggreate quality profile data from all fastq files into one plot
plotQualityProfile(fnRs[1])
# grey-scale heatmap = distribution of quality scores at each position
# green = mean
# orange solid = median
# orange dashed = 25th and 75th quantiles
# red line (only present if seqs vary in length) = % of seqs at each length

##AMPtk output #2, merged reads:
fnMerged <- sort(list.files(path, pattern = ".merged.fq", full.names = TRUE))
plotQualityProfile(fnMerged[2])

#
##AMPtk output #3, demuxed merged reads:
path <- "~/BlodgettAmpliconSequencingRawData/allITS_AMPTKpreprocessed"
FNs <- sort(list.files(path, pattern = "demux.fq", full.names = TRUE))
plotQualityProfile(FNs[1])
plotQualityProfile(FNs, aggregate = TRUE)

# Run the primerHits function on one of the filteredN files
# change the number in the brackets to change which file you're looking at!
# number in brackets associated with the row number:
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = FNs[2]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = FNs[2]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = FNs[2]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = FNs[2]))

# file types:
FNs #demuxed... AMPtk output #3
fnRawFs #Should contain primers... Raw Forward Reads
fnRawRs #Should contain primers... Raw Reverse Reads
fnFs #left primer should be removed... AMPtk output #1
fnRs #left primer should be removed... AMPtk output #1
fnMerged #merged reads... AMPtk output #2

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnRawFs[2]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRawRs[2]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnRawFs[2]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRawRs[2]))





######################################################################################################
#### STEP 3: RUN DADA2! (practice on a tiny subset of data here, but for reals, do this on Savio) ####
######################################################################################################
#
## test run on local machine with two samples, 
## then convert this into a script to run on the Savio supercomputer!
#
## Following Jon Palmer's DADA2 code, which also follows the ITS tutorial:
library(dada2)
# load the data from a folder
# we only care about the .demux.fq file, which is fully trimmed, quality controlled, and paired reads are merged
path <- "~/BlodgettAmpliconSequencingRawData/TESTdata_AMPTKpreprocessed"
fns <- list.files(path, pattern = "\\.demux.fq$")
## "list the files within the path that contain the pattern"
## The Pattern (a regular expression, a.k.a. regex):
# \\. --matches a literal .
# demux.fq --matches a literal sequence of characters (in this case, "demux.fq")
# $ --end of string
Seqs <- file.path(path, fns)
Seqs

#get sample names
sample.names <- list.files(path, pattern = "\\.demux.fq$", full.name = FALSE)
sample.names <- gsub('[.demux.fq]', '', sample.names)
sample.names

#Dereplication (remove duplicated sequences)
derepSeqs <- derepFastq(Seqs, verbose=TRUE)

#name the derep class with sample names
names(derepSeqs) <- sample.names

#Sample Inference = infer ASVs!
dadaSeqs <- dada(derepSeqs, err=NULL, selfConsist=TRUE, BAND_SIZE=32)
## "err" is normally the output from learnErrors(), which I skipped over by using AMPtk for pre-processing
## "selfConsist=TRUE" is the algorithm will alternate between sample inference and error rate estimation until convergence.
## "BAND_SIZE=32" When set, banded Needleman-Wunsch alignments are performed. Banding restricts the net cumulative number 
## of insertion of one sequence relative to the other. The default value of BAND_SIZE is 16
# multithread=FALSE (default)
# USE_QUALS=TRUE: (default) If TRUE, the dada2 error model takes into account the consensus quality score of the 
# dereplicated unique sequences. If FALSE, quality scores are ignored.

#make sequence table
seqtab <- makeSequenceTable(dadaSeqs, orderBy = "abundance")
dim(seqtab)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
dim(seqtab.nochim)

# Export ASV table:
write.csv(seqtab.nochim, "~/seqtab.nochim.csv")

## Sanity check -- make sure the chimera removal step only removed a few reads from each sample
## if chimera removal resulted in a substantial loss in reads, this is could indicate
## that primers were removed inappropriately.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
track <- cbind(sapply(dadaSeqs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("dadaSeqs", "nonchim")
rownames(track) <- sample.names
track

# Assign Taxonomy!
unite.ref <- "~/UNITE_sh_general_release_dynamic_04.02.2020.fasta"
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, verbose = TRUE)

# EXPORT TAXONOMY TABLE
write.csv(taxa, "~/ALL_ITS_taxa.csv")

#
##
#
## Post-DADA2 Sanity check Histograms
#
## sanity check to make sure the distribution of the data make sense:
library(data.table)
library(ggplot2)
otu <- fread(file.choose()) #load OTUtable
otu[1:10, 1:10] #rownames/SeqIDs are now column1!

# count the number of values in each column that are greater than zero
NumberOfSamplesWithEachOTU <- data.frame(columnsums = colSums(otu[,-1] > 0))
# [,-1] skips the first column, which is the SeqIDs...

#plot histogram:
ggplot(NumberOfSamplesWithEachOTU, aes(x=columnsums))+
  geom_histogram(bins=500)+
  theme_minimal()+
  scale_x_continuous(breaks=seq(0, 15, 1), limits=c(0, 15))+
  scale_y_continuous(breaks=seq(0, 140000, 10000), limits=c(0, 140000))+
  xlab("Binned number of Samples/OTU") +
  ylab("Count") +
  ggtitle("ITS Histogram of Samples/OTU")

# count the number of values in each row that are greater than zero
NumberOfOTUsInEachSample <- data.frame(rowsums = rowSums(otu > 0))
#plot histogram:
ggplot(NumberOfOTUsInEachSample, aes(x=rowsums))+
  geom_histogram(bins=500)+
  theme_minimal()+
  scale_x_continuous(breaks=seq(0, 1400, 100), limits=c(0, 1400))+
  scale_y_continuous(breaks=seq(0, 30, 5), limits=c(0, 30))+
  xlab("Binned number of OTUs/Sample") +
  ylab("Count") +
  ggtitle("ITS Histogram of OTUs/Sample")

# how many samples have fewer than 30 OTUs? ...play with this as a compliment to the histograms!
sum(NumberOfOTUsInEachSample < 30)

#########################################################
#### STEP 4: CREATE PHYLOSEQ OBJECT and RUN DECONTAM ####
#########################################################
##
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html#identifying-contaminants-in-marker-gene-and-metagenomics-data
##
library(phyloseq)
library(ggplot2)
library(decontam)
library(data.table)
library(tidyverse)
library(vegan)

## decontam removes ASVs that are over-represented in negative-controls
## decontam uses dada2 outputs converted into phyloseq objects, and phyloseq requires OTU and Taxa tables as matrices
#
# Load OTU table:
seqtab.nochim <- as.matrix(fread("~/20210118_allITS_seqtab.nochim.csv"),
                           rownames=1)
seqtab.nochim[1:5,1:2] 

# Load Taxonomy Table:
taxa <- as.matrix(fread("~/20210119_allITS_taxa.csv",
                        header = TRUE), rownames=1)
head(taxa)

# Extract the SeqIDs from the OTUtable
samplenames <- data.table(SeqID = row.names(seqtab.nochim))
dim(samplenames)
# add read counts per sample to the Sample Metadata Table.. (reads output by DADA2 in the post-chimera-removal sanity check)
reads <- fread("~/20210118_allITS_trackreads.csv") 
head(reads)
colnames(reads) <- c("SeqID", "dada2_reads", "nochim_reads")

class(reads$nochim_reads)
class(reads$dada2_reads)
reads$nochim_reads <- as.numeric(reads$nochim_reads)
max(reads$nochim_reads)
min(reads$nochim_reads)

# subset metadata table, add reads counts, and create the Sample Table for phyloseq
metasub1 <- merge(samplenames, SampleMetadata, by="SeqID", all.x=TRUE)
metasub2 <- merge(metasub1, reads, by="SeqID", all.x=TRUE)
sampITS <- column_to_rownames(metasub2, var="SeqID") #rownames = SeqID will be required for phyloseq

#plot the number of reads by sample to get an general idea of what the data look like
ggplot(data=metasub2, aes(x=SeqID, y=nochim_reads, color=Who)) +
  geom_point() + 
  scale_color_manual(values = c("#00ff95", "#ff9500", "#9500ff", "#ff006a")) +
  scale_y_continuous(expand= c(0,0), breaks=seq(0, 130000, 10000), limits=c(0, 130000)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))

##
# CREATE PHYLOSEQ OBJECT
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), #matrix (rownames = SeqIDs, colnames = ASVs, values = abundance of each ASV in each SeqID)
               tax_table(taxa), #matrix (rownames = ASVs, colnames = taxonomic levels, values = taxonomic classification for each ASV)
               sample_data(sampITS)) #data.frame, (rownames = SeqIDs, colnames & values = additional info for each SeqID)

# check that nothing catastrophic happened when creating the phyloseq object, table dimensions should match:
dim(tax_table(ps)) #phyloseq object
dim(taxa) #original
(otu_table(ps))[1:5, 1:5] #phyloseq object  
dim(seqtab.nochim) #original         
dim(sample_data(ps)) #phyloseq object
head(sampITS) #original
head(tax_table(ps))


## ASVs are currently identified by their DNA sequence, which is cumbersome
## Move this DNA sequence to the refseq slot of the phyloseq object
## and then give each ASV a simpler name like ASV1, ASV2, etc.
#
# move the current taxa names to a DNAStringSet object, which phyloseq will recognize at DNAseqs
dna <- Biostrings::DNAStringSet(taxa_names(ps)) 
# connect taxa_names with this DNAseq object
names(dna) <- taxa_names(ps) 
# Right now taxa_names = DNAseq, but we want both so we can manipulate taxa_names while retaining the DNAseqs in their own table
# add the DNAseq table to the phyloseq object:
ps <- merge_phyloseq(ps, dna)
ps #note the new "refseq()" part of the phyloseq object!
#
# IF YOU DON'T HAVE DNA SEQs SKIP TO HERE
# replace whatever is in the taxa_names space with ASV1, ASV2, ASV3, etc.
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# To run decontam, we will need a column ("is.neg") in the Metadata Table
# is.neg = control samples are indicated as "TRUE" and all other samples are "FALSE"
# Create the is.neg column:
is.neg <- sample_data(ps)$Who == "Control"
head(sample_data(ps)) #check that it looks correct

#
## RUN DECONTAM!
# outputs a table of stats for all taxa and decontam's decision about whether or not each taxon is a contaminant
# play around with changing the "threshold" value..
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold = 0.5)
# Summary table: how many true contaminants were identified?
table(contamdf.prev$contaminant) 
#tresh0.1 = 190contaminants, tresh0.05 = 93contaminants

# plot a histogram of the # of ASVs (y-axis) for each decontam score (a.k.a. P-value) on the x-axis
# ASVs are rows in the contamdf.prev object
ggplot(contamdf.prev, aes(x=p)) + 
  geom_histogram(binwidth=0.01, fill="purple")+
  scale_x_continuous(breaks=seq(0, 1, 0.1))
#there isn't much of a difference between 0 and 0.55...
#moving forward with 0.1, since that's the recommended threshold
# ...though I could concievably increase the threshold substantially...
# thresh0.5 = 746contaminants (14470 non-contminants)

### NOTES!
## "method=prevalence" means contaminants are identified by increased prevalence in negative controls vs. samples
# default: threshold = 0.1  ...defines the threshold for deciding what is more prevalent in the negative controls than the soil samples
# default: detailed = TRUE ...returns a data.frame with diagnostic info on the contamination decision
### Some notes about the threshold argument from Davis et al 2018:
## For P=0.5, ASVs would be classified as contaminants if present in a higher fraction of controls than samples
## The P-value is the result of 2x2 chi-squared test (many samples), or Fisher's Exact Test (few samples)
## this P-value is also refered to as a "score" or "decontam score"
## Each ASV is assigned a P-value or score
## In general, the threshold P value should be less-than-or-equal-to 0.5
## higher scores = likely not a contaminant

predecontamreads <- data.frame( SeqID = row.names(otu_table(ps)),
                                TotalReads = rowSums(otu_table(ps)))
write.csv(predecontamreads, "~/predecontamreads.csv")

sum(sample_data(ps)$Who == "Control", na.rm=TRUE) #46
sum(sample_data(ps)$Who == "Neemonika") #325
sum(sample_data(ps)$Who == "Neemonika", sample_data(ps)$Who =="Phillip", sample_data(ps)$Who =="PhillipPool") #618 ...618+46=664

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Who == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Who == c("Phillip", "PhillipPool", "Neemonika"), ps.pa)
#Warning: In sample_data(ps.pa)$Who == c("Phillip", "PhillipPool", "Neemonika") : "longer object length is not a multiple of shorter object length"

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
# Plot ASV presence in negative controls vs. positive samples and color by whether or not they were ID's as a contaminant
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant))+ 
  geom_point()+
  #geom_abline(intercept = 0, slope = 1)+ #optional line showing the threshold for equal representation in controls and samples
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
  ggtitle("Threshold = 0.1")


# REMOVE CONTAMINANTS FROM DATASET
contam <- rownames_to_column(contamdf.prev, var = "rowname") #move rownames (ASVs) to a column that data.table can work with
contam.dt <- data.table(contam) #convert to a data.table
badTaxa <- contam.dt[contaminant == "TRUE", rowname] #use data.table to generate a list of contaminant ASVs
goodTaxa <- setdiff(taxa_names(ps), badTaxa) #list of all ASVs minus the contaminant ASVs
ps2 <- prune_taxa(goodTaxa, ps) #create a new phyloseq object that only contains the good (non-contaminant) ASVs!
#check out the difference:
ps #15216 taxa
ps2 #15026 taxa
15216-15026 # =190, which is equivalent to the number of taxa that decontam identified as bad! this worked! hooray!


# Export final ps2 tables:
write.csv(tax_table(ps2), "~/ps2_TAXtable_withASVs.csv")
write.csv(otu_table(ps2), "~/ps2_OTUtable_withASVs.csv")
write.csv(refseq(ps2), "~/ps2_SEQtable_withASVs.csv")



#############################################
#### EXTRACT DATA FROM PHYLOSEQ OBJECTS #####
#############################################
# Make a copy of sample metadata table with total read counts/sample
MetadataWithReads <- copy(metasub2)
head(MetadataWithReads)
class(MetadataWithReads) #datatable and dataframe
# if you want to re-extract the sample_data() from a phyloseq object
# as.data.frame() or as.data.table() don't work for some reason, indead, use as():
SampDat.dt <- as.data.table(as(sample_data(ps2), "data.frame"))
colnames(otu.dft[1])
# Extract OTU table:
otu.df <- as.data.frame(otu_table(ps2.Neemonika))
otu.df <- as.data.frame(otu_table(ps2.NeemonTop))
otu.df[1:5,1:5]
dim(otu.df)
otu.df.rownames <- rownames_to_column(otu.df, var = "SeqID")
otu.df.rownames[1:5,1:5]
dim(otu.df)
#transpose OTU table
otu.dft <- as.data.frame(t(otu.df.rownames[,2:ncol(otu.df.rownames)]))
colnames(otu.dft) <- otu.df.rownames[,1] 
otu.dft[1:5,1:5]
otu.dft <- rownames_to_column(otu.dft, var = "ASV")
otu.dft[1:5,1:5]

# Extract taxonomy table:
tax.df <- as.data.frame(tax_table(ps2))
tax.df <- rownames_to_column(tax.df, var = "ASV")
head(tax.df)

# Extract Sequence table:
seq.df <- as.data.frame(refseq(ps2))
seq.df <- rownames_to_column(seq.df, var = "ASV")
names(seq.df)[2] <- "Sequence"

# Merge taxonomy with OTU table (= "spe" table in Numerical Ecology)
OTUtax.df <- merge(otu.dft, tax.df, by="ASV", all.x=TRUE)
dim(OTUtax.df)
OTUtax.df[1:5, 1:5]
OTUtax.df[1:5, 662:672]
unique(OTUtax.df$Phylum)
write.csv(OTUtax.df, "~/OTUtable_with_Taxonomy.csv")

#explore phyla...
OTUtax.dt <- as.data.table(OTUtax.df)
nrow(OTUtax.dt[Phylum == "NA"]) #0
nrow(OTUtax.dt[Phylum == "p__Monoblepharomycota"]) #46
nrow(OTUtax.dt[Phylum == "p__Basidiobolomycota"]) #5
nrow(OTUtax.dt[Phylum == "p__Blastocladiomycota"]) #4
nrow(OTUtax.dt[Phylum == "p__Zoopagomycota"]) #26
nrow(OTUtax.dt[Phylum == "p__Olpidiomycota"]) #6
nrow(OTUtax.dt[Phylum == "p__Kickxellomycota"]) #37
nrow(OTUtax.dt[Phylum == "p__Rozellomycota"]) #47
nrow(OTUtax.dt[Phylum == "p__Entorrhizomycota"]) #1
nrow(OTUtax.dt[Phylum == "p__Aphelidiomycota"]) #1
nrow(OTUtax.dt[Phylum == "p__Ascomycota"]) #8346
nrow(OTUtax.dt[Phylum == "p__Basidiomycota"]) #4524
nrow(OTUtax.dt[Phylum == "p__Mucoromycota"]) #285
nrow(OTUtax.dt[Phylum == "p__Glomeromycota"]) #191
nrow(OTUtax.dt[Phylum == "p__Chytridiomycota"]) #316

dim(OTUtax.df) #15026 taxa

length(unique(OTUtax.df$Kingdom)) #1
length(unique(OTUtax.df$Phylum)) #16
length(unique(OTUtax.df$Class)) #59
length(unique(OTUtax.df$Order)) #148
length(unique(OTUtax.df$Family)) #324
length(unique(OTUtax.df$Genus)) #708

# Merge Sequences with OTUtax table:
OTUtaxSeq.df <- merge(OTUtax.df, seq.df, by="ASV", all.x=TRUE)
OTUtaxSeq.df[1:5, 1:5]
OTUtaxSeq.df[1:5, 665:673]
write.csv(OTUtaxSeq.df, "~/OTUtable_withTaxonomy_and_Sequences.csv")

##########################################################
#### Nine Pyrophiles - Abundace over Time - Figure 6D ####
##########################################################

#subset for comp400 and set rownames to dateSamp
BS400otuITS16S <- ITS16Sotu[ITS16Sotu$SampID %in% Metadata400$SampID, ]
BS400otuITS16S.numdate <- merge(BS400otuITS16S, Metadata400[,c(5,26)], by="SampID", all.x=TRUE) #add dateSamp column
BS400otuITS16S.numdate$SampID <-NULL
BS400otuITS16S.numdate <- column_to_rownames(BS400otuITS16S.numdate, var="dateSamp")

#transpose, subset for taxa-of-interest
BS400otuITS16S.numdate.t <- t(BS400otuITS16S.numdate)
BS400otuITS16S.numdate.t[1:5, 1:5]


TAXtableITS16S[110130:110138, 1:7]
# subset for my favorite genera:
#generalist <- c("Bacillus", "Flavobacterium", "Massilia", "Pseudomonas", "Segetibacter")
generalist <- c("Pyronema", "Geopyxis", "Lyophyllum", "Myxomphalia", "Rhodosporidiobolus", 
                "Cellvibrio", "Bacillus", "Flavobacterium", "Massilia")
taxITS16S.pyros <- TAXtableITS16S[TAXtableITS16S$Genus %in% generalist,]
head(taxITS16S.pyros)
# subset BS400 OTUtable for only the taxa-of-interest
BS400otuITS16S.numdate.t.pyros <- BS400otuITS16S.numdate.t[rownames(BS400otuITS16S.numdate.t) %in% taxITS16S.pyros$ASV, ]
BS400otuITS16S.numdate.t.pyros[1:5,1:5]
#remove samples(columns) that sum to zero:
BS400otuITS16S.numdate.t.pyros <- BS400otuITS16S.numdate.t.pyros[ ,colSums(BS400otuITS16S.numdate.t.pyros) > 0]
BS400otuITS16S.numdate.t.pyros.m <- melt(BS400otuITS16S.numdate.t.pyros)
names(BS400otuITS16S.numdate.t.pyros.m) <- c("ASV", "dateSamp", "Value")
head(BS400otuITS16S.numdate.t.pyros.m)
#merge genus names onto table:
BS400ITS16S.Dat4Plot <- data.table(na.omit(merge(BS400otuITS16S.numdate.t.pyros.m, TAXtableITS16S[,c(1,7)], by="ASV", all.x=TRUE)))
BS400ITS16S.Dat4Plot[ ,c("numdate", "SampID"):=tstrsplit(dateSamp, "_", fixed=FALSE)]
BS400ITS16S.Dat4Plot <- BS400ITS16S.Dat4Plot[-1804, ] #remove bizarrely high Pyronema pre-fire...

#normalize each ASV to the max value for it's genus:
# print the max value for each genus:
BS400ITS16S.Dat4Plot[ , max(Value), by = Genus]
BS400ITS16S.Dat4Plot[grepl("Myxomphalia", Genus), NormVal:=Value/1512]#
BS400ITS16S.Dat4Plot[grepl("Lyophyllum", Genus), NormVal:=Value/777]#
BS400ITS16S.Dat4Plot[grepl("Pyronema", Genus), NormVal:=Value/1629]
BS400ITS16S.Dat4Plot[grepl("Rhodosporidiobolus", Genus), NormVal:=Value/319]#
BS400ITS16S.Dat4Plot[grepl("Geopyxis", Genus), NormVal:=Value/1273]#
BS400ITS16S.Dat4Plot[grepl("Massilia", Genus), NormVal:=Value/4507]
BS400ITS16S.Dat4Plot[grepl("Flavobacterium", Genus), NormVal:=Value/1607]
BS400ITS16S.Dat4Plot[grepl("Bacillus", Genus), NormVal:=Value/734]
BS400ITS16S.Dat4Plot[grepl("Cellvibrio", Genus), NormVal:=Value/230]

BS400ITS16S.Dat4Plot$numdate <- as.numeric(BS400ITS16S.Dat4Plot$numdate) 
#note that the x-axis has to be numeric for things to work like geom_smooth() or filter
BS400ITS16S.Dat4Plot$Value <- as.numeric(BS400ITS16S.Dat4Plot$Value)

BS400ITS16S.Dat4Plot.top <-  BS400ITS16S.Dat4Plot %>% group_by(.dots=c("numdate","Genus")) %>% top_n(1, NormVal)
BS400ITS16S.Dat4Plot.top.filtered <- unique(BS400ITS16S.Dat4Plot.top)
view(BS400ITS16S.Dat4Plot.top.filtered)
# line plot of Blodgett data:
ggplot(BS400ITS16S.Dat4Plotmeans, aes(x=numdate, y=NormVal, color=Genus))+
  #geom_point(aes(shape=Genus), size=3)+
  geom_line(size=1)+
  #geom_smooth(aes(fill=Genus),  size=1, alpha=0, method="loess", span=0.5)+
  #scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0, 1, 2))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                "#6de2cc", "#a65628", "#f781bf", "#999999"))+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12, angle = 90,  hjust=0, vjust=0.5))+
  theme(axis.title = element_text(size=14, face="bold"))+
  theme(legend.text = element_text(size=10, face="bold"))+
  theme(legend.title = element_text(size=10, face="bold"))+
  xlab("Time")+
  ylab("Normalized Peak Abundance")+
  ggtitle("Hi Burn (ITS & 16S)")


#normalize each ASV to the max AVERAGE value for it's genus:
## merge genus names onto table:
BS400ITS16S.Dat4Plot <- data.table(na.omit(merge(BS400otuITS16S.numdate.t.pyros.m, TAXtableITS16S[,c(1,7)], by="ASV", all.x=TRUE)))
BS400ITS16S.Dat4Plot[ ,c("numdate", "SampID"):=tstrsplit(dateSamp, "_", fixed=FALSE)]
BS400ITS16S.Dat4Plot <- BS400ITS16S.Dat4Plot[-1804, ] #remove bizarrely high Pyronema pre-fire...
# Average value for each Genus&numdate:
BS400ITS16S.Dat4Plotmeans <- BS400ITS16S.Dat4Plot[ , mean(Value), by = c("Genus", "numdate")]
names(BS400ITS16S.Dat4Plotmeans)[3] <- "mean"
view(BS400ITS16S.Dat4Plotmeans)
#normalize Averages to the max average for the genus:
BS400ITS16S.Dat4Plotmeans[ , max(mean), by = Genus]
BS400ITS16S.Dat4Plotmeans[grepl("Myxomphalia", Genus), NormVal:=mean/169.600]#
BS400ITS16S.Dat4Plotmeans[grepl("Lyophyllum", Genus), NormVal:=mean/6.2361111]#
BS400ITS16S.Dat4Plotmeans[grepl("Pyronema", Genus), NormVal:=mean/92.91666667]
BS400ITS16S.Dat4Plotmeans[grepl("Rhodosporidiobolus", Genus), NormVal:=mean/22.2000000]#
BS400ITS16S.Dat4Plotmeans[grepl("Geopyxis", Genus), NormVal:=mean/74.3333333]#
BS400ITS16S.Dat4Plotmeans[grepl("Massilia", Genus), NormVal:=mean/7.11771654]
BS400ITS16S.Dat4Plotmeans[grepl("Flavobacterium", Genus), NormVal:=mean/4.8946809]
BS400ITS16S.Dat4Plotmeans[grepl("Bacillus", Genus), NormVal:=mean/5.74418605]
BS400ITS16S.Dat4Plotmeans[grepl("Cellvibrio", Genus), NormVal:=mean/2.82857143]

view(BS400ITS16S.Dat4Plotmeans)
#
BS400ITS16S.Dat4Plotmeans$numdate <- as.numeric(BS400ITS16S.Dat4Plotmeans$numdate) 
#note that the x-axis has to be numeric for things to work like geom_smooth() or filter
BS400ITS16S.Dat4Plotmeans$Value <- as.numeric(BS400ITS16S.Dat4Plotmeans$Value)
#replace values greater than 1 with 1
BS400ITS16S.Dat4Plotmeans[,4] <- lapply(BS400ITS16S.Dat4Plotmeans[,4], function(x) ifelse(x > 1, 1, x))

# line plot of Blodgett data:
ggplot(BS400ITS16S.Dat4Plotmeans, aes(x=numdate, y=NormVal, color=Genus))+
  geom_point(shape=16, size=2)+
  geom_line(size=0.5)+
  #geom_smooth(aes(fill=Genus),  size=1, alpha=0, method="loess", span=0.5)+
  #scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0, 1, 2))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                "#6de2cc", "#a65628", "#f781bf", "#999999"))+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12, angle = 90,  hjust=0, vjust=0.5))+
  theme(axis.title = element_text(size=14, face="bold"))+
  theme(legend.text = element_text(size=14, face="bold"))+
  theme(legend.title = element_text(size=14, face="bold"))+
  xlab("Time")+
  ylab("Normalized Average Abundance")+
  ggtitle("Hi Burn (ITS & 16S)")
#export PDF = 8x4
#
##
#
##
#
## Lo Burn ITS & 16S!
#subset for comp400 and set rownames to dateSamp
BS321EotuITS16S <- ITS16Sotu[ITS16Sotu$SampID %in% Metadata321E$SampID, ]
BS321EotuITS16S.numdate <- merge(BS321EotuITS16S, Metadata321E[,c(5,26)], by="SampID", all.x=TRUE) #add dateSamp column
BS321EotuITS16S.numdate$SampID <-NULL
BS321EotuITS16S.numdate <- column_to_rownames(BS321EotuITS16S.numdate, var="dateSamp")

#transpose, subset for taxa-of-interest
BS321EotuITS16S.numdate.t <- t(BS321EotuITS16S.numdate)
BS321EotuITS16S.numdate.t[1:5, 1:5]


TAXtableITS16S[110130:110138, 1:7]
# subset for my favorite genera:
#generalist <- c("Bacillus", "Flavobacterium", "Massilia", "Pseudomonas", "Segetibacter")
generalist <- c("Pyronema", "Geopyxis", "Lyophyllum", "Myxomphalia", "Rhodosporidiobolus", 
                "Cellvibrio", "Bacillus", "Flavobacterium", "Massilia")
taxITS16S.pyros <- TAXtableITS16S[TAXtableITS16S$Genus %in% generalist,]
head(taxITS16S.pyros)
# subset BS400 OTUtable for only the taxa-of-interest
BS321EotuITS16S.numdate.t.pyros <- BS321EotuITS16S.numdate.t[rownames(BS321EotuITS16S.numdate.t) %in% taxITS16S.pyros$ASV, ]
#remove columns that sum to zero:
BS321EotuITS16S.numdate.t.pyros <- BS321EotuITS16S.numdate.t.pyros[ ,colSums(BS321EotuITS16S.numdate.t.pyros) > 0]
BS321EotuITS16S.numdate.t.pyros.m <- melt(BS321EotuITS16S.numdate.t.pyros)
names(BS321EotuITS16S.numdate.t.pyros.m) <- c("ASV", "dateSamp", "Value")
head(BS321EotuITS16S.numdate.t.pyros.m)
BS321EITS16S.Dat4Plot <- data.table(na.omit(merge(BS321EotuITS16S.numdate.t.pyros.m, TAXtableITS16S[,c(1,7)], by="ASV", all.x=TRUE)))
BS321EITS16S.Dat4Plot[ ,c("numdate", "SampID"):=tstrsplit(dateSamp, "_", fixed=FALSE)]
#BS321EITS16S.Dat4Plot <- BS321EITS16S.Dat4Plot[-1804, ] #remove bizarrely high Pyronema pre-fire...

#normalize each ASV to the max value for it's genus:
# print the max value for each genus:
BS321EITS16S.Dat4Plot[ , max(Value), by = Genus]
BS321EITS16S.Dat4Plot[grepl("Myxomphalia", Genus), NormVal:=Value/1512]#
BS321EITS16S.Dat4Plot[grepl("Lyophyllum", Genus), NormVal:=Value/777]#
BS321EITS16S.Dat4Plot[grepl("Pyronema", Genus), NormVal:=Value/1629]
BS321EITS16S.Dat4Plot[grepl("Rhodosporidiobolus", Genus), NormVal:=Value/319]#
BS321EITS16S.Dat4Plot[grepl("Geopyxis", Genus), NormVal:=Value/1273]#
BS321EITS16S.Dat4Plot[grepl("Massilia", Genus), NormVal:=Value/4507]
BS321EITS16S.Dat4Plot[grepl("Flavobacterium", Genus), NormVal:=Value/1607]
BS321EITS16S.Dat4Plot[grepl("Bacillus", Genus), NormVal:=Value/734]
BS321EITS16S.Dat4Plot[grepl("Cellvibrio", Genus), NormVal:=Value/230]
# remove taxa that fall below the 5th Percentile for all abundance values = 10

view(BS321EITS16S.Dat4Plot)
#
BS321EITS16S.Dat4Plot$numdate <- as.numeric(BS321EITS16S.Dat4Plot$numdate) 
#note that the x-axis has to be numeric for things to work like geom_smooth() or filter
BS321EITS16S.Dat4Plot$Value <- as.numeric(BS321EITS16S.Dat4Plot$Value)


BS321EITS16S.Dat4Plot.top <-  BS321EITS16S.Dat4Plot %>% group_by(.dots=c("numdate","Genus")) %>% top_n(1, NormVal)
BS321EITS16S.Dat4Plot.top.filtered <- unique(BS321EITS16S.Dat4Plot.top[,c(3:5,7)])

# line plot of Blodgett data:
ggplot(BS321EITS16S.Dat4Plot.top.filtered, aes(x=numdate, y=NormVal, color=Genus))+
  #geom_point(aes(shape=Genus), size=3)+
  geom_line(size=1)+
  #geom_smooth(aes(fill=Genus),  size=1, alpha=0, method="loess", span = 0.5)+
  #scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                "#6de2cc", "#a65628", "#f781bf", "#999999"))+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12, angle = 90,  hjust=0, vjust=0.5))+
  theme(axis.title = element_text(size=14, face="bold"))+
  theme(legend.text = element_text(size=14, face="bold"))+
  theme(legend.title = element_text(size=16, face="bold"))+
  xlab("Time")+
  ylab("Normalized Peak Abundance")+
  ggtitle("Lo Burn (ITS & 16S)")



#normalize each ASV to the max AVERAGE value for it's genus:
## merge genus names onto table:
BS321EITS16S.Dat4Plot <- data.table(na.omit(merge(BS321EotuITS16S.numdate.t.pyros.m, TAXtableITS16S[,c(1,7)], by="ASV", all.x=TRUE)))
BS321EITS16S.Dat4Plot[ ,c("numdate", "SampID"):=tstrsplit(dateSamp, "_", fixed=FALSE)]
#BS321EITS16S.Dat4Plot <- BS321EITS16S.Dat4Plot[-1804, ] #remove bizarrely high Pyronema pre-fire...
# Average value for each Genus&numdate:
BS321EITS16S.Dat4Plotmeans <- BS321EITS16S.Dat4Plot[ , mean(Value), by = c("Genus", "numdate")]
names(BS321EITS16S.Dat4Plotmeans)[3] <- "mean"
view(BS321EITS16S.Dat4Plotmeans)
#normalize Averages to the max average for the genus:
BS321EITS16S.Dat4Plotmeans[ , max(mean), by = Genus]
BS321EITS16S.Dat4Plotmeans[grepl("Myxomphalia", Genus), NormVal:=mean/169.600]#
BS321EITS16S.Dat4Plotmeans[grepl("Lyophyllum", Genus), NormVal:=mean/6.2361111]#
BS321EITS16S.Dat4Plotmeans[grepl("Pyronema", Genus), NormVal:=mean/92.91666667]
BS321EITS16S.Dat4Plotmeans[grepl("Rhodosporidiobolus", Genus), NormVal:=mean/22.2000000]#
BS321EITS16S.Dat4Plotmeans[grepl("Geopyxis", Genus), NormVal:=mean/74.3333333]#
BS321EITS16S.Dat4Plotmeans[grepl("Massilia", Genus), NormVal:=mean/7.11771654]
BS321EITS16S.Dat4Plotmeans[grepl("Flavobacterium", Genus), NormVal:=mean/4.8946809]
BS321EITS16S.Dat4Plotmeans[grepl("Bacillus", Genus), NormVal:=mean/5.74418605]
BS321EITS16S.Dat4Plotmeans[grepl("Cellvibrio", Genus), NormVal:=mean/2.82857143]

view(BS321EITS16S.Dat4Plotmeans)
#
BS321EITS16S.Dat4Plotmeans$numdate <- as.numeric(BS321EITS16S.Dat4Plotmeans$numdate) 
#note that the x-axis has to be numeric for things to work like geom_smooth() or filter
BS321EITS16S.Dat4Plotmeans$Value <- as.numeric(BS321EITS16S.Dat4Plotmeans$Value)
#replace values greater than 1 with 1
BS321EITS16S.Dat4Plotmeans[,4] <- lapply(BS321EITS16S.Dat4Plotmeans[,4], function(x) ifelse(x > 1, 1, x))

# line plot of Blodgett data:
ggplot(BS321EITS16S.Dat4Plotmeans, aes(x=numdate, y=NormVal, color=Genus))+
  geom_point(shape=16, size=2)+
  geom_line(size=0.5)+
  #geom_smooth(aes(fill=Genus),  size=1, alpha=0, method="loess", span=0.5)+
  scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0, 1, 2))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                "#6de2cc", "#a65628", "#f781bf", "#999999"))+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12, angle = 90,  hjust=0, vjust=0.5))+
  theme(axis.title = element_text(size=14, face="bold"))+
  theme(legend.text = element_text(size=14, face="bold"))+
  theme(legend.title = element_text(size=14, face="bold"))+
  xlab("Time")+
  ylab("Normalized Average Abundance")+
  ggtitle("Lo Burn (ITS & 16S)")
#export PDF = 8x4
#
##
#
##
#
## CONTROL #2 ITS & 16S!
#subset for comp400 and set rownames to dateSamp
BS321WotuITS16S <- ITS16Sotu[ITS16Sotu$SampID %in% Metadata321W$SampID, ]
BS321WotuITS16S.numdate <- merge(BS321WotuITS16S, Metadata321W[,c(5,26)], by="SampID", all.x=TRUE) #add dateSamp column
BS321WotuITS16S.numdate$SampID <-NULL
BS321WotuITS16S.numdate <- column_to_rownames(BS321WotuITS16S.numdate, var="dateSamp")

#transpose, subset for taxa-of-interest
BS321WotuITS16S.numdate.t <- t(BS321WotuITS16S.numdate)
BS321WotuITS16S.numdate.t[1:5, 1:5]


TAXtableITS16S[110130:110138, 1:7]
# subset for my favorite genera:
#generalist <- c("Bacillus", "Flavobacterium", "Massilia", "Pseudomonas", "Segetibacter")
generalist <- c("Pyronema", "Geopyxis", "Lyophyllum", "Myxomphalia", "Rhodosporidiobolus", 
                "Cellvibrio", "Bacillus", "Flavobacterium", "Massilia")
taxITS16S.pyros <- TAXtableITS16S[TAXtableITS16S$Genus %in% generalist,]
head(taxITS16S.pyros)
# subset BS400 OTUtable for only the taxa-of-interest
BS321WotuITS16S.numdate.t.pyros <- BS321WotuITS16S.numdate.t[rownames(BS321WotuITS16S.numdate.t) %in% taxITS16S.pyros$ASV, ]
#remove columns that sum to zero:
BS321WotuITS16S.numdate.t.pyros <- BS321WotuITS16S.numdate.t.pyros[ ,colSums(BS321WotuITS16S.numdate.t.pyros) > 0]
BS321WotuITS16S.numdate.t.pyros.m <- melt(BS321WotuITS16S.numdate.t.pyros)
names(BS321WotuITS16S.numdate.t.pyros.m) <- c("ASV", "dateSamp", "Value")
head(BS321WotuITS16S.numdate.t.pyros.m)
BS321WITS16S.Dat4Plot <- data.table(na.omit(merge(BS321WotuITS16S.numdate.t.pyros.m, TAXtableITS16S[,c(1,7)], by="ASV", all.x=TRUE)))
BS321WITS16S.Dat4Plot[ ,c("numdate", "SampID"):=tstrsplit(dateSamp, "_", fixed=FALSE)]
#BS321EITS16S.Dat4Plot <- BS321EITS16S.Dat4Plot[-1804, ] #remove bizarrely high Pyronema pre-fire...

#normalize each ASV to the max value for it's genus:
# print the max value for each genus:
BS321WITS16S.Dat4Plot[ , max(Value), by = Genus]
BS321WITS16S.Dat4Plot[grepl("Myxomphalia", Genus), NormVal:=Value/1512]#
BS321WITS16S.Dat4Plot[grepl("Lyophyllum", Genus), NormVal:=Value/777]#
BS321WITS16S.Dat4Plot[grepl("Pyronema", Genus), NormVal:=Value/1629]
BS321WITS16S.Dat4Plot[grepl("Rhodosporidiobolus", Genus), NormVal:=Value/319]#
BS321WITS16S.Dat4Plot[grepl("Geopyxis", Genus), NormVal:=Value/1273]#
BS321WITS16S.Dat4Plot[grepl("Massilia", Genus), NormVal:=Value/4507]
BS321WITS16S.Dat4Plot[grepl("Flavobacterium", Genus), NormVal:=Value/1607]
BS321WITS16S.Dat4Plot[grepl("Bacillus", Genus), NormVal:=Value/734]
BS321WITS16S.Dat4Plot[grepl("Cellvibrio", Genus), NormVal:=Value/230]
# remove taxa that fall below the 5th Percentile for all abundance values = 10

view(BS321WITS16S.Dat4Plot)
#
BS321WITS16S.Dat4Plot$numdate <- as.numeric(BS321WITS16S.Dat4Plot$numdate) 
#note that the x-axis has to be numeric for things to work like geom_smooth() or filter
BS321WITS16S.Dat4Plot$Value <- as.numeric(BS321WITS16S.Dat4Plot$Value)


BS321WITS16S.Dat4Plot.top <-  BS321WITS16S.Dat4Plot %>% group_by(.dots=c("numdate","Genus")) %>% top_n(1, NormVal)
BS321WITS16S.Dat4Plot.top.filtered <- unique(BS321WITS16S.Dat4Plot.top[,c(3:5,7)])

# line plot of Blodgett data:
ggplot(BS321WITS16S.Dat4Plot.top.filtered, aes(x=numdate, y=NormVal, color=Genus))+
  #geom_point(aes(shape=Genus), size=3)+
  geom_line(size=1)+
  #geom_smooth(aes(fill=Genus),  size=1, alpha=0, method="loess", span = 0.5)+
  #scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                "#6de2cc", "#a65628", "#f781bf", "#999999"))+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12, angle = 90,  hjust=0, vjust=0.5))+
  theme(axis.title = element_text(size=14, face="bold"))+
  theme(legend.text = element_text(size=14, face="bold"))+
  theme(legend.title = element_text(size=16, face="bold"))+
  xlab("Time")+
  ylab("Normalized Peak Abundance")+
  ggtitle("2C Burn (ITS & 16S)")

#normalize each ASV to the max AVERAGE value for it's genus:
## merge genus names onto table:
BS321WITS16S.Dat4Plot <- data.table(na.omit(merge(BS321WotuITS16S.numdate.t.pyros.m, TAXtableITS16S[,c(1,7)], by="ASV", all.x=TRUE)))
BS321WITS16S.Dat4Plot[ ,c("numdate", "SampID"):=tstrsplit(dateSamp, "_", fixed=FALSE)]
#BS321WITS16S.Dat4Plot <- BS321WITS16S.Dat4Plot[-1804, ] #remove bizarrely high Pyronema pre-fire...
# Average value for each Genus&numdate:
BS321WITS16S.Dat4Plotmeans <- BS321WITS16S.Dat4Plot[ , mean(Value), by = c("Genus", "numdate")]
names(BS321WITS16S.Dat4Plotmeans)[3] <- "mean"
head(BS321WITS16S.Dat4Plotmeans)
#normalize Averages to the max average for the genus:
BS321WITS16S.Dat4Plotmeans[ , max(mean), by = Genus]
BS321WITS16S.Dat4Plotmeans[grepl("Myxomphalia", Genus), NormVal:=mean/169.600]#
BS321WITS16S.Dat4Plotmeans[grepl("Lyophyllum", Genus), NormVal:=mean/6.2361111]#
BS321WITS16S.Dat4Plotmeans[grepl("Pyronema", Genus), NormVal:=mean/92.91666667]
BS321WITS16S.Dat4Plotmeans[grepl("Rhodosporidiobolus", Genus), NormVal:=mean/22.2000000]#
BS321WITS16S.Dat4Plotmeans[grepl("Geopyxis", Genus), NormVal:=mean/74.3333333]#
BS321WITS16S.Dat4Plotmeans[grepl("Massilia", Genus), NormVal:=mean/7.11771654]
BS321WITS16S.Dat4Plotmeans[grepl("Flavobacterium", Genus), NormVal:=mean/4.8946809]
BS321WITS16S.Dat4Plotmeans[grepl("Bacillus", Genus), NormVal:=mean/5.74418605]
BS321WITS16S.Dat4Plotmeans[grepl("Cellvibrio", Genus), NormVal:=mean/2.82857143]

#
BS321WITS16S.Dat4Plotmeans$numdate <- as.numeric(BS321WITS16S.Dat4Plotmeans$numdate) 
#note that the x-axis has to be numeric for things to work like geom_smooth() or filter
BS321WITS16S.Dat4Plotmeans$Value <- as.numeric(BS321WITS16S.Dat4Plotmeans$Value)
#replace values greater than 1 with 1
BS321WITS16S.Dat4Plotmeans[,4] <- lapply(BS321WITS16S.Dat4Plotmeans[,4], function(x) ifelse(x > 1, 1, x))

# line plot of Blodgett data:
ggplot(BS321WITS16S.Dat4Plotmeans, aes(x=numdate, y=NormVal, color=Genus))+
  geom_point(shape=16, size=2)+
  geom_line(size=0.5)+
  #geom_smooth(aes(fill=Genus),  size=1, alpha=0, method="loess", span=0.5)+
  scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0, 1, 2))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                "#6de2cc", "#a65628", "#f781bf", "#999999"))+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12, angle = 90,  hjust=0, vjust=0.5))+
  theme(axis.title = element_text(size=14, face="bold"))+
  theme(legend.text = element_text(size=14, face="bold"))+
  theme(legend.title = element_text(size=14, face="bold"))+
  xlab("Time")+
  ylab("Normalized Average Abundance")+
  ggtitle("Control #2 (ITS & 16S)")
#export PDF = 8x4


#
##
#
##
#
## CONTROL #1 ITS & 16S!
#subset for comp400 and set rownames to dateSamp
BS240otuITS16S <- ITS16Sotu[ITS16Sotu$SampID %in% Metadata240$SampID, ]
BS240otuITS16S.numdate <- merge(BS240otuITS16S, Metadata240[,c(5,26)], by="SampID", all.x=TRUE) #add dateSamp column
BS240otuITS16S.numdate$SampID <-NULL
BS240otuITS16S.numdate <- column_to_rownames(BS240otuITS16S.numdate, var="dateSamp")

#transpose, subset for taxa-of-interest
BS240otuITS16S.numdate.t <- t(BS240otuITS16S.numdate)
BS240otuITS16S.numdate.t[1:5, 1:5]


TAXtableITS16S[110130:110138, 1:7]
# subset for my favorite genera:
#generalist <- c("Bacillus", "Flavobacterium", "Massilia", "Pseudomonas", "Segetibacter")
generalist <- c("Pyronema", "Geopyxis", "Lyophyllum", "Myxomphalia", "Rhodosporidiobolus", 
                "Cellvibrio", "Bacillus", "Flavobacterium", "Massilia")
taxITS16S.pyros <- TAXtableITS16S[TAXtableITS16S$Genus %in% generalist,]
head(taxITS16S.pyros)
# subset BS400 OTUtable for only the taxa-of-interest
BS240otuITS16S.numdate.t.pyros <- BS240otuITS16S.numdate.t[rownames(BS240otuITS16S.numdate.t) %in% taxITS16S.pyros$ASV, ]
#remove columns that sum to zero:
BS240otuITS16S.numdate.t.pyros <- BS240otuITS16S.numdate.t.pyros[ ,colSums(BS240otuITS16S.numdate.t.pyros) > 0]
BS240otuITS16S.numdate.t.pyros.m <- melt(BS240otuITS16S.numdate.t.pyros)
names(BS240otuITS16S.numdate.t.pyros.m) <- c("ASV", "dateSamp", "Value")
head(BS240otuITS16S.numdate.t.pyros.m)
BS240ITS16S.Dat4Plot <- data.table(na.omit(merge(BS240otuITS16S.numdate.t.pyros.m, TAXtableITS16S[,c(1,7)], by="ASV", all.x=TRUE)))
BS240ITS16S.Dat4Plot[ ,c("numdate", "SampID"):=tstrsplit(dateSamp, "_", fixed=FALSE)]
#BS321EITS16S.Dat4Plot <- BS321EITS16S.Dat4Plot[-1804, ] #remove bizarrely high Pyronema pre-fire...

#normalize each ASV to the max value for it's genus:
# print the max value for each genus:
BS240ITS16S.Dat4Plot[ , max(Value), by = Genus]
BS240ITS16S.Dat4Plot[grepl("Myxomphalia", Genus), NormVal:=Value/1512]#
BS240ITS16S.Dat4Plot[grepl("Lyophyllum", Genus), NormVal:=Value/777]#
BS240ITS16S.Dat4Plot[grepl("Pyronema", Genus), NormVal:=Value/1629]
BS240ITS16S.Dat4Plot[grepl("Rhodosporidiobolus", Genus), NormVal:=Value/319]#
BS240ITS16S.Dat4Plot[grepl("Geopyxis", Genus), NormVal:=Value/1273]#
BS240ITS16S.Dat4Plot[grepl("Massilia", Genus), NormVal:=Value/4507]
BS240ITS16S.Dat4Plot[grepl("Flavobacterium", Genus), NormVal:=Value/1607]
BS240ITS16S.Dat4Plot[grepl("Bacillus", Genus), NormVal:=Value/734]
BS240ITS16S.Dat4Plot[grepl("Cellvibrio", Genus), NormVal:=Value/230]
# remove taxa that fall below the 5th Percentile for all abundance values = 10

view(BS240ITS16S.Dat4Plot)
#
BS240ITS16S.Dat4Plot$numdate <- as.numeric(BS240ITS16S.Dat4Plot$numdate) 
#note that the x-axis has to be numeric for things to work like geom_smooth() or filter
BS321WITS16S.Dat4Plot$Value <- as.numeric(BS240ITS16S.Dat4Plot$Value)


BS240ITS16S.Dat4Plot.top <-  BS240ITS16S.Dat4Plot %>% group_by(.dots=c("numdate","Genus")) %>% top_n(1, NormVal)
BS240ITS16S.Dat4Plot.top.filtered <- unique(BS240ITS16S.Dat4Plot.top[,c(3:5,7)])

# line plot of Blodgett data:
ggplot(BS240ITS16S.Dat4Plot.top.filtered, aes(x=numdate, y=NormVal, color=Genus))+
  #geom_point(aes(shape=Genus), size=3)+
  geom_line(size=1)+
  #geom_smooth(aes(fill=Genus),  size=1, alpha=0, method="loess", span = 0.5)+
  #scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                "#6de2cc", "#a65628", "#f781bf", "#999999"))+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12, angle = 90,  hjust=0, vjust=0.5))+
  theme(axis.title = element_text(size=14, face="bold"))+
  theme(legend.text = element_text(size=14, face="bold"))+
  theme(legend.title = element_text(size=16, face="bold"))+
  xlab("Time")+
  ylab("Normalized Peak Abundance")+
  ggtitle("1C (ITS & 16S)")
#export PDF = 8x4


#normalize each ASV to the max AVERAGE value for it's genus:
## merge genus names onto table:
BS240ITS16S.Dat4Plot <- data.table(na.omit(merge(BS240otuITS16S.numdate.t.pyros.m, TAXtableITS16S[,c(1,7)], by="ASV", all.x=TRUE)))
BS240ITS16S.Dat4Plot[ ,c("numdate", "SampID"):=tstrsplit(dateSamp, "_", fixed=FALSE)]
#BS321WITS16S.Dat4Plot <- BS321WITS16S.Dat4Plot[-1804, ] #remove bizarrely high Pyronema pre-fire...
# Average value for each Genus&numdate:
BS240ITS16S.Dat4Plotmeans <- BS240ITS16S.Dat4Plot[ , mean(Value), by = c("Genus", "numdate")]
names(BS240ITS16S.Dat4Plotmeans)[3] <- "mean"
head(BS240ITS16S.Dat4Plotmeans)
#normalize Averages to the max average for the genus:
BS240ITS16S.Dat4Plotmeans[ , max(mean), by = Genus]
BS240ITS16S.Dat4Plotmeans[grepl("Myxomphalia", Genus), NormVal:=mean/169.600]#
BS240ITS16S.Dat4Plotmeans[grepl("Lyophyllum", Genus), NormVal:=mean/6.2361111]#
BS240ITS16S.Dat4Plotmeans[grepl("Pyronema", Genus), NormVal:=mean/92.91666667]
BS240ITS16S.Dat4Plotmeans[grepl("Rhodosporidiobolus", Genus), NormVal:=mean/22.2000000]#
BS240ITS16S.Dat4Plotmeans[grepl("Geopyxis", Genus), NormVal:=mean/74.3333333]#
BS240ITS16S.Dat4Plotmeans[grepl("Massilia", Genus), NormVal:=mean/7.11771654]
BS240ITS16S.Dat4Plotmeans[grepl("Flavobacterium", Genus), NormVal:=mean/4.8946809]
BS240ITS16S.Dat4Plotmeans[grepl("Bacillus", Genus), NormVal:=mean/5.74418605]
BS240ITS16S.Dat4Plotmeans[grepl("Cellvibrio", Genus), NormVal:=mean/2.82857143]

#
BS240ITS16S.Dat4Plotmeans$numdate <- as.numeric(BS240ITS16S.Dat4Plotmeans$numdate) 
#note that the x-axis has to be numeric for things to work like geom_smooth() or filter
BS240ITS16S.Dat4Plotmeans$Value <- as.numeric(BS240ITS16S.Dat4Plotmeans$Value)
#replace values greater than 1 with 1
BS240ITS16S.Dat4Plotmeans[,4] <- lapply(BS240ITS16S.Dat4Plotmeans[,4], function(x) ifelse(x > 1, 1, x))

# line plot of Blodgett data:
ggplot(BS240ITS16S.Dat4Plotmeans, aes(x=numdate, y=NormVal, color=Genus))+
  geom_point(shape=16, size=2)+
  geom_line(size=0.5)+
  #geom_smooth(aes(fill=Genus),  size=1, alpha=0, method="loess", span=0.5)+
  scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0, 1, 2))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
                                "#6de2cc", "#a65628", "#f781bf", "#999999"))+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12, angle = 90,  hjust=0, vjust=0.5))+
  theme(axis.title = element_text(size=14, face="bold"))+
  theme(legend.text = element_text(size=14, face="bold"))+
  theme(legend.title = element_text(size=14, face="bold"))+
  xlab("Time")+
  ylab("Normalized Average Abundance")+
  ggtitle("Control #1 (ITS & 16S)")
#export PDF = 8x4
######################################
#### Blodgett Temp & Precip Data! ####
######################################
library(data.table)
library(ggplot2)
BStime <- fread("~/TOTALBlodgettTempPrecipData_SeqDatesOnly.csv")
head(BStime)
# Create a new column that is a merge of Date and Hour columns:
BStime[ ,DateHour:=do.call(paste, c(.SD, sep="/")), .SDcols=c(1,5)]
head(BStime)
BStime <- na.omit(BStime)

BStime$DateHour <- as.Date(BStime$DateHour, "%m/%d/%y/%H")
head(BStime)

ggplot(data = BStime, mapping=aes(x=DateHour, y=Prec_in, group=1))+ 
  geom_bar(stat = "identity", color="blue", fill="blue", width = 0.5)+ 
  geom_line(mapping = aes(y = (TempC+30)/7), color="red", size=0.5)+ #transform Temp
  scale_x_date(date_breaks = "2 weeks", date_labels = "%e %b %Y")+
  scale_y_continuous(limits = c(0, 10), 
                     "Precipitation [inches]", 
                     sec.axis = sec_axis(~ .*7-30, name = "Temperature [°C]"))+ #undo the Temp transformation above
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color="grey90"),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.3), 
        axis.ticks = element_line(color="black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=12, color="black"),
        axis.text.x = element_text(size=12, hjust=1, vjust=0.5, angle=90),
        axis.text.y = element_text(size=12, color="black"))+
  xlab("Date")

### Summarize temp & precip data by date
# mean temp/day
BStemp.meansd <- BStime[ ,list(mean(TempC), 
                               mean(TempC)+sd(TempC), 
                               mean(TempC)-sd(TempC)),
                         by=DateHour, ]
colnames(BStemp.meansd) <- c("Date", "MeanTemp", "TempTopSD", "TempBottomSD")
BStemp.meansd <- na.omit(BStemp.meansd)

ggplot(BStemp.meansd, aes(x=Date, y=MeanTemp))+
  geom_point(shape=16, size=1)+
  geom_errorbar(aes(ymax=TempTopSD, ymin=TempBottomSD), width=0.1, alpha=0.2)+
  WhiteThemes

# sum of all precip/day
BSprecip <- BStime[ ,list(sum(Prec_in)), by=DateHour, ]

colnames(BSprecip) <- c("Date", "Precip")

ggplot(BSprecip, aes(Date, Precip))+
  geom_col()
WhiteThemes

ggplot(BSprecip, aes(Date, Precip))+
  geom_bar(stat="identity")

####################################
######### DIVERSITY METRICS ########
####################################
library(vegan)

## IS THERE A DIFFERENCE IN DIVERSITY BY DEPTH? -- Phillip's samples only
#sample_data(ps2)$Who == "Phillip" & sample_data(ps2)$Who == "PhillipPool" &
#sample_data(ps2)$Date== "17-Feb-20" & sample_data(ps2)$Date== 17-Feb-20

#subset for the two dates that are pre- and post- fire for Phillip's samples
unique(MetadataWithReads$Date)
MetadataPP <- MetadataWithReads[Date == "2/17/20" | Date == "10/8/19", ]
MetadataPP <- MetadataPP[Who == "Phillip" | Who == "PhillipPool"]
MetadataPP <- MetadataPP[Plot == "321west"]

otu.df <- as.data.frame(otu_table(ps2))
otu.df.rownames <- rownames_to_column(otu.df, var="SeqID")
dim(otu.df.rownames)

# Subset OTUtax table for only the column names that match the SeqIDs in MetadataPP
otu.df.rownames[1:5, 1:5]
OTUtaxPP <- otu.df.rownames[otu.df.rownames$SeqID %in% MetadataPP$SeqID, ]
OTUtaxPP[1:3, 1:5]
rownames(OTUtaxPP) <- c() #clear rownames
OTUtaxPP <- column_to_rownames(OTUtaxPP, var="SeqID") #move SeID column into rownames
unique(ShanDivPP$depth)
ShanDivPP <- as.data.table(diversity(OTUtaxPP, index="shannon"))
ShanDivPP[ ,numdate:=MetadataPP$numdate]
ShanDivPP[ ,depth:=MetadataPP$depth]
ggplot(ShanDivPP[depth == "p1-10" | depth == "1" | depth == "2" | depth == "3" |
                   depth == "4" | depth == "5" | depth == "10" ], aes(x=V1, y=depth, color=factor(numdate)))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75),aes(group=numdate))+
  scale_y_discrete(limits=rev(c("1", "2", "3", "4", "5", "10", "p1-10")))+
  scale_color_manual(values=c("#00BFC4", "#F8766D"), labels=c("Pre-Fire", "Post-Fire"))+
  xlab("Shannon Diversity Index")+
  ylab("Depth (cm)")+
  theme_light()+
  theme(axis.text = element_text(size = 12))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_blank())+
  ggtitle("Phillip's Samples - 321west")

dat4test <- ShanDivPP[depth == "p1-3" | depth == "p1-10"]
t.test(V1 ~ depth, data=dat4test[numdate == "18309"]) #numdate=18309 or Date == "17-Feb-20" ...p=0.5872
t.test(V1 ~ depth, data=dat4test[numdate == "18177"])  #numdate=18177 or Date == "8-Oct-19" ...p=0.7613
t.test(V1 ~ numdate, data=dat4test[depth == "p1-3"])  #p=0.2152
t.test(V1 ~ numdate, data=dat4test[depth == "p1-10"]) #p=0.3074
unique(ShanDivPP$depth)

## ANOVAs
dat4test <- ShanDivPP[depth == "p1-10" | depth == "1" | depth == "2" | depth == "3" |
                        depth == "4" | depth == "5" | depth == "10" ]
anova <- aov(V1 ~ depth, data=dat4test[numdate == "18309"])
summary(anova)  #p=0.948
anova <- aov(V1 ~ numdate, data=dat4test)
summary(anova)  #p=0.0736

dat4test <- ShanDivPP[depth == "p1-10" | depth == "1" | depth == "2" | depth == "3" |
                        depth == "4" | depth == "5" | depth == "10" ]
anova <- aov(V1 ~ depth, data=dat4test[numdate == "18177"])
summary(anova)  #p=0.995

dat4test <- ShanDivPP[depth == "p1-3" | depth == "1" | depth == "2" | depth == "3" ]
anova <- aov(V1 ~ depth, data=dat4test[numdate == "18309"])
summary(anova)  #p=0.611
anova <- aov(V1 ~ numdate, data=dat4test)
summary(anova)  #p=0.187

ShanDivPP <- ShanDivPP[V1 != 0]


dat4test <- ShanDivPP[depth == "p1-3" | depth == "1" | depth == "2" | depth == "3"]
anova <- aov(V1 ~ depth, data=dat4test[numdate == "18177"])
summary(anova)  #p=0.742

#
## COMP240!!!
#

#subset for the two dates that are pre- and post- fire for Phillip's samples
unique(MetadataWithReads$Date)
MetadataPP <- MetadataWithReads[Date == "2/17/20" | Date == "10/8/19", ]
MetadataPP <- MetadataPP[Who == "Phillip" | Who == "PhillipPool"]
MetadataPP <- MetadataPP[Plot == "240"]

# Subset OTUtax table for only the column names that match the SeqIDs in MetadataPP
OTUtaxPP <- otu.df.rownames[otu.df.rownames$SeqID %in% MetadataPP$SeqID, ]
OTUtaxPP[1:3, 1:5]
rownames(OTUtaxPP) <- c() #clear rownames
OTUtaxPP <- column_to_rownames(OTUtaxPP, var="SeqID") #move SeID column into rownames
unique(ShanDivPP$depth)

#Shannon Diversity
ShanDivPP <- as.data.table(diversity(OTUtaxPP, index="shannon"))
ShanDivPP[ ,numdate:=MetadataPP$numdate]
ShanDivPP[ ,depth:=MetadataPP$depth]
ggplot(ShanDivPP[depth == "p1-10" | depth == "1" | depth == "2" | depth == "3" |
                   depth == "4" | depth == "5" | depth == "10" ], aes(x=V1, y=depth, color=factor(numdate)))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75),aes(group=numdate))+
  scale_y_discrete(limits=rev(c("1", "2", "3", "4", "5", "10", "p1-10")))+
  scale_color_manual(values=c("#00BFC4", "#F8766D"), labels=c("Pre-Fire", "Post-Fire"))+
  xlab("Shannon Diversity Index")+
  ylab("Depth (cm)")+
  theme_light()+
  theme(axis.text = element_text(size = 12))+
  theme(legend.text = element_text(size = 12))+
  theme(legend.title = element_blank())+
  ggtitle("Phillip's Samples - 240")
# PDF export: 6x5 (1-10), 5x4 (1-3), 4.5x4 (pools only)

dat4test <- ShanDivPP[depth == "p1-3" | depth == "p1-10"]
t.test(V1 ~ depth, data=dat4test[numdate == "18309"]) #numdate=18309 or Date == "17-Feb-20" ...p=0.4787
t.test(V1 ~ depth, data=dat4test[numdate == "18177"])  #numdate=18177 or Date == "8-Oct-19" ...p=0.7495
t.test(V1 ~ numdate, data=dat4test[depth == "p1-3"])  #p=0.2454
t.test(V1 ~ numdate, data=dat4test[depth == "p1-10"]) #p=0.7693
unique(ShanDivPP$depth)

## ANOVAs
dat4test <- ShanDivPP[depth == "p1-10" | depth == "1" | depth == "2" | depth == "3" |
                        depth == "4" | depth == "5" | depth == "10" ]
anova <- aov(V1 ~ depth, data=dat4test[numdate == "18309"])
summary(anova)  #p=0.104

dat4test <- ShanDivPP[depth == "p1-10" | depth == "1" | depth == "2" | depth == "3" |
                        depth == "4" | depth == "5" | depth == "10" ]
anova <- aov(V1 ~ depth, data=dat4test[numdate == "18177"])
summary(anova)  #p=0.562

dat4test <- ShanDivPP[depth == "p1-3" | depth == "1" | depth == "2" | depth == "3" ]
anova <- aov(V1 ~ depth, data=dat4test[numdate == "18309"])
summary(anova)  #p=0.673

dat4test <- ShanDivPP[depth == "p1-3" | depth == "1" | depth == "2" | depth == "3"]
anova <- aov(V1 ~ depth, data=dat4test[numdate == "18177"])
summary(anova)  #p=0.444

# Simpson Diversity
SimpDivPP <- as.data.table(diversity(OTUtaxPP, index="simpson"))
SimpDivPP[ ,Date:=MetadataPP$Date]
SimpDivPP[ ,depth:=MetadataPP$depth]
ggplot(SimpDivPP, aes(x=V1, y=depth, color=Date))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75),aes(group=Date))+
  scale_y_discrete(limits=rev(c("1", "2", "3", "p1-3", "4", "5", "10", "p1-10", "20")))+
  xlab("Simpson Diversity Index")+
  ylab("Depth (cm)")+
  theme_light()+
  theme(axis.text = element_text(size = 12))+
  theme(legend.text = element_text(size = 12))+
  ggtitle("Phillip's Samples - 321west")

dat4test <- SimpDivPP[depth == "p1-3" | depth == "p1-10"]
t.test(V1 ~ depth, data=dat4test[Date == "17-Feb-20"])
t.test(V1 ~ depth, data=dat4test[Date == "8-Oct-19"])
t.test(V1 ~ Date, data=dat4test[depth == "p1-3"])
t.test(V1 ~ Date, data=dat4test[depth == "p1-10"])


#calculate species RICHNESS
# Where ASV abundance value is greater than 0, give it a 1, then sum the 1's by sample:
OTUtaxPP.rich <- apply(OTUtaxPP > 0,1,sum)
OTUtaxPP.rich.dt <- as.data.table(OTUtaxPP.rich)
OTUtaxPP.rich.dt[ ,depth:=MetadataPP$depth]
OTUtaxPP.rich.dt[ ,Date:=MetadataPP$Date]
names(OTUtaxPP.rich.dt)[1] <- "Richness"
ggplot(OTUtaxPP.rich.dt, aes(x=Richness, y=depth, color=Date))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75),aes(group=Date))+
  scale_y_discrete(limits=rev(c("1", "2", "3", "p1-3", "4", "5", "10", "p1-10", "20")))+
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+
  xlab("Richness")+
  ylab("Depth (cm)")+
  theme_light()+
  theme(axis.text = element_text(size = 12))+
  theme(legend.text = element_text(size = 12))+
  ggtitle("Phillip's Samples - 321west")

dat4test <- OTUtaxPP.t.rich.dt[depth == "p1-3" | depth == "p1-10"]
t.test(Richness ~ depth, data=dat4test[Date == "17-Feb-20"])
t.test(Richness ~ depth, data=dat4test[Date == "8-Oct-19"])
t.test(Richness ~ Date, data=dat4test[depth == "p1-3"])
t.test(Richness ~ Date, data=dat4test[depth == "p1-10"])

dat4test <- OTUtaxPP.t.rich.dt[Date == "17-Feb-20"]
anova <- aov(Richness ~ depth, data=dat4test)
summary(anova) #no sig difs
tukey <- TukeyHSD(anova, "depth")
result <- data.frame(tukey$depth)
result["p.adj"]

dat4test <- OTUtaxPP.t.rich.dt[Date == "8-Oct-19"]
anova <- aov(Richness ~ depth, data=dat4test)
summary(anova)#no sig difs#no sig difs
tukey <- TukeyHSD(anova, "depth")
result <- data.frame(tukey$depth)
result["p.adj"]

#calculate species EVENNESS
OTUtaxPP.even <- as.data.table(diversity(OTUtaxPP, index="simpson"))/log(OTUtaxPP.rich)
OTUtaxPP.even[ ,Date:=MetadataPP$Date]
OTUtaxPP.even[ ,depth:=MetadataPP$depth]
names(OTUtaxPP.even)[1] <- "Evenness"
ggplot(OTUtaxPP.even, aes(x=Evenness, y=depth, color=Date))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75),aes(group=Date))+
  scale_y_discrete(limits=rev(c("1", "2", "3", "p1-3", "4", "5", "10", "p1-10", "20")))+
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+
  xlab("Evenness")+
  ylab("Depth (cm)")+
  theme_light()+
  theme(axis.text = element_text(size = 12))+
  theme(legend.text = element_text(size = 12))+
  ggtitle("Phillip's Samples - 321west")

dat4test <- OTUtaxPP.t.even[depth == "p1-3" | depth == "p1-10"]
t.test(Evenness ~ depth, data=dat4test[Date == "17-Feb-20"])
t.test(Evenness ~ depth, data=dat4test[Date == "8-Oct-19"])
t.test(Evenness ~ Date, data=dat4test[depth == "p1-3"])
t.test(Evenness ~ Date, data=dat4test[depth == "p1-10"])
##

#
##
#
##
#
##
#

### DIVERSTIY METRICS on main "Neemonika" data
#
# Calculate diversity metric with the vegan package
ShanDivAll <- as.data.table(diversity(otu.df, index="shannon")) #index options = shannon, simpson, or invsimpson
head(ShanDivAll)
ShanDivAll$SeqID <- rownames(otu.df)
ShanDivAll <- merge(ShanDivAll, MetaNeemonika[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(ShanDivAll)
View(ShanDivAll)
ShanDivAll$Plot <- factor(ShanDivAll$Plot, levels=c("240", "321west", "321east", "400"))
ShanDivAll$Fire <- factor(ShanDivAll$Fire, levels=c("Pre-fire", "Post-fire"))
library(ggplot2)
ggplot(ShanDivAll, aes(x=numdate, y=V1, color=Plot))+
  geom_point(aes(shape=Fire), size=1)+
  geom_smooth(aes(fill=Plot),  size=1, alpha=0.1, method="loess")+
  theme_classic()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black", size=14))+
  theme(axis.text.y = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  scale_x_continuous(breaks=seq(17815, 18325, 30), limits=c(17815, 18325))+
  scale_y_continuous(breaks=seq(2, 5, 1), limits=c(2, 5))+
  scale_shape_manual(values=c(5, 19))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  xlab("Time")+
  ylab("Shannon Diversity")+
  ggtitle("All Neemonika Shannon Diversity (3-6cm samples omitted)")
###export as .png or .tiff: 600x350


#calculate species RICHNESS
# Where ASV abundance value is greater than 0, give it a 1, then sum the 1's by sample:
otu.df.rich <- apply(otu.df > 0, 1, sum)
head(otu.df.rich)
otu.df.rich.dt <- as.data.table(otu.df.rich)
otu.df.rich.dt[1:5,]
dim(ot.df.rich.dt)
head(MetaNeemonika)
otu.df.rich.dt[ ,SeqID:=rownames(otu.df)]
otu.df.rich.dt <- merge(otu.df.rich.dt, MetaNeemonika[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(otu.df.rich.dt)
names(otu.df.rich.dt)[2] <- "Richness"

otu.df.rich.dt$Plot <- factor(otu.df.rich.dt$Plot, levels=c("240", "321west", "321east", "400"))
otu.df.rich.dt$Fire <- factor(otu.df.rich.dt$Fire, levels=c("Pre-fire", "Post-fire"))

ggplot(otu.df.rich.dt, aes(x=numdate, y=Richness, color=Plot, fill=Plot))+
  geom_point(aes(shape=Fire), size=1)+
  geom_smooth(aes(fill=Plot),  size=1, alpha=0.1)+
  theme_classic()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black", size=14))+
  theme(axis.text.y = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  scale_x_continuous(breaks=seq(17815, 18325, 30), limits=c(17815, 18325))+
  scale_y_continuous(breaks=seq(0, 800, 100), limits=c(0, 800))+
  scale_shape_manual(values=c(5, 19))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  xlab("Time")+
  ylab("Richness")+
  ggtitle("All Neemonika Richness")
###export as .png or .tiff: 600x350, PDF: 6x4


#calculate species EVENNESS
otu.df.even <- as.data.table(diversity(otu.df, index="shannon"))/log(otu.df.rich)
otu.df.even[ ,SeqID:=rownames(otu.df)]
otu.df.even <- merge(otu.df.even, MetaNeemonika[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(otu.df.even)
names(otu.df.even)[2] <- "Evenness"

otu.df.even$Plot <- factor(otu.df.even$Plot, levels=c("240", "321west", "321east", "400"))
otu.df.even$Fire <- factor(otu.df.even$Fire, levels=c("Pre-fire", "Post-fire"))

ggplot(otu.df.even, aes(x=numdate, y=Evenness, color=Plot, fill=Plot))+
  geom_point(aes(shape=Fire), size=1)+
  geom_smooth(aes(fill=Plot),  size=1, alpha=0.1)+
  theme_classic()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black", size=14))+
  theme(axis.text.y = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  scale_x_continuous(breaks=seq(17815, 18325, 30), limits=c(17815, 18325))+
  #scale_y_continuous(breaks=seq(0.3, 0.9, 0.1), limits=c(0.3, 0.9))+
  scale_shape_manual(values=c(5, 19))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  xlab("Time")+
  ylab("Evenness")+
  ggtitle("All Neemonika Evenness (3-6cm samples omitted)")
###export as .png or .tiff: 600x350












######################################################
### HELLLINGER TRANSFORMATION + PCA + PERMANOVA  #####
######################################################
library(vegan)
library(phyloseq)
library(ggplot2)
library(data.table)
library(ggforce) #for geom_mark_hull() on PCA plots
library(concaveman) #for geom_mark_hull() on PCA plots
#
##
#
### Main ITS data - "Neemonika" samples
#
# subset data:
ps2.Neemonika <- prune_samples(sample_data(ps2)$Who == "Neemonika" & 
                                 sample_data(ps2)$Plot != "380" &
                                 sample_data(ps2)$numdate < 18309, ps2)

ps2.Neemonika #15026 taxa, 285 samples
# copy phyloseq object for Hellinger Transformation
ps2.Neemonika.hel <- ps2.Neemonika
# Hellinger Transformation
otu_table(ps2.Neemonika.hel) <-otu_table(decostand(otu_table(ps2.Neemonika.hel), method = "hellinger"), taxa_are_rows=FALSE)

# PCA
ps2.Neemonika.hel.pca <- rda(otu_table(ps2.Neemonika.hel))
#
# Broken Stick Model
screeplot(ps2.Neemonika.hel.pca, bstick = TRUE, 
          npcs = length(ps2.Neemonika.hel.pca$CA$eig),
          main = "Hellinger PCA - All Neemonika Data, 10cm and 0-3cm only") #1500x600
colnames(MetaNeemonika)
sample_data(ps2.Neemonika.hel)$Plot <- factor(sample_data(ps2.Neemonika.hel)$Plot, levels=c("240", "321west", "321east", "400"))
### PCA plot with elipses
library(ggplot2)
plot_ordination(ps2.Neemonika.hel, ps2.Neemonika.hel.pca, color="Fire")+
  stat_ellipse(geom="polygon", aes(fill=Fire), alpha=0.05)+
  #scale_color_gradient2(midpoint=40, low="green", mid="orange", high="blue")+
  #scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
  #                   labels=c("1c", "2c", "lo", "hi"))+
  #scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
  #                   labels=c("1c", "2c", "lo", "hi"))+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Hellinger PCA - All Neemonika") #600x500


#PERMANOVA to test for a statistically significant differences between Plot:
adonis(otu_table(ps2.Neemonika.hel) ~ sample_data(ps2.Neemonika)$Fire, permutations = 1000) 
#                                  Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# sample_data(ps2.Neemonika)$Fire   1     5.423  5.4225  24.443 0.0795 0.000999 ***
# Residuals                       283    62.782  0.2218         0.9205           
# Total                           284    68.204                 1.0000           
# p=0.001 for OTUs or ASVs


### PCA plot with gradient colors
plot_ordination(ps2.Neemonika.hel, ps2.Neemonika.hel.pca, color="numdate")+
  #stat_ellipse(geom="polygon", aes(fill=pH), alpha=0.2)+
  scale_color_gradient2(midpoint=18100, low="green", mid="orange", high="blue")+
  geom_point(size=2) +
  #theme_linedraw()+
  theme_minimal()+
  ggtitle("Hellinger PCA - All Neemonika Data") #600x500


#
##
#
##
#
##
#
### DOES DEPTH MATTER? Phillip's samples only
#
####### POST-FIRE comp240
ps2.F240pp <- prune_samples(sample_data(ps2)$Date== "2020-02-17" & 
                              sample_data(ps2)$Plot== "240" & 
                              sample_data(ps2)$Who == "Phillip" | 
                              sample_data(ps2)$Date== "2020-02-17" & 
                              sample_data(ps2)$Plot== "240" & 
                              sample_data(ps2)$Who == "PhillipPool", ps2)
ps2.F240pp #54 samples, 15026 taxa

# copy phyloseq object for Hellinger Transformation
ps2.F240pp.hel <- ps2.F240pp
# Hellinger Transformation
otu_table(ps2.F240pp.hel) <-otu_table(decostand(otu_table(ps2.F240pp.hel), method = "hellinger"), taxa_are_rows=FALSE)
otu_table(ps2.F240pp.hel)[1:5,1:5]
mean(otu_table(ps2.F240pp.hel))
sd(otu_table(ps2.F240pp.hel))
hist(ps2.F240pp.hel)

# PCA
ps2.F240pp.hel.pca <- rda(otu_table(ps2.F240pp.hel))
summary(ps2.F240pp.hel.pca)
# Broken Stick Model
screeplot(ps2.F240pp.hel.pca, bstick = TRUE, 
          npcs = length(ps2.F240pp.hel.pca$CA$eig),
          main = "Hellinger PCA - February (no burn control) comp240") #1500x600

#### PCA plot by depth
plot_ordination(ps2.F240pp, ps2.F240pp.hel.pca, color="depth")+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=depth), alpha=0.1)+
  theme_classic()+
  geom_point(size=2) + 
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  ggtitle("Hellinger PCA - February (no burn control) comp240") 
#TIF export: 600x500 or PDF export: 5x4

#PERMANOVA to test for a statistically significant differences between reps:
adonis(otu_table(ps2.F240pp.hel) ~ sample_data(ps2.F240pp)$depth) #p=0.289

#
#
######### POST-FIRE comp321
#
# subset data
unique(sample_data(ps2)$Date)
ps2.F321pp <- prune_samples(sample_data(ps2)$Date== "2020-02-17" & 
                              sample_data(ps2)$Plot== "321west" & 
                              sample_data(ps2)$Who == "Phillip" | 
                              sample_data(ps2)$Date== "2020-02-17" & 
                              sample_data(ps2)$Plot== "321west" & 
                              sample_data(ps2)$Who == "PhillipPool", ps2)

ps2.F321pp #53 samples, 7622 taxa
# copy phyloseq object for Hellinger Transformation
ps2.F321.hel <- ps2.F321pp
# Hellinger Transformation
otu_table(ps2.F321.hel) <-otu_table(decostand(otu_table(ps2.F321.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps2.F321.hel.pca <- rda(otu_table(ps2.F321.hel))
# Broken Stick Modle:
screeplot(ps2.F321.hel.pca, bstick = TRUE, 
          npcs = length(ps2.F321.hel.pca$CA$eig),
          main = "Hellinger PCA - February (post-burn) comp321") #1500x600

# PERMANOVA to test for a statistically significant differences between reps:
adonis(otu_table(ps2.F321.hel) ~ sample_data(ps2.F321pp)$Rep) 
# p=0.001
#
# PCA plot by depth:
plot_ordination(ps2.F321.hel, ps2.F321.hel.pca, color="depth")+
  geom_mark_hull(concavity = 5,expand=0,radius=0,alpha=0.1,aes(fill=depth))+
  geom_point(size=2) + 
  theme_light()+
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  ggtitle("Hellinger PCA - February (post-burn) comp321") 
#TIF export: 600x500 or PDF export: 5x4
# PERMANOVA to test for a statistically significant differences by depth:
adonis(otu_table(ps2.F321.hel) ~ sample_data(ps2.F321pp)$depth) #p=0.5

#
#
#
###### PRE-FIRE 240
#
# subset data
ps2.O240pp <- prune_samples(sample_data(ps2)$Date== "8-Oct-19" & 
                              sample_data(ps2)$Plot== "240" & 
                              sample_data(ps2)$Who == "Phillip" | 
                              sample_data(ps2)$Date== "8-Oct-19" & 
                              sample_data(ps2)$Plot== "240" & 
                              sample_data(ps2)$Who == "PhillipPool", ps2)
ps2.O240.sub <- subset_samples(ps2.O240pp, dada2_reads != "1") #remove samples with only one read
ps2.O240pp #52 samples, 15026 taxa
ps2.O240.sub #50 samples, 15026 taxa
sample_data(ps2.O240pp)
# copy phyloseq object for Hellinger Transformation
ps2.O240.hel <- ps2.O240pp
# Hellinger Transformation
otu_table(ps2.O240.hel) <-otu_table(decostand(otu_table(ps2.O240.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps2.O240.hel.pca <- rda(otu_table(ps2.O240.hel))
# Broken Stick Model
screeplot(ps2.O240.hel.pca, bstick = TRUE, 
          npcs = length(ps2.O240.hel.pca$CA$eig),
          main = "Hellinger PCA - October (no burn control) comp240") #1500x600

# PCA plot by depth
plot_ordination(ps2.O240.hel, ps2.O240.hel.pca, color="depth")+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=depth))+ 
  geom_point(size=2) + 
  theme_classic()+
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  ggtitle("Hellinger PCA - October (no burn control) comp240") 
#TIF export: 600x500 or PDF export: 5x4
#PERMANOVA to test for a statistically significant differences by depths:
adonis(otu_table(ps2.O240.hel) ~ sample_data(ps2.O240pp)$depth) 
# p-value=0.183

#
#
#
######### PRE-FIRE comp321
#
# subset data
ps2.O321pp <- prune_samples(sample_data(ps2)$Date== "8-Oct-19" & 
                              sample_data(ps2)$Plot== "321west" & 
                              sample_data(ps2)$Who == "Phillip" | 
                              sample_data(ps2)$Date== "8-Oct-19" & 
                              sample_data(ps2)$Plot== "321west" & 
                              sample_data(ps2)$Who == "PhillipPool", ps2)
ps2.O321pp #51 samples, 7622 taxa

# copy phyloseq object for Hellinger Transformation
ps2.O321.hel <- ps2.O321pp
# Hellinger Transformation
otu_table(ps2.O321.hel) <-otu_table(decostand(otu_table(ps2.O321.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps2.O321.hel.pca <- rda(otu_table(ps2.O321.hel))
# Broken Stick plot
screeplot(ps2.O321.hel.pca, bstick = TRUE, 
          npcs = length(ps2.O321.hel.pca$CA$eig),
          main = "Hellinger PCA - October (pre-burn) comp321") #1500x600
#
# PCA by depth
plot_ordination(ps2.O321.hel, ps2.O321.hel.pca, color="depth")+
  #geom_polygon(aes(fill=depth), alpha=0.2) + 
  geom_mark_hull(concavity = 5, expand=0, radius=0, aes(fill=depth))+
  geom_point(size=2) + 
  theme_classic()+
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  ggtitle("Hellinger PCA - October (pre-burn) comp321") 
#TIF export: 600x500 or PDF export: 5x4
#PERMANOVA to test for a statistically significant differences by depth:
adonis(otu_table(ps2.O321.hel) ~ sample_data(ps2.O321pp)$depth)
# p=0.114

##
#
##
#
#
##
#



########################################
#### Venn Diagrams of TITAN Results ####
########################################
library(data.table)
library(reshape2)
library(eulerr)
library(venn) #okay if this gives an error! ..."there is no package called ‘admisc’"...? still works!
#
# Function that will convert an input dataframe into a presence/absense table
vennfun <- function(x) { 
  x$id <- seq(1, nrow(x)) #add a column of numbers that is required by Reshape
  xm <- melt(x, id.vars="id", na.rm=TRUE) #melt table into two columns: variables and values
  xc <- dcast(xm, value~variable, fun.aggregate=length) #list presence/absence of each value for each variable (1 or 0)
  rownames(xc) <- xc$value #make the value column the rownames
  xc$value <- NULL #remove redundent value column
  xc[xc > 0] <- 1 #replace all values greater than 0 with "1"
  xc #output the new dataframe
}

# Input dataframe should look something like this:
### Variable1   Variable2   Variable3
### value2      value5      value2
### value7      value8      value8
### etc.        etc.        etc.
# For example, the variables are treatments and the values are genes  
# that are differentially expressed between conditions
#

#load table of taxa-of-interest in the format as above:
TaxList <- fread(file.choose(), na.strings=c("", "NA")) # fill empty cells with NA
head(TaxList)
colSums(!is.na(TaxList)) #count the number of non-NA values in each column

#filter each column for unique OTUs:
temp <- sapply(TaxList, unique)
TaxList.unique <- as.data.frame(sapply(temp, "length<-", max(lengths(temp)))) 
colSums(!is.na(TaxList.unique))

#subset for columns the ones we want to build a venn diagram with:
TaxList.unique.sub <- TaxList.unique[,c(1:2,9:11)]

# create presence/absence table
TaxLists4venn <- vennfun(TaxList.unique.sub)

# Plot Venn Diagram
TheVenn <- venn(TaxLists4venn)
plot(TheVenn, quantities = TRUE, labels=c("BurnSCNM", "ControlSCNM", "Module1", "Module2", "Module3"))


# output details for specific sections of the Venn Diagram:
VennOverlap <- rownames(TaxLists4venn[TaxLists4venn$TITANburned == 1 & 
                                    TaxLists4venn$TITANcontrol == 1  & 
                                    TaxLists4venn$SCNMcontrol == 1  &
                                    TaxLists4venn$SCNMburned == 1, ])

write.csv(VennOverlap, "~/VennOverlap.csv")


######################################
#### Concatenate ITS and 16S data ####
######################################
# 
# Combine ITS and 16S OTU tables!
#
# Rename samples in all datasets to SampleIDs: "BS1, BS2, etc."
# Then match SampleIDs to add 16S ASVs to ITS ASVs
#
# create SampID variable for ITS dataset and ust it to replace SeqID in the OTU table
head(MetadataWithReads)
#create SampID column
MetadataWithReads <- unite(MetadataWithReads, SampID, c(SampleID_Type, SampleID_number), sep="", remove=FALSE)
#merge metadata SeqID and SampID columns with OTU table
otuITS.SampIDs <- merge(MetadataWithReads[,c(1,6)], otu.df.rownames, by="SeqID") 
otuITS.SampIDs[1:5, 1:5]
otuITS.SampIDs[1:5, 15025:15028]
otuITS.SampIDs$SampID #only BS samples!
otuITS.SampIDs$SeqID <- NULL #remove SeqID column, we'll use SampID from here onward
dim(otuITS.SampIDs) #285 x 15027
dim(otu.df.rownames) #285 x 15027
dim(MetadataWithReads[Who=="Neemonika" & numdate<18308 & Plot != "380",]) #284 x 25

# create SampID variable for 16S dataset and ust it to replace SeqID in the OTU table
tail(SAMPtable16S)
OTUtable16S[1:5, 1:5]#check that SeqIDs = rownames

SAMPtable16S <- unite(SAMPtable16S, SampID, c(SampleID_Type, SampleID_number), sep="", remove=FALSE)
otu16S.SampIDs <- merge(data.frame(SeqID=rownames(SAMPtable16S), SampID=SAMPtable16S[,"SampID"]), OTUtable16S, by.x="SeqID", by.y="row.names")
otu16S.SampIDs[31:35, 1:5]
otu16S.SampIDs[1:5, 15025:15028]
otu16S.SampIDs$SampID #everything!
otu16S.SampIDs$SeqID <- NULL

#subset for only Neemonika samples prior to the 321W burn, match SampIDs:
otu16S.NMSampIDs <- otu16S.SampIDs[otu16S.SampIDs$SampID %in% MetadataWithReads[Who=="Neemonika" & numdate<18308 & Plot != "380",]$SampID, ]
otu16S.NMSampIDs$SampID
otuITS.NMSampIDs <- otuITS.SampIDs[otuITS.SampIDs$SampID %in% MetadataWithReads[Who=="Neemonika" & numdate<18308 & Plot != "380",]$SampID, ]
otuITS.NMSampIDs$SampID

# merge 16S and ITS data together!
ITS16Sotu <- merge(otuITS.NMSampIDs, otu16S.NMSampIDs, by="SampID")
write.csv(ITS16Sotu, "~/ITS16Sotu.csv")
# 285 Neemonika samples in total, but only 255 are shared between 16S and ITS datasets...
dim(ITS16Sotu) #255x110139
ITS16Sotu[1:5,1:5] #ASVs
ITS16Sotu[1:5,110130:110139] #OTUs

#
##
#
##
#

# Add FUNGuild info to OTU-taxonomy table:
# $ python FUNGuild.py taxa -otu ps2_OTUtax_forFunGuild.txt -format tsv -column taxonomy -classifier unite
# $ python FUNGuild.py guild -taxa ps2_OTUtax_forFunGuild.taxa.txt
TaxGuild <- fread("~/ps2_OTUtax_forFunGuild.taxa.txt")

# combine ITS and 16S taxonomy tables!
colnames(TAXtable16S) # 16S
colnames(TaxGuild) # ITS
# ensure that both tables are data.frame's and the first column lists ASVs
TAXtable16S.df <- data.frame(TAXtable16S)
TAXtable16S.df <- rownames_to_column(TAXtable16S.df, var="ASV")
head(TAXtable16S.df)

# make taxonomic names easier to deal with:
TAXtable16S.df$Domain <- gsub("d__", "", TAXtable16S.df$Domain) #remove the "d__" at the beginning of each taxa name
TAXtable16S.df$Phylum <- gsub("p__", "", TAXtable16S.df$Phylum) #remove the "p__" at the beginning of each taxa name
TAXtable16S.df$Class <- gsub("c__", "", TAXtable16S.df$Class) #remove the "c__" at the beginning of each taxa name
TAXtable16S.df$Order <- gsub("o__", "", TAXtable16S.df$Order) #remove the "o__" at the beginning of each taxa name
TAXtable16S.df$Family <- gsub("f__", "", TAXtable16S.df$Family) #remove the "f__" at the beginning of each taxa name
TAXtable16S.df$Genus <- gsub("g__", "", TAXtable16S.df$Genus) #remove the "g__" at the beginning of each taxa name
TAXtable16S.df$Species <- gsub("s__", "", TAXtable16S.df$Species) #remove the "s__" at the beginning of each taxa name

# Subset TAX tables for these columns: ASV, Domain, Phylum, Class, Order, Family, Genus, Species
TAXtable16S.sub <- TAXtable16S.df[,c(1,3:9)]
TAXtableITS.sub <- TaxGuild[,c(1,2,3,5,4,6,7,8)]
colnames(TAXtable16S.sub)
colnames(TAXtableITS.sub)
# Rename the second column so that it matches for both ITS and 16S tables:
names(TAXtableITS.sub)[2] <- "Domain"

# concatenate ITS and 16S TAX tables:
TAXtableITS16S <- rbind(TAXtableITS.sub, TAXtable16S.sub)


############################################
### sncm.ft() funciton from Burns, et al ###
############################################
## From Burns et al:
install.packages("minpack.lm")
install.packages("Hmisc")
install.packages("stats4")
library(minpack.lm)
library(Hmisc)
library(stats4)

#Adam Burns - 2/10/2015
#aburns2@uoregon.edu
#From Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development
#Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity when stats=FALSE. For use in R.
#spp: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
#pool: A community table for defining source community (optional; Default=NULL).
#taxon: A table listing the taxonomic calls for each otu, with OTU ids as row names and taxonomic classifications as columns.
#If stats=TRUE the function will return fitting statistics.
#If stats=FALSE the function will return a table of observed and predicted values for each otu.


sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  library(minpack.lm)
  library(Hmisc)
  library(stats4)
  
  options(warn=-1)
  
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m.mle <- mle(sncm.LL, method="Nelder-Mead", start=list(m=0.1, sigma=0.1), nobs=length(p))
  ### Changed the Optimization method to Nelder-Mead from the default
  ### ...I think the default is L-BFGS-B? "method = if(!useLim) "BFGS" else "L-BFGS-B""
  ### L-BFGS-B is a quasi-Newton method that allows box constraints...
  ### Nelder-Mead is the default method for the optim function (which is used by the mle function)
  ### Nelder-Mean is "robust and relatively slow"
  ### optim() has three method options: quasi-Newton, Nelder-Mead, or Conjugate Gradient
  ### ...unclear if there's a compelling reason to use or avoid any of these options? They all seem very similar...
  
  
  ##Calculate Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, N, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  
  ### START OF OPTIONAL BINOMIAL MODEL
  ### This model breaks with our Blodgett Data, which causes the entire script to fail...
  ### thus, we omit this part of the original script.
  #
  #bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  #aic.bino <- AIC(bino.mle, k=2)
  #bic.bino <- BIC(bino.mle)
  ##Goodness of fit for binomial model
  #bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
  #Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  #RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  #bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  #### END OF BINOMIAL MODEL!
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  if(stats==TRUE){
    fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), 
                           #binoLL=numeric(), Rsqr.bino=numeric(), RMSE.bino=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), ## OPTIONAL BINOMIAL OUTPUTS!
                           poisLL=numeric(), Rsqr=numeric(), Rsqr.pois=numeric(), 
                           RMSE=numeric(),  RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), 
                           AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), 
                           Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, 
                      #bino.mle@details$value, Rsqr.bino, RMSE.bino, aic.bino, bic.bino, ## OPTIONAL BINOMIAL OUTPUTS!
                      pois.mle@details$value, Rsqr, Rsqr.pois, RMSE, RMSE.pois, 
                      aic.fit, bic.fit, aic.pois, bic.pois, N, nrow(spp), length(p), d)
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], pois.pred, pois.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'pois.pred', 'pois.pred.lwr', 'pois.pred.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    return(B)
  }
}

########################################
#### SNCM Neutral Model (Figure 5B) ####
########################################
#
## RUN SNCM.FIT!!! ...run the code in the previous section to define the sncm.fit() function
#

## New FireResponsive Indicator table from previous lit/observations!
FireTaxa <- fread("~/FireResponsiveTaxa_MASTER_TABLE.csv")
FireTaxaGenera <- FireTaxa[is.na(FireTaxa$Genus) == FALSE] #remove rows where Genus = NA
head(FireTaxaGenera)
FireTaxaGenera.pos <- FireTaxaGenera[FireTaxaGenera$ResponseToFire == "positive", ]
TAXtableITS16S[110130:110138, 1:7]

# Add fire-responsive genera info for taxonomy table
FireTaxTable <- merge(TAXtableITS16S, FireTaxaGenera.pos[,c(4:6)], by="Genus", all.y=TRUE) #2754 taxa!
# USE THIS TABLE TO HIGHLIGHT FIRE-RESPONSIVE TAXA!
head(FireTaxTable)

#
##
#
##
#
# 
## RUN SCNM.FIT ON COMBINED ITS & 16S DATA!
#
# add Plot information to ITS16Sotu table:
colnames(MetadataWithReads)
dim(ITS16Sotu) # 255 x 110139 (110138 taxa + SampID column)
dim(TAXtableITS16S) # 110138 x 8

ITS16Sotu[1:5,1:5] #ASVs=Fungi
ITS16Sotu[1:5,110130:110139] #OTUs=Bacteria (technically ASVs)
dim(ITS16Sotu)
# ITS16Sotu table comes from generating the input table for FastSpar!
# Merge Plot column onto the OTU table by matching the SampID columns:
ITS16Sotu.sites <- merge(ITS16Sotu, MetadataWithReads[,c(5,8)], by="SampID")
ITS16Sotu.sites[1:5, 1:5]
#convert ASV/OTU columns to numeric (somehow they became integers...)
ITS16Sotu.sites.num <- data.frame(lapply(ITS16Sotu.sites[,2:110139], as.numeric))
# add back Plot column
ITS16Sotu.Plot <- cbind(ITS16Sotu.sites$Plot, ITS16Sotu.sites.num)
ITS16Sotu.Plot[1:5, 1:5]
names(ITS16Sotu.Plot)[1] <- "Plot"
unique(ITS16Sotu.Plot$Plot)
# Convert to a data.table to make future manipulations easier:
ITS16Sotu.Plot <- as.data.table(ITS16Sotu.Plot)

#
##
#
## comp400! = Hi
#
#
# subset for Hi plot data
ITS16Sotu400 <- ITS16Sotu.Plot[Plot == "400",] #67 samples
class(ITS16Sotu400$OTU25) #numeric
ITS16Sotu400[1:5, 1:5]
# remove plot column and convert to a matrix
ITS16Sotu400$Plot <- NULL 
ITS16Sotu400.m <- as.matrix(ITS16Sotu400)
ITS16Sotu400.m[1:5, 1:5]
# keep only columns that sum to a value greater than zero:
ITS16Sotu400.nonzero <- ITS16Sotu400.m[ ,colSums(ITS16Sotu400.m) > 0] 
ITS16Sotu400.nonzero[1:5, 1:5]

# Fit the Neutral Model
BS400fit <- sncm.fit(spp=ITS16Sotu400.nonzero, stats=FALSE) #outputs model fit values
BS400fit.stats <- sncm.fit(spp=ITS16Sotu400.nonzero, stats=TRUE) #outputs general stats like R-squared, m, and AIC!

# move rownames to a column
BS400fit <- tibble::rownames_to_column(BS400fit, "OTU_ID")
# add a "group" column that can be used ot assign colors during plotting
# "if freq > pred.upr color the dot red, otherwise, if freq < pred.lwr, color the point yellow, otherwise black"
BS400fit_grouped <- BS400fit %>%
  mutate(group = ifelse(freq > pred.upr, "Above", ifelse(freq<pred.lwr,"Below","Neutral")))
class(BS400fit_grouped$p)

# add additional info to Neutral Model output table:
BS400fit_grouped.tax <- merge(BS400fit_grouped, TAXtableITS16S, by.x="OTU_ID", by.y="ASV", all.x=TRUE) # add taxonomy
BS400fit_grouped.tax.sub <- BS400fit_grouped.tax[!(BS400fit_grouped.tax$Domain =="Eukaryota"), ] # remove 16S Eukaryotes
BS400fit_grouped.taxFUN <- merge(BS400fit_grouped.tax, TaxGuild, by.x="OTU_ID", by.y="ASV", all.x=TRUE) #add FUNGuilds
write.csv(BS400fit_grouped.taxFUN, "~/BS400fit_grouped.taxFUN.csv") #export table

# load table with presence/absence info from Venn Diagrams and TITAN data
Inds <- fread(file.choose()) #TITAN_SCNM_ALL_DATA_forVenn.csv
head(Inds)
head(FireTaxTable)

# merge SCNM data onto Inds to highlight TITAN indicator taxa
BS400fit_TITAN <- merge(Inds, BS400fit_grouped.tax.sub, by.x="TITAN_comp400", by.y="OTU_ID", all.x=TRUE)
#BS400fit_Pyrophiles <- merge(Inds, BS400fit_grouped.tax.sub, by.x="TITAN400PyrophileOTU", by.y="OTU_ID", all.x=TRUE) #ASVs from TITAN of known pyrophiles
BS400fit_FireResponsive <- merge(FireTaxTable[,c(2,10)], BS400fit_TITAN, by.x="ASV", by.y="TITAN_comp400") #all ASVs representing previoulsy documented fire-responsive fungi & bacteria
BS400fit_FireResponsiveBACTERIA <- BS400fit_FireResponsive[Domain == "Bacteria"]
BS400fit_FireResponsiveFUNGI <- BS400fit_FireResponsive[Domain == "Fungi"]

# Plot Neutral Model fit:
ggplot(BS400fit_grouped.tax.sub)+
  geom_point(aes(x = log(p), y = freq, color=group, shape=Domain), size=1, stroke=0.7)+
  #circle TITAN indicators:
  geom_point(data=BS400fit_TITAN, shape=1, stroke=0.5, size=1.5, aes(x = log(p), y = freq), color="black")+ 
  geom_point(data=BS400fit_FireResponsiveFUNGI, shape=1, stroke=1, size=1.5, aes(x = log(p), y = freq), color="#20f038")+
  geom_point(data=BS400fit_FireResponsiveBACTERIA, shape=1, stroke=1, size=1.5, aes(x = log(p), y = freq), color="#00ffff")+ 
  guides(color = guide_legend(override.aes = list(size=2)))+
  guides(shape = guide_legend(override.aes = list(size=2)))+
  geom_line(aes(x = log(p), y = freq.pred), color="grey30")+
  geom_line(aes(x = log(p), y = pred.lwr), color="grey30", linetype=2)+
  geom_line(aes(x = log(p), y = pred.upr), color="grey30", linetype=2)+
  theme_classic()+
  scale_color_manual(values= c("orange", "purple", "grey"))+
  scale_shape_manual(values=c(15,17,4))+
  ggtitle("Hi-Burn Neutral Model fit, all time-points, ITS & 16S")+
  xlab("log(Mean Relative Abundance)")+
  ylab("Frequency")+
  annotate(geom="text", x=min(log(BS400fit_grouped.tax.sub$p)), y=max(freq), 
           label=paste0("R-squared = ", formatC(BS400fit.stats$Rsqr, digits=4)), 
           size=4, vjust = "inward", hjust = "inward")+
  annotate(geom="text", x=min(log(BS400fit_grouped.tax.sub$p)), y=0.9, 
           label=paste0("m = ", formatC(BS400fit.stats$m, digits=2)), # "digits" are non-zero numbers... digits=2: 0.0031, digits=4: 0.003192
           size=4, vjust = "inward", hjust = "inward")
#export tiff: 600x400
#export PDF: 5x3

# 
##
#
## comp321E! = Lo
#
#
# Subset for Lo burn plot data:
ITS16Sotu321E <- ITS16Sotu.Plot[Plot == "321east"] #94 samples
# remove plot column and convert to a matrix:
ITS16Sotu321E$Plot <- NULL
ITS16Sotu321E[1:5, 1:5]
ITS16Sotu321E.m <- as.matrix(ITS16Sotu321E)
ITS16Sotu321E.m[1:5, 1:5]

# keep only columns that sum to a value greater than zero:
ITS16Sotu321E.nonzero <- ITS16Sotu321E.m[ ,colSums(ITS16Sotu321E.m) > 0]
ITS16Sotu321E.nonzero[1:5, 1:5]

# Fit the Neutral Model
BS321Efit <- sncm.fit(spp=ITS16Sotu321E.nonzero, stats=FALSE) #outputs model fit values
BS321Efit.stats <- sncm.fit(spp=ITS16Sotu321E.nonzero, stats=TRUE) #outputs general stats like R-squared, m, and AIC!

# move rownames to a column
BS321Efit <- tibble::rownames_to_column(BS321Efit, "OTU_ID")
# add a "group" column that can be used ot assign colors during plotting
# "if freq > pred.upr color the dot red, otherwise, if freq < pred.lwr, color the point yellow, otherwise black"
BS321Efit_grouped <- BS321Efit %>%
  mutate(group = ifelse(freq > pred.upr, "Above", ifelse(freq<pred.lwr,"Below","Neutral")))

# Add taxonomy and FUNGuilds to Neutral Model output table:
BS321Efit_grouped.tax <- merge(BS321Efit_grouped, TAXtableITS16S, by.x="OTU_ID", by.y="ASV", all.x=TRUE)
BS321Efit_grouped.tax.sub <- BS321Efit_grouped.tax[!(BS321Efit_grouped.tax$Domain =="Eukaryota"), ]

# Add TITAN indicator info and presence/absence info from Venn Diagrams:
Inds <- fread(file.choose()) #TITAN_SCNM_ALL_DATA_forVenn.csv
head(Inds)
# merge SCNM data onto Inds to highlight TITAN indicators:
BS321Efit_TITAN <- merge(Inds, BS321Efit_grouped.tax.sub, by.x="TITAN_comp321E", by.y="OTU_ID", all.x=TRUE)
BS321Efit_FireResponsive <- merge(FireTaxTable[,c(2,10)], BS321Efit_TITAN, by.x="ASV", by.y="TITAN_comp321E") #all ASVs representing previoulsy documented fire-responsive fungi & bacteria
BS321Efit_FireResponsiveBACTERIA <- BS321Efit_FireResponsive[Domain == "Bacteria"]
BS321Efit_FireResponsiveFUNGI <- BS321Efit_FireResponsive[Domain == "Fungi"]

# Plot Neutral Model fit:
ggplot(BS321Efit_grouped.tax.sub)+
  geom_point(aes(x = log(p), y = freq, color=group, shape=Domain), size=1, stroke=0.7)+
  #circle TITAN indicators:
  geom_point(data=BS321Efit_TITAN, shape=1, stroke=0.5, size=1.5, aes(x = log(p), y = freq), color="black")+
  geom_point(data=BS321Efit_FireResponsiveFUNGI, shape=1, stroke=1, size=1.5, aes(x = log(p), y = freq), color="#20f038")+
  geom_point(data=BS321Efit_FireResponsiveBACTERIA, shape=1, stroke=1, size=1.5, aes(x = log(p), y = freq), color="#00ffff")+ 
  guides(color = guide_legend(override.aes = list(size=2)))+
  guides(shape = guide_legend(override.aes = list(size=2)))+
  geom_line(aes(x = log(p), y = freq.pred), color="grey30")+
  geom_line(aes(x = log(p), y = pred.lwr), color="grey30", linetype=2)+
  geom_line(aes(x = log(p), y = pred.upr), color="grey30", linetype=2)+
  theme_classic()+
  scale_color_manual(values= c("orange", "purple", "grey"))+
  scale_shape_manual(values=c(15,17,4))+
  ggtitle("Lo-Burn Neutral Model fit, all time-points, ITS & 16S")+
  xlab("log(Mean Relative Abundance)")+
  ylab("Frequency")+
  annotate(geom="text", x=min(log(BS321Efit_grouped.tax.sub$p)), y=max(freq), 
           label=paste0("R-squared = ", formatC(BS321Efit.stats$Rsqr, digits=4)), 
           size=4, vjust = "inward", hjust = "inward")+
  annotate(geom="text", x=min(log(BS321Efit_grouped.tax.sub$p)), y=0.9, 
           label=paste0("m = ", formatC(BS321Efit.stats$m, digits=2)), 
           size=4, vjust = "inward", hjust = "inward")
#export tiff/png: 600x400
#export PDF: 5x3

#
##
#
## comp321W! = 2c
#
#
# subset for 2c samples only:
ITS16Sotu321W <- ITS16Sotu.Plot[Plot == "321west"] #46 samples
# Remove plot column and convert to matrix:
ITS16Sotu321W$Plot <- NULL
ITS16Sotu321W[1:5, 1:5]
ITS16Sotu321W.m <- as.matrix(ITS16Sotu321W)
ITS16Sotu321W.m[1:5, 1:5]

# keep only columns that sum to a value greater than zero:
ITS16Sotu321W.nonzero <- ITS16Sotu321W.m[ ,colSums(ITS16Sotu321W.m) > 0]
dim(ITS16Sotu321W.nonzero) #22430 taxa x 46 samples
dim(ITS16Sotu321W.m) #110138 taxa
ITS16Sotu321W.nonzero[1:5, 1:5]

# Fit the Neutral Model
BS321Wfit <- sncm.fit(spp=ITS16Sotu321W.nonzero, stats=FALSE) #outputs model fit values
BS321Wfit.stats <- sncm.fit(spp=ITS16Sotu321W.nonzero, stats=TRUE) #outputs general stats like R-squared, m, and AIC!

# move rownames to a column
BS321Wfit <- tibble::rownames_to_column(BS321Wfit, "OTU_ID")
# add a "group" column that can be used ot assign colors during plotting
# "if freq > pred.upr color the dot red, otherwise, if freq < pred.lwr, color the point yellow, otherwise black"
BS321Wfit_grouped <- BS321Wfit %>%
  mutate(group = ifelse(freq > pred.upr, "Above", ifelse(freq<pred.lwr,"Below","Neutral")))

# Add Taxonomy info to Neutral Model output table:
BS321Wfit_grouped.tax <- merge(BS321Wfit_grouped, TAXtableITS16S, by.x="OTU_ID", by.y="ASV", all.x=TRUE)
BS321Wfit_grouped.tax.sub <- BS321Wfit_grouped.tax[!(BS321Wfit_grouped.tax$Domain =="Eukaryota"), ] #omit 16S Eukaryota...

# Add TITAN indicator info and presence/absence info from venn diagrams:
Inds <- fread(file.choose()) #TITAN_SCNM_ALL_DATA_forVenn.csv
head(Inds)
# merge SCNM data onto Inds to highlight taxa ID's by TITAN
BS321Wfit_TITAN <- merge(Inds, BS321Wfit_grouped.tax.sub, by.x="TITAN_comp321W", by.y="OTU_ID", all.x=TRUE)
BS321Wfit_FireResponsive <- merge(FireTaxTable[,c(2,10)], BS321Wfit_TITAN, by.x="ASV", by.y="TITAN_comp321W") #all ASVs representing previoulsy documented fire-responsive fungi & bacteria
BS321Wfit_FireResponsiveBACTERIA <- BS321Wfit_FireResponsive[Domain == "Bacteria"]
BS321Wfit_FireResponsiveFUNGI <- BS321Wfit_FireResponsive[Domain == "Fungi"]

# plot Neutral Model fit:
ggplot(BS321Wfit_grouped.tax.sub)+
  geom_point(aes(x = log(p), y = freq, color=group, shape=Domain), size=1, stroke=0.7)+
  #circle TITAN indicators:
  geom_point(data=BS321Wfit_TITAN, shape=1, stroke=0.5, size=1.5, aes(x = log(p), y = freq), color="black")+ 
  geom_point(data=BS321Wfit_FireResponsiveFUNGI, shape=1, stroke=1, size=1.5, aes(x = log(p), y = freq), color="#20f038")+
  geom_point(data=BS321Wfit_FireResponsiveBACTERIA, shape=1, stroke=1, size=1.5, aes(x = log(p), y = freq), color="#00ffff")+ 
  guides(color = guide_legend(override.aes = list(size=2)))+
  guides(shape = guide_legend(override.aes = list(size=2)))+
  geom_line(aes(x = log(p), y = freq.pred), color="grey30")+
  geom_line(aes(x = log(p), y = pred.lwr), color="grey30", linetype=2)+
  geom_line(aes(x = log(p), y = pred.upr), color="grey30", linetype=2)+
  theme_classic()+
  scale_color_manual(values= c("orange", "purple", "grey"))+
  scale_shape_manual(values=c(15,17,4))+
  ggtitle("Control2 Neutral Model fit, all time-points, ITS & 16S")+
  xlab("log(Mean Relative Abundance)")+
  ylab("Frequency")+
  annotate(geom="text", x=min(log(BS321Wfit_grouped.tax.sub$p)), y=max(freq), 
           label=paste0("R-squared = ", formatC(BS321Wfit.stats$Rsqr, digits=4)), 
           size=4, vjust = "inward", hjust = "inward")+
  annotate(geom="text", x=min(log(BS321Wfit_grouped.tax.sub$p)), y=0.9, 
           label=paste0("m = ", formatC(BS321Wfit.stats$m, digits=2)), 
           size=4, vjust = "inward", hjust = "inward")
#export tiff: 600x400
#export PDF: 5x3


#
##
#
## comp240! = 1c
#
#
# Subset for 1c samples:
ITS16Sotu240 <- ITS16Sotu.Plot[Plot == "240"] #48 samples
# remove Plot column and convert to a matrix:
ITS16Sotu240$Plot <- NULL
ITS16Sotu240[1:5, 1:5]
ITS16Sotu240.m <- as.matrix(ITS16Sotu240)
ITS16Sotu240.m[1:5, 1:5]

# keep only columns that sum to a value greater than zero:
ITS16Sotu240.nonzero <- ITS16Sotu240.m[ ,colSums(ITS16Sotu240.m) > 0]
dim(ITS16Sotu240.nonzero) #22906 taxa x 48 samples
dim(ITS16Sotu240.m) #110138 taxa
ITS16Sotu240.nonzero[1:5, 1:5]

# Fit the Neutral Model
BS240fit <- sncm.fit(spp=ITS16Sotu240.nonzero, stats=FALSE) #outputs model fit values
BS240fit.stats <- sncm.fit(spp=ITS16Sotu240.nonzero, stats=TRUE) #outputs general stats like R-squared, m, and AIC!

# move rownames to a column
BS240fit <- tibble::rownames_to_column(BS240fit, "OTU_ID")
# add a "group" column that can be used ot assign colors during plotting
# "if freq > pred.upr color the dot red, otherwise, if freq < pred.lwr, color the point yellow, otherwise black"
BS240fit_grouped <- BS240fit %>%
  mutate(group = ifelse(freq > pred.upr, "Above", ifelse(freq<pred.lwr,"Below","Neutral")))

# Add taxonomy data to neutral model output table:
BS240fit_grouped.tax <- merge(BS240fit_grouped, TAXtableITS16S, by.x="OTU_ID", by.y="ASV", all.x=TRUE)
BS240fit_grouped.tax.sub <- BS240fit_grouped.tax[!(BS240fit_grouped.tax$Domain =="Eukaryota"), ] #omit 16S Eukaryota

# Add TITAN indicators and presence/absence info from Venn diagrams:
Inds <- fread(file.choose()) #TITAN_SCNM_ALL_DATA_forVenn.csv
head(Inds)

# merge SCNM data onto Inds to highlight taxa ID's by TITAN
BS240fit_TITAN <- merge(Inds, BS240fit_grouped.tax.sub, by.x="TITAN_comp240", by.y="OTU_ID", all.x=TRUE)
BS240fit_FireResponsive <- merge(FireTaxTable[,c(2,10)], BS240fit_TITAN, by.x="ASV", by.y="TITAN_comp240") #all ASVs representing previoulsy documented fire-responsive fungi & bacteria
BS240fit_FireResponsiveBACTERIA <- BS240fit_FireResponsive[Domain == "Bacteria"]
BS240fit_FireResponsiveFUNGI <- BS240fit_FireResponsive[Domain == "Fungi"]

# Plot Neutral Model fit:
ggplot(BS240fit_grouped.tax.sub)+
  geom_point(aes(x = log(p), y = freq, color=group, shape=Domain), size=1, stroke=0.7)+
  #circle TITAN indicators:
  geom_point(data=BS240fit_TITAN, shape=1, stroke=0.5, size=1.5, aes(x = log(p), y = freq), color="black")+ 
  geom_point(data=BS240fit_FireResponsiveFUNGI, shape=1, stroke=1, size=1.5, aes(x = log(p), y = freq), color="#20f038")+
  geom_point(data=BS240fit_FireResponsiveBACTERIA, shape=1, stroke=1, size=1.5, aes(x = log(p), y = freq), color="#00ffff")+ 
  guides(color = guide_legend(override.aes = list(size=2)))+
  guides(shape = guide_legend(override.aes = list(size=2)))+
  geom_line(aes(x = log(p), y = freq.pred), color="grey30")+
  geom_line(aes(x = log(p), y = pred.lwr), color="grey30", linetype=2)+
  geom_line(aes(x = log(p), y = pred.upr), color="grey30", linetype=2)+
  theme_classic()+
  scale_color_manual(values= c("orange", "purple", "grey"))+
  scale_shape_manual(values=c(15,17,4))+
  ggtitle("Control1 Neutral Model fit, all time-points, ITS & 16S")+
  xlab("log(Mean Relative Abundance)")+
  ylab("Frequency")+
  annotate(geom="text", x=min(log(BS240fit_grouped.tax.sub$p)), y=max(freq), 
           label=paste0("R-squared = ", formatC(BS240fit.stats$Rsqr, digits=4)), 
           size=4, vjust = "inward", hjust = "inward")+
  annotate(geom="text", x=min(log(BS240fit_grouped.tax.sub$p)), y=0.9, 
           label=paste0("m = ", formatC(BS240fit.stats$m, digits=2)), 
           size=4, vjust = "inward", hjust = "inward")
#export tiff: 600x400
#export PDF: 5x3


######################################################
##### Stacked Bar Plots of SCNCM fits (Figure 5A) ####
######################################################
#
## Basic Stacked Barplots summarizing the content of each SCNM graph:
#
# use the table() function to count how many of each varibale is in the group column (above, below, or neutral)
BS400fit.sub <- as.data.frame(table(BS400fit_grouped.tax.sub$group)) 
BS400fit.sub$Plot <- "hi"
BS400fit.sub.prop <- mutate_if(BS400fit.sub, is.numeric, funs(./sum(.)*100)) #convert counts to proportions by dividing each value by the column sum

BS240fit.sub <- as.data.frame(table(BS240fit_grouped.tax.sub$group)) 
BS240fit.sub$Plot <- "1c"
BS240fit.sub.prop <- mutate_if(BS240fit.sub, is.numeric, funs(./sum(.)*100)) 

BS321Wfit.sub <- as.data.frame(table(BS321Wfit_grouped.tax.sub$group)) 
BS321Wfit.sub$Plot <- "2c"
BS321Wfit.sub.prop <- mutate_if(BS321Wfit.sub, is.numeric, funs(./sum(.)*100)) 

BS321Efit.sub <- as.data.frame(table(BS321Efit_grouped.tax.sub$group)) 
BS321Efit.sub$Plot <- "lo"
BS321Efit.sub.prop <- mutate_if(BS321Efit.sub, is.numeric, funs(./sum(.)*100)) 

fitProp <- as.data.frame(rbind(BS321Efit.sub.prop, BS321Wfit.sub.prop, BS240fit.sub.prop, BS400fit.sub.prop))
fitProp$Plot <- factor(fitProp$Plot, levels=c("1c", "2c", "lo", "hi"))
fitCount <- as.data.frame(rbind(BS321Efit.sub, BS321Wfit.sub, BS240fit.sub, BS400fit.sub))
fitCount$Plot <- factor(fitCount$Plot, levels=c("1c", "2c", "lo", "hi"))

ggplot(fitCount)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text.y = element_text(color="black", size=11))+
  theme(axis.text.x = element_text(color="black", size=11))+
  theme(legend.text = element_text(size=11))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("All Data")+
  ylab("Count (number of taxa)")+
  labs(fill="Group")

ggplot(fitProp)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=30))+
  theme(axis.title = element_text(color="black", size=30))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("All Data")+
  ylab("Frequency (%)")+
  xlab("")+
  labs(fill="Group")
#PDF = 5x5


#
##
#
##
#

## Barplots showing where previously described Fire Responsive Taxa fall on teh Neutral Model:
#
# use the table() function to count how many of each varibale is in the group column (above, below, or neutral)
BS400fit_TITAN.sub <- as.data.frame(table(BS400fit_TITAN$group)) 
BS400fit_TITAN.sub$Plot <- "hi"
BS400fit_TITAN.sub.prop <- mutate_if(BS400fit_TITAN.sub, is.numeric, funs(./sum(.)*100)) 

BS240fit_TITAN.sub <- as.data.frame(table(BS240fit_TITAN$group)) 
BS240fit_TITAN.sub$Plot <- "1c"
BS240fit_TITAN.sub.prop <- mutate_if(BS240fit_TITAN.sub, is.numeric, funs(./sum(.)*100)) 

BS321Wfit_TITAN.sub <- as.data.frame(table(BS321Wfit_TITAN$group)) 
BS321Wfit_TITAN.sub$Plot <- "2c"
BS321Wfit_TITAN.sub.prop <- mutate_if(BS321Wfit_TITAN.sub, is.numeric, funs(./sum(.)*100)) 

BS321Efit_TITAN.sub <- as.data.frame(table(BS321Efit_TITAN$group)) 
BS321Efit_TITAN.sub$Plot <- "lo"
BS321Efit_TITAN.sub.prop <- mutate_if(BS321Efit_TITAN.sub, is.numeric, funs(./sum(.)*100)) 

fitTITANprop <- as.data.frame(rbind(BS321Efit_TITAN.sub.prop,BS321Wfit_TITAN.sub.prop,BS240fit_TITAN.sub.prop,BS400fit_TITAN.sub.prop))
fitTITANprop$Plot <- factor(fitTITANprop$Plot, levels=c("1c", "2c", "lo", "hi"))
fitTITAN <- as.data.frame(rbind(BS400fit_TITAN.sub,BS240fit_TITAN.sub,BS321Wfit_TITAN.sub,BS321Efit_TITAN.sub))
fitTITAN$Plot <- factor(fitTITAN$Plot, levels=c("1c", "2c", "lo", "hi"))

ggplot(fitTITAN)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text.y = element_text(color="black", size=11))+
  theme(axis.text.x = element_text(color="black", size=11))+
  theme(legend.text = element_text(size=11))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("Fire Responsive Taxa")+
  ylab("Count (number of taxa)")+
  labs(fill="Group")

ggplot(fitTITANprop)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=30))+
  theme(axis.title = element_text(color="black", size=30))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("TITAN")+
  ylab("Frequency (%)")+
  labs(fill="Group")

#
##
#
##
#

## Barplots showing the connection between Network Modules and SCNM!
#
# use the table() function to count how many of each varibale is in the group column (above, below, or neutral)
ModuleTaxa <- fread(file.choose()) #ModuleTaxa.csv
head(ModuleTaxa)

#
## MODULE 1 !!
#
BS400fit_Mod1 <- merge(ModuleTaxa, BS400fit_grouped.tax.sub, by.x="Module1", by.y="OTU_ID")
BS240fit_Mod1 <- merge(ModuleTaxa, BS240fit_grouped.tax.sub, by.x="Module1", by.y="OTU_ID")
BS321Wfit_Mod1 <- merge(ModuleTaxa, BS321Wfit_grouped.tax.sub, by.x="Module1", by.y="OTU_ID")
BS321Efit_Mod1 <- merge(ModuleTaxa, BS321Efit_grouped.tax.sub, by.x="Module1", by.y="OTU_ID")

BS400fit_Mod1.sub <- as.data.frame(table(BS400fit_Mod1$group)) 
BS400fit_Mod1.sub$Plot <- "hi"
BS400fit_Mod1.sub.prop <- mutate_if(BS400fit_Mod1.sub, is.numeric, funs(./sum(.)*100)) 

BS240fit_Mod1.sub <- as.data.frame(table(BS240fit_Mod1$group)) 
BS240fit_Mod1.sub$Plot <- "1c"
BS240fit_Mod1.sub.prop <- mutate_if(BS240fit_Mod1.sub, is.numeric, funs(./sum(.)*100)) 

BS321Wfit_Mod1.sub <- as.data.frame(table(BS321Wfit_Mod1$group)) 
BS321Wfit_Mod1.sub$Plot <- "2c"
BS321Wfit_Mod1.sub.prop <- mutate_if(BS321Wfit_Mod1.sub, is.numeric, funs(./sum(.)*100)) 

BS321Efit_Mod1.sub <- as.data.frame(table(BS321Efit_Mod1$group)) 
BS321Efit_Mod1.sub$Plot <- "lo"
BS321Efit_Mod1.sub.prop <- mutate_if(BS321Efit_Mod1.sub, is.numeric, funs(./sum(.)*100)) 

fitMOD1prop <- as.data.frame(rbind(BS321Efit_Mod1.sub.prop, BS321Wfit_Mod1.sub.prop, BS240fit_Mod1.sub.prop, BS400fit_Mod1.sub.prop))
fitMOD1prop$Plot <- factor(fitMOD1prop$Plot, levels=c("1c", "2c", "lo", "hi"))
fitMOD1 <- as.data.frame(rbind(BS321Efit_Mod1.sub, BS321Wfit_Mod1.sub, BS240fit_Mod1.sub, BS400fit_Mod1.sub))
fitMOD1$Plot <- factor(fitMOD1$Plot, levels=c("1c", "2c", "lo", "hi"))

ggplot(fitMOD1)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text.y = element_text(color="black", size=11))+
  theme(axis.text.x = element_text(color="black", size=11))+
  theme(legend.text = element_text(size=11))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("Network Module 1")+
  ylab("Count (number of taxa)")+
  labs(fill="Group")

ggplot(fitMOD1prop)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=30))+
  theme(axis.title = element_text(color="black", size=30))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("Network Module 1")+
  ylab("Frequency (%)")+
  labs(fill="Group")

#
## MODULE 2 !!
#
BS400fit_Mod2 <- merge(ModuleTaxa, BS400fit_grouped.tax.sub, by.x="Module2", by.y="OTU_ID")
BS240fit_Mod2 <- merge(ModuleTaxa, BS240fit_grouped.tax.sub, by.x="Module2", by.y="OTU_ID")
BS321Wfit_Mod2 <- merge(ModuleTaxa, BS321Wfit_grouped.tax.sub, by.x="Module2", by.y="OTU_ID")
BS321Efit_Mod2 <- merge(ModuleTaxa, BS321Efit_grouped.tax.sub, by.x="Module2", by.y="OTU_ID")

BS400fit_Mod2.sub <- as.data.frame(table(BS400fit_Mod2$group)) 
BS400fit_Mod2.sub$Plot <- "hi"
BS400fit_Mod2.sub.prop <- mutate_if(BS400fit_Mod2.sub, is.numeric, funs(./sum(.)*100)) 

BS240fit_Mod2.sub <- as.data.frame(table(BS240fit_Mod2$group)) 
BS240fit_Mod2.sub$Plot <- "1c"
BS240fit_Mod2.sub.prop <- mutate_if(BS240fit_Mod2.sub, is.numeric, funs(./sum(.)*100)) 

BS321Wfit_Mod2.sub <- as.data.frame(table(BS321Wfit_Mod2$group)) 
BS321Wfit_Mod2.sub$Plot <- "2c"
BS321Wfit_Mod2.sub.prop <- mutate_if(BS321Wfit_Mod2.sub, is.numeric, funs(./sum(.)*100)) 

BS321Efit_Mod2.sub <- as.data.frame(table(BS321Efit_Mod2$group)) 
BS321Efit_Mod2.sub$Plot <- "lo"
BS321Efit_Mod2.sub.prop <- mutate_if(BS321Efit_Mod2.sub, is.numeric, funs(./sum(.)*100)) 

fitMOD2prop <- as.data.frame(rbind(BS321Efit_Mod2.sub.prop, BS321Wfit_Mod2.sub.prop, BS240fit_Mod2.sub.prop, BS400fit_Mod2.sub.prop))
fitMOD2prop$Plot <- factor(fitMOD2prop$Plot, levels=c("1c", "2c", "lo", "hi"))
fitMOD2 <- as.data.frame(rbind(BS321Efit_Mod2.sub, BS321Wfit_Mod2.sub, BS240fit_Mod2.sub, BS400fit_Mod2.sub))
fitMOD2$Plot <- factor(fitMOD2$Plot, levels=c("1c", "2c", "lo", "hi"))

ggplot(fitMOD2)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text.y = element_text(color="black", size=11))+
  theme(axis.text.x = element_text(color="black", size=11))+
  theme(legend.text = element_text(size=11))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("Network Module 1")+
  ylab("Count (number of taxa)")+
  labs(fill="Group")

ggplot(fitMOD2prop)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=30))+
  theme(axis.title = element_text(color="black", size=30))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("Network Module 2")+
  ylab("Frequency (%)")+
  labs(fill="Group")

#
## MODULE 3 !!
#
BS400fit_Mod3 <- merge(ModuleTaxa, BS400fit_grouped.tax.sub, by.x="Module3", by.y="OTU_ID")
BS240fit_Mod3 <- merge(ModuleTaxa, BS240fit_grouped.tax.sub, by.x="Module3", by.y="OTU_ID")
BS321Wfit_Mod3 <- merge(ModuleTaxa, BS321Wfit_grouped.tax.sub, by.x="Module3", by.y="OTU_ID")
BS321Efit_Mod3 <- merge(ModuleTaxa, BS321Efit_grouped.tax.sub, by.x="Module3", by.y="OTU_ID")

BS400fit_Mod3.sub <- as.data.frame(table(BS400fit_Mod3$group)) 
BS400fit_Mod3.sub$Plot <- "hi"
BS400fit_Mod3.sub.prop <- mutate_if(BS400fit_Mod3.sub, is.numeric, funs(./sum(.)*100)) 

BS240fit_Mod3.sub <- as.data.frame(table(BS240fit_Mod3$group)) 
BS240fit_Mod3.sub$Plot <- "1c"
BS240fit_Mod3.sub.prop <- mutate_if(BS240fit_Mod3.sub, is.numeric, funs(./sum(.)*100)) 

BS321Wfit_Mod3.sub <- as.data.frame(table(BS321Wfit_Mod3$group)) 
BS321Wfit_Mod3.sub$Plot <- "2c"
BS321Wfit_Mod3.sub.prop <- mutate_if(BS321Wfit_Mod3.sub, is.numeric, funs(./sum(.)*100)) 

BS321Efit_Mod3.sub <- as.data.frame(table(BS321Efit_Mod3$group)) 
BS321Efit_Mod3.sub$Plot <- "lo"
BS321Efit_Mod3.sub.prop <- mutate_if(BS321Efit_Mod3.sub, is.numeric, funs(./sum(.)*100)) 

fitMOD3prop <- as.data.frame(rbind(BS321Efit_Mod3.sub.prop, BS321Wfit_Mod3.sub.prop, BS240fit_Mod3.sub.prop, BS400fit_Mod3.sub.prop))
fitMOD3prop$Plot <- factor(fitMOD3prop$Plot, levels=c("1c", "2c", "lo", "hi"))
fitMOD3 <- as.data.frame(rbind(BS321Efit_Mod3.sub, BS321Wfit_Mod3.sub, BS240fit_Mod3.sub, BS400fit_Mod3.sub))
fitMOD3$Plot <- factor(fitMOD3$Plot, levels=c("1c", "2c", "lo", "hi"))

ggplot(fitMOD3)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text.y = element_text(color="black", size=11))+
  theme(axis.text.x = element_text(color="black", size=11))+
  theme(legend.text = element_text(size=11))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("Network Module 3")+
  ylab("Count (number of taxa)")+
  labs(fill="Group")

ggplot(fitMOD3prop)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=30))+
  theme(axis.title = element_text(color="black", size=30))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("Network Module 3")+
  ylab("Frequency (%)")+
  labs(fill="Group")







#########################################################################
######  Subset tables to create FastSpar inputs, and run FastSpar #######
#########################################################################
# FastSpar = newer, faster SparCC
# https://academic.oup.com/bioinformatics/article/35/6/1064/5086389
# github.com/scwatts/FastSpar
#...recommended by Jonathan Friedman who was the first author on the original SparCC paper

## INSTALLATION (in terminal)
# $ conda install -c bioconda -c conda-forge fastspar
# $ conda install -y parallel

#
## Input table is a simple OTU table with ASVs as rows and sampleIDs as columns in a .tsv format
#
##
# transpose concatenated OTU table (exclude the SampID column which will confuse the t() function)
ITS16Sotu.t <- as.data.frame(t(ITS16Sotu[,2:ncol(ITS16Sotu)]))
ITS16Sotu.t[1:5, 1:5] #note that the column/SeqID names dissappeared, so we have to add back the SeqIDs:
colnames(ITS16Sotu.t) <- ITS16Sotu$SampID 
ITS16Sotu.t[1:5, 1:5]
# remove rows that sum to zero (meaning those ASVs are not present in *any* of these samples)
# said another way, "keep rows that sum to a number greater than 0":
ITS16Sotu.t <- ITS16Sotu.t[rowSums(ITS16Sotu.t) > 0, ] 
ITS16Sotu.t[1:5, 1:5]
#
# Keep only taxa that are present in >10 samples (>10 samples is the minimum required by FastSpar)
ITS16Sotu.t <- ITS16Sotu.t[rowSums(ITS16Sotu.t>0) > 10, ] 

#check that there are no NAs... if this returns a zero, then that means that zero columns have an NA in them
sum(colSums(is.na(ITS16Sotu.t))>0) 
# add a column at the beginning of the table that contains the ASV# info, 
# with the column name that FastSpar will be expecting "#OTU ID" (this is the BIOM format...)
ITS16Sotu.t <- cbind(`#OTU ID`= rownames(ITS16Sotu.t), ITS16Sotu.t)
# export .tsv to use with FastSpar
write.table(ITS16Sotu.t, file = "~/AllNeemonika_ITS16Sotu_FastSparINPUT.tsv", row.names=FALSE, sep="\t")
# in theory you should immediately be able to use this table with FastSpar... 
# however the formatting of the first column name as it's exported from R is somehow incompatible with FastSpar
## Try running FastSpar with this dataset,
# if it complains that the "#OTU ID" column doesn't exist, 
# the open the table in Excel and copy & paste "#OTU ID" from Terminal into the header in Excel and re-save
# (it wont look like anything changes, but somehow this makes a difference)
# Then try re-running, and it should work!

#
##
#
# RUN FASTSPAR (in terminal):
#
# $ fastspar --otu_table AllNeemonika_ITS16Sotu_FastSparINPUT.tsv --correlation AllNeemonika_ITS16S_median_correlation.tsv --covariance AllNeemonika_ITS16S_median_covariance.tsv
#
##
#





########################################################
#########  Wrangling the outputs of FastSpar!  #########
####### Histograms, Test for Modularity, Zi & Pi #######
########################################################
library(igraph)
# Three output tables from FastSpar: covariance, correlation, and correlation pvalues...
# Focus on Correlation and P-value tables!
#
# Output tables are in the format of ASVxASV... melt, merge, and attach taxonomy!
#
# load the median_correlation.tsv and pvalue.tsv tables output but FastSpar
FastSparCor <- fread("~/AllNeemonika_median_correlation.tsv")
FastSparP <- fread("~/AllNeemonika_pvalues.tsv") 

FastSparCor[1:5, 1:5]
FastSparP[1:5, 1:5]
# rename the first column to something that's easier to work with in R
names(FastSparCor)[1] <- "OTU_ID"
names(FastSparP)[1] <- "OTU_ID"

# FastSpar calculated all pairwise correlations, in all directions
# but since directionality doesn't matter for these data, we have duplicate correlation values, 
# for example: ASV258-ASV13 and ASV13-ASV258 each have unique correlation values that are exactly the same. 
# Remove all duplicate values!
#
# remove the "lower triangle" of duplicate correlation values
FastSparCor <- column_to_rownames(FastSparCor, var="OTU_ID")
FastSparCor[lower.tri(FastSparCor,diag=TRUE)] <- NA #fill the "lower triangle with NAs
FastSparCor <- rownames_to_column(FastSparCor, var="OTU_ID")

# remove the "lower triangle" of duplicate p-values
FastSparP <- column_to_rownames(FastSparP, var="OTU_ID")
FastSparP[lower.tri(FastSparP,diag=TRUE)] <- NA #fill the "lower triangle with NAs
FastSparP <- rownames_to_column(FastSparP, var="OTU_ID")

# melt Correlation table:
FastSparCor.m <- melt(FastSparCor)
# specify column names:
ColnamesList <- c("Taxon1", "Taxon2", "Correlation")
names(FastSparCor.m) <- ColnamesList
#remove rows with NA's (generated when removing the "lower triangle")
FastSparCor.m <- FastSparCor.m[rowSums(is.na(FastSparCor.m)) == 0, ]
#
# melt p-value table:
FastSparP.m <- melt(FastSparP)
# specify column names:
ColnamesList <- c("Taxon1", "Taxon2", "p_value")
names(FastSparP.m) <- ColnamesList
#remove rows with NA's (generated when removing the "lower triangle")
FastSparP.m <- FastSparP.m[rowSums(is.na(FastSparP.m)) == 0, ]

#
##
#
##
#
# Plot a HISTOGRAM of correlation values
#
FastSparCorP <- merge(FastSparP.m, FastSparCor.m, by=c("Taxon1", "Taxon2")) #takes a couple minutes...
head(FastSparCorP)

range(FastSparCorP$Correlation) 
#-0.8238 to 0.8628 
mean(FastSparCorP$Correlation) 
# 0.0001545147 
sd(FastSparCorP$Correlation) 
# 0.04662308 

## Histogram of Correlation Values:
## note that plotting *all* the data without filtering by p-value first takes several minutes...
ggplot(FastSparCorP, aes(x=Correlation))+
  geom_histogram(bins=500, fill="grey25")+
  theme_bw()+
  theme(axis.title = element_text(angle=0, color="black", size=12))+
  theme(axis.text = element_text(angle=0, color="black", size=12))+
  #scale_x_continuous(breaks=seq(-1, 1, 0.1), limits=c(-1, 1), expand=c(0,0))+
  #scale_y_continuous(breaks=seq(0, 900000, 50000), limits=c(0, 900000), expand=c(0,0))+
  ggtitle("Histogram of correlation values")+
  xlab("Correlation Value")+
  ylab("Count (number of taxon-taxon pairs)")

# Scatterplot
ggplot(FastSparCorP, aes(x=Correlation, y=p_value))+
  geom_point()+
  theme_bw()+
  theme(axis.title = element_text(angle=0, color="black", size=12))+
  theme(axis.text = element_text(angle=0, color="black", size=12))

#
###
#
##
#
# IS THERE A MODULAR SUB-STRUSTRUCTURE TO THE NETWORK?
#
# code adapted from Whitman et al, 2019
# https://github.com/TheaWhitman/WoodBuffalo/blob/master/Paper_Analyses_Figures/Mega-Network_Analysis.pub.ipynb
# 
# keep only rows with a p-value < 0.01
FastSparP.msub <- FastSparP.m[FastSparP.m$p_value < 0.01,] #takes a minute or two..
#keep only rows with a correlation value >0.1 
dim(FastSparP.msub) # 421041 x 3
# Merge correlation value table and p-value table
FastSparCorP <- merge(FastSparP.msub, FastSparCor.m, by=c("Taxon1", "Taxon2")) #takes a couple minutes...
dim(FastSparCorP) # 421041 x  4
#generate a non-repeating list of all taxa that I want to be nodes in my network:
#stack the first two columns (Taxon1 and Taxon2) on top of eachother (base R function similar to melt)
#then keep only the first column and remove any repeating values
FastSparCorP$Taxon2 <- as.character(FastSparCorP$Taxon2) #both Taxon columns should be characters...
nodes <- unique(stack(FastSparCorP[1:2])[1])
names(nodes)[1] <- "ASV"
head(nodes)
#5812 nodes (p<0.01, >10samples)
#
#add a column that gives each ASV an ID, which we'll use in the edge list...
nodes <- rowid_to_column(nodes, "id")

# Create "from" (i.e. Taxon1) and "to" (i.e. Taxon2) columns that contain the ID number for each ASV in my nodes list
edges <- FastSparCorP %>% 
  left_join(nodes, by = c("Taxon1" = "ASV")) %>% 
  rename(from = id)
edges <- edges %>% 
  left_join(nodes, by = c("Taxon2" = "ASV")) %>% 
  rename(to = id)
dim(edges)

library(igraph)
network_properties <- function(edge_list){
  N = graph_from_data_frame(edge_list,directed=FALSE)
  # Creates the network
  
  m = ecount(N)
  # Number of edges
  n = vcount(N)
  # Number of nodes
  k = 2*m/n
  # Average degree
  apl = mean_distance(N,directed=FALSE)
  # Average path length
  c = transitivity(N, type="global")
  cAve = transitivity(N, type = "average")
  # Clustering coefficient of whole graph - 
  # Transitivity measures the probability that the adjacent vertices of a vertex are connected.
  cl.mean = mean(closeness(N)) #note that this step fails when subsetting my data for 2SD or 4SD from the mean...
  cl.sd = sd(closeness(N)) #note that this step fails when subsetting my data for 2SD or 4SD from the mean...
  # closeness of graph
  # "how many steps are required to access every other vertex from a given vertex"
  ed = edge_density(N)
  # edge density of graph
  # ratio of the number of edges vs. the number of all possible edges
  d = diameter(N)
  # diameter of graph
  # length of the longest path between two nodes
  
  # Turn it into a dataframe
  df = data.frame(m,n,k,apl,c,cAve,cl.mean,cl.sd,ed,d)
  return(df)
}

# save summary table of network properties
network.properties_Allp0.01 <- network_properties(edges) 

# create an igraph object:
igraph <- graph_from_data_frame(edges, directed = FALSE) 
##
#
# Identify community structure (modules) within the correlation network:
Modules <- cluster_fast_greedy(graph=igraph, merges = TRUE, modularity = TRUE, membership = TRUE)
# fast greedy modularity optimization algorithm for finding community structure in very large networks (part of igraph)
# ref: Clauset, Newman, and Moore, 2004

# Explore the output:
head(membership(Modules)) #lists which module each taxa belongs to
sizes(Modules) #lists the number of taxa in each module
# 19 modules, most taxa are in the first three

# Modularity (Q) is an index measuring the extent to which a network is divided into modules, 
# Q > 0.3 is "a good indicator of significant community structure in a network" (Clauset 2004)
# Q > 0.4 is "exceptional" (Newman 2006)
modularity(Modules)
# Q = 0.3291829 

#
##
#
##
#

## FURTHER EXPLORING NETWORK TOPOLOGY

#
# Connectivity of each node can be determined based on its: 
# within-module connectivity (Zi), and
# among-module connectivity (Pi) 
# (Guimera & Amaral 2005: defined Zi, Pi, and seven different node topologies based on by Zi x Pi thresholds)

# Zi and Pi can be used to classify the nodes based on the topological roles they play in the network
# Node topologies are organized into four categories: 
# MODULE HUBS (highly connected nodes within modules, Zi > 2.5),
# NETWORK HUBS (highly connected nodes within entire network, Zi > 2.5 and Pi > 0.62),
# CONNECTORS (nodes that connect modules, Pi > 0.62), and
# PERIPHERALS (nodes connected in modules with few outside connections, Zi < 2.5 and Pi < 0.62). 
# Olesen et al. 2007 adapted Guimera & Amaral's method, 
# simplifying the interpretation to the four categories listed above.

# To calculate the Zi and Pi of each node:
# First, find out which module it is in
# Make a list of all the other nodes in that module
# Calculate the connectivity of that node to all those other nodes
# Do this for each node
# Then, Zi is calculated as:
# (number of links from a given node to other nodes in the module - the average number for nodes in this module)
# Divided by the standard deviation of this value for nodes in this module.

### Calculate Pi for all nodes
adding_Pi <- function(nodes,Modules,igraph){
  
  Pi=data.frame(Name=rep(nodes$ASV, dim(Modules[])),
                HomeModuleNumber=rep(0,length(nodes$ASV)),
                OtherModuleNumber=rep(0,length(nodes$ASV)),
                TotalCON=rep(0,length(nodes$ASV)),
                CON=rep(0,length(nodes$ASV)))
  # Establish empty dataframe
  
  for (i in 1:length(nodes$ASV)){
    node = paste(nodes$ASV[i])
    HomeModuleNumber = Position(function(x) node %in% x, Modules[], nomatch = 0)
    ModuleNumbers = 1:dim(Modules[])
    n = length(ModuleNumbers)
    TotalCON = as.numeric(lengths(adjacent_vertices(igraph,node)))
    lowend = (i-1)*n+1
    highend = n*i
    Pi$Name[lowend:highend]=nodes$ASV[i]
    Pi$HomeModuleNumber[lowend:highend]=HomeModuleNumber
    Pi$OtherModuleNumber[lowend:highend]=ModuleNumbers
    Pi$TotalCON[lowend:highend]=TotalCON
    if(HomeModuleNumber !=0){    
      for (j in ModuleNumbers){
        OtherModuleNumber = j
        NodesInOtherModule = Modules[[OtherModuleNumber]]
        modgraph = induced_subgraph(graph=igraph, v=c(node,NodesInOtherModule))
        CON=as.numeric(lengths(adjacent_vertices(modgraph,node)))
        Pi$CON[Pi$HomeModuleNumber==HomeModuleNumber & Pi$OtherModuleNumber==OtherModuleNumber & Pi$Name==node]=CON
      }
    }
  }
  Pi$kk2 = (Pi$CON/Pi$TotalCON)^2
  return(Pi)
}
Pi <- adding_Pi(nodes, Modules, igraph) 
# ALL DATA: start @ 2:43 ...end @ 3:31....
# only took about 20min with the 3-6cm omitted data!

head(Pi)

Pifinal <- Pi %>%
  group_by(Name, HomeModuleNumber, TotalCON) %>%
  summarise(Sum=sum(kk2)) %>%
  mutate(Pi=1-Sum)

head(Pifinal)

### Calculate Zi for all nodes from the Pi data
ZiNEW <- Pi %>% 
  filter(HomeModuleNumber==OtherModuleNumber) %>% 
  mutate(MeanCON=mean(CON), SdCON=sd(CON), Zi=((CON-MeanCON)/SdCON))
head(ZiNEW)


# Bringing module data together with Pi and Zi thresholds for defining hubs, connectors, etc.
Making_module_data = function(Pifinal, ZiNEW){
  Pthresh = 0.62
  Zthresh = 2.5
  ModuleData=data.frame(Name=Pifinal$Name, 
                        Module=Pifinal$HomeModuleNumber,
                        TotalCON=Pifinal$TotalCON,
                        ModuleCON=ZiNEW$MeanCON,
                        Pi=Pifinal$Pi,
                        Zi=ZiNEW$Zi)
  ModuleData$Class = ifelse(ModuleData$Zi>Zthresh & ModuleData$Pi>Pthresh,"Network Hub",
                            ifelse(ModuleData$Zi>Zthresh & ModuleData$Pi<Pthresh,"Module Hub",
                                   ifelse(ModuleData$Zi<Zthresh & ModuleData$Pi>Pthresh,"Connector", "Peripheral")))
  return(ModuleData)
}
ModuleData <- Making_module_data(Pifinal, ZiNEW)
head(ModuleData)
hist(ModuleData$TotalCON, breaks=100, xlab="Total Connectivity", 
     ylab="Number of Nodes", main=paste("Histogram of how many nodes each node is connected to"))

## Generate a table of taxa (a.k.a. nodes) with all the Network Topology info from above (i.e. Pi, Zi, etc.)
add_modInfo = function(nodes,ModuleData){
  nodes$Pi=c()
  for (i in 1:dim(nodes)[1]){
    ASV = nodes$ASV[i]
    x = ifelse(ASV %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==ASV,]$Pi!=(-Inf),ModuleData[ModuleData$Name==ASV,]$Pi, NA))
    nodes$Pi[i] = x
  }
  
  nodes$Zi=c()
  for (i in 1:dim(nodes)[1]){
    ASV = nodes$ASV[i]
    x = ifelse(ASV %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==ASV,]$Zi!=(-Inf),ModuleData[ModuleData$Name==ASV,]$Zi, NA))
    nodes$Zi[i] = x
  }
  
  nodes$NetworkRole=c()
  for (i in 1:dim(nodes)[1]){
    ASV = nodes$ASV[i]
    x = ifelse(ASV %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==ASV,]$Class!=(-Inf),ModuleData[ModuleData$Name==ASV,]$Class, NA))
    nodes$NetworkRole[i] = x
  }
  
  nodes$Module=c()
  for (i in 1:dim(nodes)[1]){
    ASV = nodes$ASV[i]
    x = ifelse(ASV %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==ASV,]$Module!=(-Inf),ModuleData[ModuleData$Name==ASV,]$Module, NA))
    nodes$Module[i] = x
  }
  
  
  nodes$NetworkRoleColour = ifelse(nodes$NetworkRole=="Connector","red",
                                   ifelse(nodes$NetworkRole=="Module Hub","navy","white"))
  
  return(nodes)
}
nodes.mod <- add_modInfo(nodes, ModuleData) #takes a minute or two
colnames(nodes.mod)


# plot Zi vs. Pi, and highlight the four categories defined by Olesen et al 2007, based on Guimera & Amaral 2005
ggplot(ModuleData)+
  geom_point(aes(x=Pi,y=Zi,color=Class))+
  theme_bw()+
  geom_hline(yintercept=2.5, color="black", linetype="dashed")+
  geom_vline(xintercept=0.62, color="black", linetype="dashed")+
  scale_color_manual(values = c("#ff9506", "#0670ff", "#ff0670", "#4fb900"))




##################################################
######  Wrangling the outputs of FastSpar!  ######
####### Add additional info to node table ########
######### create tables for CYTOSCAPE! ###########
##################################################
#
# EDGES TABLE:
edges <- copy(FastSparCorP) 
head(edges)
# add a column for designating positive vs. negative correlations:
edges$Sign <- ifelse(edges$Correlation<0, "neg", "pos")
# now that we have a Sign column, we can take the absolute value of all coreelation values:
edges$Correlation <- abs(edges$Correlation)
sum(is.na(edges)) #check to make sure there are no NAs (should sum to zero)
head(edges)

edges.sub <- edges[edges$Correlation > 0.3, ]
dim(edges.sub) #13889 for >0.3 (omit 407152 edges)
dim(edges) #421041

#export Edge table
write.csv(edges.sub, "~/edges_ALLDATA_ForCytoscape.csv")

#double check that the taxa in my edge taxa and node table match:
EdgeTaxa <- unique(c(edges$Taxon1, edges$Taxon2))
length(EdgeTaxa) #5812 (all) or 1274 (>0.3)
head(EdgeTaxa)

head(nodes.mod) #from Zi & Pi section
NodeTaxa <- unique(nodes.mod$String)
length(NodeTaxa) #5812 or 4165

setequal(EdgeTaxa, NodeTaxa) #TRUE! 

#
##
#

# NODE TABLE --PART 1--
# ...copy & pasted from previous section...
# keep only rows with a p-value < 0.01
FastSparP.msub <- FastSparP.m[FastSparP.m$p_value < 0.01,] #takes a minute or two..
dim(FastSparP.msub) # 351807 x 3
# Merge correlation value table and p-value table
FastSparCorP <- merge(FastSparP.msub, FastSparCor.m, by=c("Taxon1", "Taxon2")) #takes a couple minutes...
dim(FastSparCorP)
#generate a non-repeating list of all taxa that I want to be nodes in my network:
#stack the first two columns (WAxon1 and Taxon2) on top of eachother (base R function similar to melt)
#then keep only the first column and remove any repeating values
FastSparCorP$Taxon2 <- as.character(FastSparCorP$Taxon2) #both Taxon columns should be characters...
nodes <- unique(stack(FastSparCorP[1:2])[1])
names(nodes)[1] <- "ASV"
head(nodes)

# NODES TABLE --Part 2--
# Copied from previous section... nodes object was used to calculate Modularity, Zi, Pi, etc..
# Add columns to this nodes.mod table that you'll want to use for network visualization!
colnames(nodes.mod)
head(nodes.mod)

# add taxonomy
head(TAXtableITS16S)
tail(nodes)
nodes.mod.tax <- merge(nodes.mod.OTU, TAXtableITS16S, by.x="OTU", by.y="ASV", all.x=TRUE)
head(nodes.mod.tax)
names(nodes.mod.tax)[1] <- "ASV"

# add Neutral Model Results:
nodes.mod.tax.SCNM <- merge(nodes.mod.tax, BS400fit_grouped[,c(1,10)], by.x="ASV", by.y="OTU_ID", all.x=TRUE)
nodes.mod.tax.SCNM <- merge(nodes.mod.tax.SCNM, BS240fit_grouped[,c(1,10)], by.x="ASV", by.y="OTU_ID", all.x=TRUE)
nodes.mod.tax.SCNM <- merge(nodes.mod.tax.SCNM, BS321Efit_grouped[,c(1,10)], by.x="ASV", by.y="OTU_ID", all.x=TRUE)
nodes.mod.tax.SCNM <- merge(nodes.mod.tax.SCNM, BS321Wfit_grouped[,c(1,10)], by.x="ASV", by.y="OTU_ID", all.x=TRUE)
colnames(nodes.mod.tax.SCNM)
names(nodes.mod.tax.SCNM)[13] <- "comp400"
names(nodes.mod.tax.SCNM)[14] <- "comp240"
names(nodes.mod.tax.SCNM)[15] <- "comp321E"
names(nodes.mod.tax.SCNM)[16] <- "comp321W"
dim(nodes.mod.tax.SCNM)

# add TITAN results: 
head(TITANresults)
nodes.mod.tax.SCNM.TITAN <- merge(nodes.mod.tax.SCNM, TITANresults[ ,c(1:2)], by.x="ASV", by.y="TITAN400_OTUid", all.x=TRUE)
nodes.mod.tax.SCNM.TITAN <- merge(nodes.mod.tax.SCNM.TITAN, TITANresults[ ,c(7:8)], by.x="ASV", by.y="TITAN321E_OTUid", all.x=TRUE)
nodes.mod.tax.SCNM.TITAN <- merge(nodes.mod.tax.SCNM.TITAN, TITANresults[ ,c(5:6)], by.x="ASV", by.y="TITAN240_OTUid", all.x=TRUE)
nodes.mod.tax.SCNM.TITAN <- merge(nodes.mod.tax.SCNM.TITAN, TITANresults[ ,c(3:4)], by.x="ASV", by.y="TITAN321W_OTUid", all.x=TRUE)
dim(nodes.mod.tax.SCNM.TITAN)
head(nodes.mod.tax.SCNM.TITAN)

# add columns that are essentially the three parts of the Venn Diagram between Burn vs. No Burn
#
head(nodes.mod.tax.SCNM.TITAN)
head(Dat4Venn.pa)
scnmBURN <- Dat4Venn.pa[Dat4Venn.pa$BurnNonNeutral == 1  & Dat4Venn.pa$ControlNonNeutral == 0, ]  #2521
head(scnmBURN)
names(scnmBURN)[1] <- "ASV"
names(scnmBURN)[2] <- "scnmBURN"
scnmBURN[,3] <- NULL
sum(scnmBURN$scnmBURN) #2521
scnmOVERLAP <- Dat4Venn.pa[Dat4Venn.pa$BurnNonNeutral == 1  & Dat4Venn.pa$ControlNonNeutral == 1, ] #1076
head(scnmOVERLAP)
names(scnmOVERLAP)[1] <- "ASV"
names(scnmOVERLAP)[2] <- "scnmOVERLAP"
scnmOVERLAP[,3] <- NULL
sum(scnmOVERLAP$scnmOVERLAP) #1076
scnmNOBURN <- Dat4Venn.pa[Dat4Venn.pa$BurnNonNeutral == 0  & Dat4Venn.pa$ControlNonNeutral == 1, ] #1503
head(scnmNOBURN)
names(scnmNOBURN)[1] <- "ASV"
names(scnmNOBURN)[3] <- "scnmNOBURN"
scnmNOBURN[,2] <- NULL
sum(scnmNOBURN$scnmNOBURN) #1503

nodes.mod.tax.SCNMpa.TITAN <- merge(nodes.mod.tax.SCNM.TITAN, scnmBURN, by="ASV", all.x=TRUE)
nodes.mod.tax.SCNMpa.TITAN <- merge(nodes.mod.tax.SCNMpa.TITAN, scnmOVERLAP, by="ASV", all.x=TRUE)
nodes.mod.tax.SCNMpa.TITAN <- merge(nodes.mod.tax.SCNMpa.TITAN, scnmNOBURN, by="ASV", all.x=TRUE)
head(nodes.mod.tax.SCNMpa.TITAN)
dim(nodes.mod.tax.SCNMpa.TITAN)

head(FireTaxTable.sub)
nodes.mod.tax.SCNMpa.TITAN.Pyros <- merge(nodes.mod.tax.SCNMpa.TITAN, FireTaxTable.sub[,c(2,9)], by="ASV", all.x=TRUE)

#allTITANsites <- fread(file.choose())
head(allTITANsites)
dim(allTITANsites)
nodes.mod.tax.SCNMpa.allTITAN.Pyros <- merge(nodes.mod.tax.SCNMpa.TITAN.Pyros, allTITANsites, by="ASV", all.x=TRUE)
nodes.mod.tax.SCNMpa.allTITAN.Pyros.venn <- merge(nodes.mod.tax.SCNMpa.allTITAN.Pyros, FireResponsiveTaxa, by.x="ASV", by.y="OTU", all.x=TRUE)
dim(nodes.mod.tax.SCNMpa.allTITAN.Pyros.venn)

head(ITS16S.PAtreat)
nodes.mod.tax.SCNMpa.allTITAN.Pyros.venn.PAtreat[nodes.mod.tax.SCNMpa.allTITAN.Pyros.venn.PAtreat$ASV == "ASV663",]
names(ITS16S.PAtreat)[1] <- "ASV"
nodes.mod.tax.SCNMpa.allTITAN.Pyros.venn.PAtreat <- merge(nodes.mod.tax.SCNMpa.allTITAN.Pyros.venn, ITS16S.PAtreat, by="ASV", all.y=TRUE)


write.csv(nodes.mod.tax.SCNMpa.allTITAN.Pyros.venn.PAtreat, "/Users/monikafischer/Desktop/Nodes_ALLDATA_ForCytoscape_0.2sub.csv")


nodes.mod.tax.SCNMpa.TITAN.venn <- merge(nodes.mod.tax.SCNMpa.allTITAN.Pyros.venn.PAtreat, FireResponsiveTaxa, by.x="ASV", by.y="OTU", all.x=TRUE)
head(TITANscnm_VennMiddleBlob) #256 taxa that are in the center bits of the TITANvSCNMvCONTROLvBURN venn
nodes.mod.tax.SCNMpa.TITAN.venn2 <- merge(nodes.mod.tax.SCNMpa.TITAN.venn, TITANscnm_VennMiddleBlob, by="ASV", all.x=TRUE)
dim(nodes.mod.tax.SCNMpa.TITAN.venn)
write.csv(nodes.mod.tax.SCNMpa.TITAN.venn2, "/Users/monikafischer/Desktop/Nodes_ALLDATA_0.2sub_ForCytoscape_update.csv")

head(nodes.mod.tax.SCNMpa.TITAN.venn2)

### Add update TITAN results and results from making a VennDiagram of TITAN vs. SCNM vs. Burn vs. NoBurn
FireResponsiveTaxa1 # overlaping taxa that were identified by TITAN and fell outside the SCNM (excluding overlap with any unburned taxa)
FireResponsiveTaxa2 # Burned TITAN taxa only (excluding overlap with any unburned taxa)
FireResponsiveTaxa3 # Burned SCNM taxa only (excluding overlap with any unburned taxa)

#update this table:
head(nodes.mod.OTU.tax.SCNM)

# load TITAN results:
TITANresults <- fread(file.choose(), na.strings=c("", "NA"))
head(TITANresults)

nodes.mod.OTU.tax.SCNM.TITAN <- merge(nodes.mod.OTU.tax.SCNM, TITANresults[ ,c(1:2)], by.x="OTU", by.y="TITAN400_OTUid", all.x=TRUE)
nodes.mod.OTU.tax.SCNM.TITAN <- merge(nodes.mod.OTU.tax.SCNM.TITAN, TITANresults[ ,c(3:4)], by.x="OTU", by.y="TITAN321W_OTUid", all.x=TRUE)
nodes.mod.OTU.tax.SCNM.TITAN <- merge(nodes.mod.OTU.tax.SCNM.TITAN, TITANresults[ ,c(5:6)], by.x="OTU", by.y="TITAN240_OTUid", all.x=TRUE)
nodes.mod.OTU.tax.SCNM.TITAN <- merge(nodes.mod.OTU.tax.SCNM.TITAN, TITANresults[ ,c(7:8)], by.x="OTU", by.y="TITAN321E_OTUid", all.x=TRUE)

dim(nodes.mod.OTU.tax.SCNM.TITAN)
head(nodes.mod.OTU.tax.SCNM.TITAN)
view(nodes.mod.OTU.tax.SCNM.TITAN)


FireResponsiveTaxa1.df <- data.table(OTU = FireResponsiveTaxa1,
                                     VennSection = "TITANSCNM")
FireResponsiveTaxa2.df <- data.table(OTU = FireResponsiveTaxa2,
                                     VennSection = "TITANonly")
FireResponsiveTaxa3.df <- data.table(OTU = FireResponsiveTaxa3,
                                     VennSection = "SCNMonly")

FireResponsiveTaxa <- data.table(rbind(FireResponsiveTaxa1.df,FireResponsiveTaxa2.df,FireResponsiveTaxa3.df))
nodes.mod.tax.SCNMpa.TITAN.venn <- merge(nodes.mod.tax.SCNMpa.TITAN, FireResponsiveTaxa, by.x="ASV", by.y="OTU", all.x=TRUE)

head(TITANscnm_VennMiddleBlob) #256 taxa that are in the center bits of the TITANvSCNMvCONTROLvBURN venn
nodes.mod.tax.SCNMpa.TITAN.venn2 <- merge(nodes.mod.tax.SCNMpa.TITAN.venn, TITANscnm_VennMiddleBlob, by="ASV", all.x=TRUE)
dim(nodes.mod.tax.SCNMpa.TITAN.venn)
write.csv(nodes.mod.tax.SCNMpa.TITAN.venn2, "~/Nodes_ALLDATA_0.3sub_ForCytoscape_update.csv")

#
##
#
nodes.mod.OTU.tax.SCNM.TITAN.venn.FUNGUILD <- merge(nodes.mod.OTU.tax.SCNM.TITAN.venn, TaxGuild, by.x="OTU", by.y="ASV", all.x=TRUE)
write.csv(nodes.mod.OTU.tax.SCNM.TITAN.venn.FUNGUILD, "~/Nodes_ALLDATA_ForCytoscape_withFUNGUILD.csv")

##
#
##
#

#######################################################
#### Timeline Plot of Blodgett Temp & Precip Data! ####
#######################################################
library(data.table)
library(ggplot2)
BStime <- fread("~/TOTALBlodgettTempPrecipData_SeqDatesOnly.csv")
head(BStime)
# Create a new column that is a merge of Date and Hour columns:
BStime[ ,DateHour:=do.call(paste, c(.SD, sep="/")), .SDcols=c(1,5)]
head(BStime)
BStime <- na.omit(BStime)

BStime$DateHour <- as.Date(BStime$DateHour, "%m/%d/%y/%H")
head(BStime)

ggplot(data = BStime, mapping=aes(x=DateHour, y=Prec_in, group=1))+ 
  geom_bar(stat = "identity", color="blue", fill="blue", width = 0.5)+ 
  geom_line(mapping = aes(y = (TempC+30)/7), color="red", size=0.5)+ #transform Temp
  scale_x_date(date_breaks = "2 weeks", date_labels = "%e %b %Y")+
  scale_y_continuous(limits = c(0, 10), 
                     "Precipitation [inches]", 
                     sec.axis = sec_axis(~ .*7-30, name = "Temperature [°C]"))+ #undo the Temp transformation above
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color="grey90"),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.3), 
        axis.ticks = element_line(color="black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size=12, color="black"),
        axis.text.x = element_text(size=12, hjust=1, vjust=0.5, angle=90),
        axis.text.y = element_text(size=12, color="black"))+
  xlab("Date")

### Summarize temp & precip data by date
# mean temp/day
BStemp.meansd <- BStime[ ,list(mean(TempC), 
                               mean(TempC)+sd(TempC), 
                               mean(TempC)-sd(TempC)),
                         by=DateHour, ]
colnames(BStemp.meansd) <- c("Date", "MeanTemp", "TempTopSD", "TempBottomSD")
BStemp.meansd <- na.omit(BStemp.meansd)

ggplot(BStemp.meansd, aes(x=Date, y=MeanTemp))+
  geom_point(shape=16, size=1)+
  geom_errorbar(aes(ymax=TempTopSD, ymin=TempBottomSD), width=0.1, alpha=0.2)+
  WhiteThemes

# sum of all precip/day
BSprecip <- BStime[ ,list(sum(Prec_in)), by=DateHour, ]

colnames(BSprecip) <- c("Date", "Precip")

ggplot(BSprecip, aes(Date, Precip))+
  geom_col()
WhiteThemes

ggplot(BSprecip, aes(Date, Precip))+
  geom_bar(stat="identity")