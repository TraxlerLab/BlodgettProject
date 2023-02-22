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

####################################
######### DIVERSITY METRICS ########
############# Figure 2 #############
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












##########################################################
### Rarefaction & Diversity Metrics with Rarefied data ###
###################### ITS only ##########################
##########################################################
# plot rarefaction curves using a phyloseq function for ITS data:
rarecurve(otu_table(ps2.Neemonika), step=50, cex=0.5, label = FALSE)
#start @ 11:50... end ~1pm ish...
# Export PDF 20x25
#
# this curve is the procedure of "rarefaction"... 
# ...rarefying is the normalization procedure!
#
# Rarefaction Curve = number of ASVs as a function of number of reads (a.k.a. "sample size")
# ...one line for each sample showing how many ASVs are detected for each number of reads

# rarefy the data:
ps2.Neemonika.rarefied = rarefy_even_depth(ps2.Neemonika, 
                                           rngseed=1, # set.seed for random number generation
                                           sample.size=min(sample_sums(ps2.Neemonika)), #reads will be removed such that all samples have the same number of reads as the sample with the fewest reads
                                           replace=FALSE, #resample with or without replacement, FALSE = more memory intensive!, TRUE = instantaneous!
                                           verbose=TRUE,
                                           trimOTUs=TRUE) #remove OTUs that have zero reads left after rarefaction

#sampling without replacement took a few minutes and resulted in 6220 OTUS being removed because they were no longer present
#sample with replacement was basically instantaneous and resulted in 6379 OTUs being removed because they were no longer present

#sample depths pre- and post- rarefying:
range(sample_sums(ps2.Neemonika)) #13751-84409
sample_sums(ps2.Neemonika.rarefied) #all the same! 13751

# Extract OTU table:
otu.df.rare <- as.data.frame(otu_table(ps2.Neemonika.rarefied))
otu.df.rare[1:5,1:5]
dim(otu.df.rare)

unique(otu.df.rich.dt$numdate)
otu.df.rich.dt[numdate==18086]
#ANOVA
anova <- aov(Richness~Plot, data=otu.df.rich.dt[numdate==17898 | numdate==17904])
summary(anova)
tukey <- TukeyHSD(anova, "Plot")
result <- data.frame(tukey$Plot)
result["p.adj"]

# 17816 & 17820 & 17833 p=0.89 (lo pre- and post fire, 1c, 2c)
# 17855 p=0.593
# 17862 p=0.485
# 17868 p=0.021 (posthoc TukeyHSD --> lo vs. 1c p=0.0190, all others ns)
# 17878 p=0.372
# 17898 p=0.174 (hi pre-fire, lo post-fire, 1c, 2c)
# 17898 & 17904 p=0.00233 (lo post-fire, hi pre- and post-fire, 1c, 2c)
# 17988 p=0.113
# 18024 p=0.0409 (hi vs. 1c p=0.298, all others ns)
# 18057 p=0.218
# 18086 p=0.016 (hi vs. 2c p=0.01291, all others ns)
# 18114 p=0.516
# 18140 p=0.00732 (hi vs. lo p=0.00428, all others ns)
# 18177 p=7.01e-07 (hi vs. lo, hi vs. 1c, and hi vs. 2c are all less than 1e-5, the rest are ns)
# 18206 p=0.726
# 18239 p=0.0398 (all p>0.05 after TukeyHSD...)
# 18267 p=0.000315 (hi vs. all other plots p<0.02, the rest are ns)
# 18305 p=0.129

### DIVERSTIY METRICS CALCULATED ON ALL DATA AT ONCE
# Calculate diversity metric with the vegan package
ShanDivAll.rare <- as.data.table(diversity(otu.df.rare, index="shannon")) #index options = shannon, simpson, or invsimpson
head(ShanDivAll.rare)
ShanDivAll.rare$SeqID <- rownames(otu.df.rare)
ShanDivAll.rare <- merge(ShanDivAll.rare, MetaNeemonika[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(ShanDivAll.rare)
ShanDivAll.rare$Plot <- factor(ShanDivAll.rare$Plot, levels=c("240", "321west", "321east", "400"))
ShanDivAll.rare$Fire <- factor(ShanDivAll.rare$Fire, levels=c("Pre-fire", "Post-fire"))

ggplot(ShanDivAll.rare, aes(x=numdate, y=V1, color=Plot))+
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
  ggtitle("All Neemonika ITS Shannon Diversity RAREFIED")
###export as .png or .tiff: 600x350


#calculate species RICHNESS
# Where ASV abundance value is greater than 0, give it a 1, then sum the 1's by sample:
otu.df.rare.rich <- apply(otu.df.rare > 0, 1, sum)
otu.df.rare[1:5, 1:5]
head(otu.df.rare.rich)
otu.df.rare.rich.dt <- as.data.table(otu.df.rare.rich)
otu.df.rare.rich.dt[1:5,]
dim(otu.df.rare.rich.dt)
head(MetaNeemonika)
otu.df.rare.rich.dt[ ,SeqID:=rownames(otu.df.rare)]
otu.df.rare.rich.dt <- merge(otu.df.rare.rich.dt, MetaNeemonika[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(otu.df.rare.rich.dt)
names(otu.df.rare.rich.dt)[2] <- "Richness"

otu.df.rare.rich.dt$Plot <- factor(otu.df.rare.rich.dt$Plot, levels=c("240", "321west", "321east", "400"))
otu.df.rare.rich.dt$Fire <- factor(otu.df.rare.rich.dt$Fire, levels=c("Pre-fire", "Post-fire"))

ggplot(otu.df.rare.rich.dt, aes(x=numdate, y=Richness, color=Plot, fill=Plot))+
  geom_point(aes(shape=Fire), size=1)+
  #stat_summary(fun = "median", size = 10, geom = "point", shape = "_")+
  #stat_summary(aes(y = Richness, group = Plot), fun = median, geom="line", size = 0.5) +
  #stat_smooth_func(geom="text", method="lm", parse=TRUE) +
  #geom_smooth(aes(fill=Plot),  method="lm", size=1, alpha=0.1, span = 0.95)+
  geom_smooth(aes(fill=Plot),  method="loess", linewidth=1, alpha=0.1)+
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
  ggtitle("ITS Richness - RAREFIED")
###export as .png or .tiff: 600x350, PDF: 6x4


#calculate species EVENNESS
otu.df.rare.even <- as.data.table(diversity(otu.df.rare, index="shannon"))/log(otu.df.rare.rich)
otu.df.rare.even[ ,SeqID:=rownames(otu.df.rare)]
otu.df.rare.even <- merge(otu.df.rare.even, MetaNeemonika[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(otu.df.rare.even)
names(otu.df.rare.even)[2] <- "Evenness"

otu.df.rare.even$Plot <- factor(otu.df.rare.even$Plot, levels=c("240", "321west", "321east", "400"))
otu.df.rare.even$Fire <- factor(otu.df.rare.even$Fire, levels=c("Pre-fire", "Post-fire"))

ggplot(otu.df.rare.even, aes(x=numdate, y=Evenness, color=Plot, fill=Plot))+
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
  ggtitle("All Neemonika Evenness RAREFIED")
###export as .png or .tiff: 600x350


#
##
#
##
#
# copy phyloseq object for Hellinger Transformation
ps2.Neemonika.rare.hel <- ps2.Neemonika.rarefied
# Hellinger Transformation
otu_table(ps2.Neemonika.rare.hel)[1:10, 1:10]
otu_table(ps2.Neemonika.rare.hel) <-otu_table(decostand(otu_table(ps2.Neemonika.rare.hel), method = "hellinger"), taxa_are_rows=FALSE)

# PCA
ps2.Neemonika.rare.hel.pca <- rda(otu_table(ps2.Neemonika.rare.hel))

# Broken Stick Model
screeplot(ps2.Neemonika.rare.hel.pca, bstick = TRUE, 
          npcs = length(ps2.Neemonika.rare.hel.pca$CA$eig),
          main = "Hellinger PCA - All Neemonika Data - RAREFIED") #1500x600
colnames(MetaNeemonika)
sample_data(ps2.Neemonika.rare.hel)$Plot <- factor(sample_data(ps2.Neemonika.rare.hel)$Plot, levels=c("240", "321west", "321east", "400"))
### PCA plot with elipses
plot_ordination(ps2.Neemonika.rare.hel, ps2.Neemonika.rare.hel.pca, color="Plot")+
  stat_ellipse(geom="polygon", aes(fill=Plot), alpha=0.05)+
  #scale_color_gradient2(midpoint=10, low="green", mid="orange", high="blue")+
  scale_fill_manual(values=c("#009acd",   "#00cd9a","#ff7f01","#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01","#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  scale_shape_manual(values=c(16,1, 4, 8))+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Hellinger PCA - All Neemonika - RAREFIED") #600x500 TIFF or 6x4 PDF

#PERMANOVA to test for a statistically significant differences between Plot:
adonis(otu_table(ps2.Neemonika.rare.hel) ~ sample_data(ps2.Neemonika.rarefied)$Fire, permutations = 1000) 
#                                  Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# sample_data(ps2.Neemonika)$Fire   1     5.407  5.4068  24.243 0.0789 0.000999 ***
# Residuals                       283    63.117  0.2230         0.9211           
# Total                           284    68.524                 1.0000           
# p=0.001 for OTUs or ASVs
#
## ORIGINAL RESULT WITHOUT RAREFYING!::
#                                  Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# sample_data(ps2.Neemonika)$Fire   1     5.423  5.4225  24.443 0.0795 0.000999 ***
# Residuals                       283    62.782  0.2218         0.9205           
# Total                           284    68.204                 1.0000           
# p=0.001 for OTUs or ASVs





##########################################################
### Rarefaction & Diversity Metrics with Rarefied data ###
###################### 16S only ##########################
##########################################################

# 16S data:
ps16S
ps16S.Neemonika <- prune_samples(sample_data(ps16S)$Who == "Neemonika" &  #only Neemonika data (no Phillip or Controls)
                                   sample_data(ps16S)$Plot != "380" &       #exclude comp380 data
                                   sample_data(ps16S)$numdate < 18308,      #exclude dates on/after 321W burn
                                 ps16S)  
#
range(sample_sums(ps16S.Neemonika)) #318-438483
range(sample_sums(ps2.Neemonika)) #13751-84409
# samples with suspiciously low read counts (<1000):
# BS92.321ET.16S, BS306.321W.16S, BS261.400B.16S, BS21.321ET.16S, BS113.400B.16S
# remove samples with greater than 300000 reads(only for generating Rarefaction curve...)
ps16S.Neemonika = prune_samples(sample_sums(ps16S.Neemonika)<=300000, ps16S.Neemonika)
# remove samples with less than 10000 reads because most blanks and KC's has <10000reads
ps16S.Neemonika = prune_samples(sample_sums(ps16S.Neemonika)>=10000, ps16S.Neemonika)
range(sample_sums(ps16S.Neemonika)) #10215-135722 for rarecurve, 10215-438483 otherwise

#write.csv(otu_table(ps16S.Neemonika), "/Users/monikafischer/Desktop/OTUtable_16Sneemonika_SampleReadsGreaterThan1000.csv")
#write.csv(sample_data(ps16S.Neemonika), "/Users/monikafischer/Desktop/SAMPtable_16Sneemonika_SampleReadsGreaterThan1000.csv")
#write.csv(tax_table(ps16S.Neemonika), "/Users/monikafischer/Desktop/TAXtable_16Sneemonika_SampleReadsGreaterThan1000.csv")

# 16S rarefaction curves:
rarecurve(otu_table(ps16S.Neemonika), step=50, cex=0.5, label = FALSE)
# Export PDF 20x25


# rarefy the data:
ps16S.Neemonika.rarefied = rarefy_even_depth(ps16S.Neemonika, 
                                             rngseed=1, # set.seed for random number generation
                                             sample.size=min(sample_sums(ps16S.Neemonika)), #reads will be removed such that all samples have the same number of reads as the sample with the fewest reads
                                             replace=TRUE, #resample with or without replacement, FALSE = more memory intensive!, TRUE = instantaneous!
                                             verbose=TRUE,
                                             trimOTUs=TRUE) #remove OTUs that have zero reads left after rarefaction

#sample with replacement was basically instantaneous and resulted in 65644 OTUs being removed because they were no longer present

#sample depths pre- and post- rarefying:
range(sample_sums(ps16S.Neemonika)) #1097 - 438483
sample_sums(ps16S.Neemonika.rarefied) #all the same! 1097

# Extract OTU table:
otu16S.df.rare <- as.data.frame(otu_table(ps16S.Neemonika.rarefied))
otu16S.df.rare[1:5,1:5]
dim(otu16S.df.rare)

#calculate species RICHNESS with RAREFIED data
# Where ASV abundance value is greater than 0, give it a 1, then sum the 1's by sample:
otu16S.df.rare.rich <- apply(otu16S.df.rare > 0, 1, sum)
head(otu16S.df.rare.rich)
otu16S.df.rare.rich.dt <- as.data.table(otu16S.df.rare.rich)
otu16S.df.rare.rich.dt[1:5,]
dim(otu16S.df.rare.rich.dt)
colnames(SAMPtable16S) 
otu16S.df.rare.rich.dt[ ,SeqID:=rownames(otu16S.df.rare)]
otu16S.df.rare[1:5, 1:5]
#View(SAMPtable16S)
SAMPtable16S.rn <- rownames_to_column(SAMPtable16S, var = "SeqID")


otu16S.df.rare.rich.dt <- merge(otu16S.df.rare.rich.dt, SAMPtable16S.rn[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(otu16S.df.rare.rich.dt)
names(otu16S.df.rare.rich.dt)[2] <- "Richness"

otu16S.df.rare.rich.dt$Plot <- factor(otu16S.df.rare.rich.dt$Plot, levels=c("240", "321west", "321east", "400"))
otu16S.df.rare.rich.dt$Fire <- factor(otu16S.df.rare.rich.dt$Fire, levels=c("Pre-fire", "Post-fire"))


ggplot(otu16S.df.rare.rich.dt, aes(x=numdate, y=Richness, color=Plot, fill=Plot))+
  geom_point(aes(shape=Fire), size=1)+
  #stat_summary(fun = "median", size = 10, geom = "point", shape = "_")+
  #stat_summary(aes(y = Richness, group = Plot), fun = median, geom="line", size = 0.5) +
  #stat_smooth_func(geom="text", method="lm", parse=TRUE) +
  geom_smooth(aes(fill=Plot),  method="loess", size=1, alpha=0.1, span = 0.95)+
  #geom_smooth(aes(fill=Plot),  method="lm", size=1, alpha=0.1)+
  theme_classic()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black", size=14))+
  theme(axis.text.y = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  scale_x_continuous(breaks=seq(17815, 18325, 30), limits=c(17815, 18325))+
  #scale_y_continuous(breaks=seq(0, 350, 50), limits=c(0, 350))+
  scale_shape_manual(values=c(5, 19))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  xlab("Time")+
  ylab("Richness")+
  ggtitle("16S Richness - RAREFIED")
###export as .png or .tiff: 600x350, PDF: 6x4



#
##
#
##
#
#calculate species RICHNESS with original, un-rarefied data
# Extract OTU table:
otu16S.df <- as.data.frame(otu_table(ps16S.Neemonika))
otu16S.df[1:5,1:5]
dim(otu16S.df)

# Where ASV abundance value is greater than 0, give it a 1, then sum the 1's by sample:
otu16S.df.rich <- apply(otu16S.df > 0, 1, sum)
head(otu16S.df.rich)
otu16S.df.rich.dt <- as.data.table(otu16S.df.rich)
otu16S.df.rich.dt[1:5,]
range(otu16S.df.rich.dt$otu16S.df.rich)

dim(otu16S.df.rich.dt)
colnames(SAMPtable16S)
otu16S.df.rich.dt[ ,SeqID:=rownames(otu16S.df)]
otu16S.df[1:5, 1:5]
#View(SAMPtable16S)
SAMPtable16S.rn <- rownames_to_column(SAMPtable16S, var = "SeqID")

otu16S.df.rich.dt <- merge(otu16S.df.rich.dt, SAMPtable16S.rn[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(otu16S.df.rich.dt)
names(otu16S.df.rich.dt)[2] <- "Richness"
range(otu16S.df.rich.dt$Richness)

otu16S.df.rich.dt$Plot <- factor(otu16S.df.rich.dt$Plot, levels=c("240", "321west", "321east", "400"))
otu16S.df.rich.dt$Fire <- factor(otu16S.df.rich.dt$Fire, levels=c("Pre-fire", "Post-fire"))



ggplot(otu16S.df.rich.dt, aes(x=numdate, y=Richness, color=Plot, fill=Plot))+
  geom_point(aes(shape=Fire), size=1)+
  #stat_summary(fun = "median", size = 10, geom = "point", shape = "_")+
  #stat_summary(aes(y = Richness, group = Plot), fun = median, geom="line", size = 0.5) +
  #stat_smooth_func(geom="text", method="lm", parse=TRUE) +
  geom_smooth(aes(fill=Plot),  method="lm", size=1, alpha=0.1, span = 0.95)+
  stat_regline_equation(label.y = 350, aes(label = ..rr.label..))+ #library(ggpubr)
  #geom_smooth(aes(fill=Plot),  method="loess", size=1, alpha=0.1)+
  theme_classic()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black", size=14))+
  theme(axis.text.y = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  scale_x_continuous(breaks=seq(17815, 18325, 30), limits=c(17815, 18325))+
  scale_y_continuous(breaks=seq(0, 3500, 500), limits=c(0, 3500))+
  scale_shape_manual(values=c(5, 19))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  xlab("Time")+
  ylab("Richness")+
  ggtitle("16S Richness")
###export as .png or .tiff: 600x350, PDF: 6x4




#
##
#
##
#

### DIVERSTIY METRICS CALCULATED ON ALL DATA AT ONCE
# Calculate diversity metric with the vegan package
ShanDivAll16S.rare <- as.data.table(diversity(otu16S.df.rare, index="shannon")) #index options = shannon, simpson, or invsimpson
head(ShanDivAll16S.rare)
ShanDivAll16S.rare$SeqID <- rownames(otu16S.df.rare)
ShanDivAll16S.rare <- merge(ShanDivAll16S.rare, SAMPtable16S[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(ShanDivAll16S.rare)
ShanDivAll16S.rare$Plot <- factor(ShanDivAll16S.rare$Plot, levels=c("240", "321west", "321east", "400"))
ShanDivAll16S.rare$Fire <- factor(ShanDivAll16S.rare$Fire, levels=c("Pre-fire", "Post-fire"))

ggplot(ShanDivAll16S.rare, aes(x=numdate, y=V1, color=Plot))+
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
  ggtitle("All Neemonika 16S Shannon Diversity RAREFIED")
###export as .png or .tiff: 600x350



#calculate species EVENNESS
otu16S.df.rare.even <- as.data.table(diversity(otu16S.df.rare, index="shannon"))/log(otu16S.df.rare.rich)
otu16S.df.rare.even[ ,SeqID:=rownames(otu16S.df.rare)]
otu16S.df.rare.even <- merge(otu16S.df.rare.even, SAMPtable16S[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(otu16S.df.rare.even)
names(otu16S.df.rare.even)[2] <- "Evenness"

otu16S.df.rare.even$Plot <- factor(otu16S.df.rare.even$Plot, levels=c("240", "321west", "321east", "400"))
otu16S.df.rare.even$Fire <- factor(otu16S.df.rare.even$Fire, levels=c("Pre-fire", "Post-fire"))

ggplot(otu16S.df.rare.even, aes(x=numdate, y=Evenness, color=Plot, fill=Plot))+
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
  ggtitle("16S Evenness RAREFIED")
###export as .png or .tiff: 600x350


#
##
#
##
#
# copy phyloseq object for Hellinger Transformation
ps2.Neemonika.rare.hel <- ps2.Neemonika.rarefied
# Hellinger Transformation
otu_table(ps2.Neemonika.rare.hel)[1:10, 1:10]
otu_table(ps2.Neemonika.rare.hel) <-otu_table(decostand(otu_table(ps2.Neemonika.rare.hel), method = "hellinger"), taxa_are_rows=FALSE)

# PCA
ps2.Neemonika.rare.hel.pca <- rda(otu_table(ps2.Neemonika.rare.hel))

# Broken Stick Model
screeplot(ps2.Neemonika.rare.hel.pca, bstick = TRUE, 
          npcs = length(ps2.Neemonika.rare.hel.pca$CA$eig),
          main = "Hellinger PCA - All Neemonika Data - RAREFIED") #1500x600
colnames(MetaNeemonika)
sample_data(ps2.Neemonika.rare.hel)$Plot <- factor(sample_data(ps2.Neemonika.rare.hel)$Plot, levels=c("240", "321west", "321east", "400"))
### PCA plot with elipses
plot_ordination(ps2.Neemonika.rare.hel, ps2.Neemonika.rare.hel.pca, color="Plot")+
  stat_ellipse(geom="polygon", aes(fill=Plot), alpha=0.05)+
  #scale_color_gradient2(midpoint=10, low="green", mid="orange", high="blue")+
  scale_fill_manual(values=c("#009acd",   "#00cd9a","#ff7f01","#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01","#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  scale_shape_manual(values=c(16,1, 4, 8))+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Hellinger PCA - All Neemonika - RAREFIED") #600x500 TIFF or 6x4 PDF

#PERMANOVA to test for a statistically significant differences between Plot:
adonis(otu_table(ps2.Neemonika.rare.hel) ~ sample_data(ps2.Neemonika.rarefied)$Fire, permutations = 1000) 
#                                  Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# sample_data(ps2.Neemonika)$Fire   1     5.407  5.4068  24.243 0.0789 0.000999 ***
# Residuals                       283    63.117  0.2230         0.9211           
# Total                           284    68.524                 1.0000           
# p=0.001 for OTUs or ASVs
#
## ORIGINAL RESULT WITHOUT RAREFYING!::
#                                  Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# sample_data(ps2.Neemonika)$Fire   1     5.423  5.4225  24.443 0.0795 0.000999 ***
# Residuals                       283    62.782  0.2218         0.9205           
# Total                           284    68.204                 1.0000           
# p=0.001 for OTUs or ASVs




######################################################
### HELLLINGER TRANSFORMATION + PCA + PERMANOVA  #####
############### ITS Only - Figure 2 ################## 
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
### Main ITS data - "Neemonika" samples = core sample set for manuscript
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
          main = "Hellinger PCA - All Neemonika Data") #1500x600
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


#PERMANOVA to test for a statistically significant differences between Burned vs. No-Burn Controls:
adonis(otu_table(ps2.Neemonika.hel) ~ sample_data(ps2.Neemonika)$Fire, permutations = 1000) 
#                                  Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# sample_data(ps2.Neemonika)$Fire   1     5.423  5.4225  24.443 0.0795 0.000999 ***
# Residuals                       283    62.782  0.2218         0.9205           
# Total                           284    68.204                 1.0000           
# p=0.001 for OTUs or ASVs


### PCA plot with gradient colors
plot_ordination(ps2.Neemonika.hel, ps2.Neemonika.hel.pca, color="numdate")+
  scale_color_gradient2(midpoint=18100, low="green", mid="orange", high="blue")+
  geom_point(size=2) +
  theme_minimal()+
  ggtitle("Hellinger PCA - All Neemonika Data") #600x500


#
##
#
##
#
##
#
### DOES DEPTH MATTER? PCAs of Phillip's samples only
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



######################################################
### HELLLINGER TRANSFORMATION + PCA + PERMANOVA  #####
############### 16S Only - Figure 2 ################## 
######################################################
ps16S.Neemonika

# copy phyloseq object for Hellinger Transformation
ps16S.Neemonika.hel <- ps16S.Neemonika
# Hellinger Transformation
otu_table(ps16S.Neemonika.hel)[1:10, 1:10]
otu_table(ps16S.Neemonika.hel) <-otu_table(decostand(otu_table(ps16S.Neemonika.hel), method = "hellinger"), taxa_are_rows=FALSE)

# PCA
ps16S.Neemonika.hel.pca <- rda(otu_table(ps16S.Neemonika.hel))

# Broken Stick Model
screeplot(ps16S.Neemonika.hel.pca, bstick = TRUE, 
          npcs = length(ps16S.Neemonika.hel.pca$CA$eig),
          main = "Hellinger PCA - All Neemonika Data, 10cm and 0-3cm only") #1500x600

sample_data(ps16S.Neemonika.hel)$Plot <- factor(sample_data(ps16S.Neemonika.hel)$Plot, levels=c("240", "321west", "321east", "400"))
### PCA plot with elipses
library(ggplot2)
plot_ordination(ps16S.Neemonika.hel, ps16S.Neemonika.hel.pca, color="Plot", shape="Season")+
  #stat_ellipse(geom="polygon", aes(fill=Plot), alpha=0.05)+
  #scale_color_gradient2(midpoint=10, low="green", mid="orange", high="blue")+
  scale_fill_manual(values=c("#009acd",   "#00cd9a","#ff7f01","#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01","#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  scale_shape_manual(values=c(16,1, 4, 8))+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Hellinger PCA - All Neemonika - 16S") #600x500 TIFF or 6x4 PDF

#PERMANOVA to test for a statistically significant differences between fire treatment:
adonis2(otu_table(ps16S.Neemonika.hel) ~ sample_data(ps16S.Neemonika)$Fire, permutations = 1000) 
#                                     Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# sample_data(ps16S.Neemonika)$Fire   1     1.752  1.7517  4.7089 0.01827 0.000999 ***
# Residuals                         253    94.116  0.3720         0.98173             
# Total                             254    95.868                 1.00000
#
## DATA FILTERED FOR SAMPLES WITH >15000reads!
#                                     Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# sample_data(ps16S.Neemonika)$Fire   1     1.611 1.61084  4.5719 0.02518 0.000999 ***
# Residuals                         177    62.363 0.35233         0.97482             
# Total                             178    63.974                 1.00000

#
##
#
##
#
##
#

# Extract OTU table:
otu16S.df <- as.data.frame(otu_table(ps16S.Neemonika))
otu16S.df[1:5,1:5]


#calculate species RICHNESS
# Where ASV abundance value is greater than 0, give it a 1, then sum the 1's by sample:
otu16S.df.rich <- apply(otu16S.df > 0, 1, sum)
head(otu16S.df.rich)
otu16S.df.rich.dt <- as.data.table(otu16S.df.rich)
otu16S.df.rich.dt[1:5,]
dim(otu16S.df.rich.dt)


head(sample_data(ps16S.Neemonika))
MetaNeemonika16S <- as(sample_data(ps16S.Neemonika), "data.frame")
head(MetaNeemonika16S)
MetaNeemonika16S.dt <- as.data.table(MetaNeemonika16S)
head(MetaNeemonika16S.dt)
MetaNeemonika16S.dt[ ,SeqID:=rownames(MetaNeemonika16S)]

otu16S.df.rich.dt[ ,SeqID:=rownames(otu16S.df)]
otu16S.df.rich.dt <- merge(otu16S.df.rich.dt, MetaNeemonika16S.dt[,c("SeqID", "depth", "numdate", "Plot", "Fire")], by="SeqID")
head(otu16S.df.rich.dt)
names(otu16S.df.rich.dt)[2] <- "Richness"

otu16S.df.rich.dt$Plot <- factor(otu16S.df.rich.dt$Plot, levels=c("240", "321west", "321east", "400"))
otu16S.df.rich.dt$Fire <- factor(otu16S.df.rich.dt$Fire, levels=c("Pre-fire", "Post-fire"))

ggplot(otu16S.df.rich.dt, aes(x=numdate, y=Richness, color=Plot, fill=Plot))+
  geom_point(aes(shape=Fire), size=1)+
  #stat_summary(fun = "median", size = 10, geom = "point", shape = "_")+
  #stat_summary(aes(y = Richness, group = Plot), fun = median, geom="line", size = 0.5) +
  #stat_smooth_func(geom="text", method="lm", parse=TRUE) +
  #geom_smooth(aes(fill=Plot),  method="lm", size=1, alpha=0.1, span = 0.95)+
  geom_smooth(aes(fill=Plot),  method="lm", size=1, alpha=0.1)+
  theme_classic()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black", size=14))+
  theme(axis.text.y = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  scale_x_continuous(breaks=seq(17815, 18325, 30), limits=c(17815, 18325))+
  #scale_y_continuous(breaks=seq(0, 800, 100), limits=c(0, 800))+
  scale_shape_manual(values=c(5, 19))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  xlab("Time")+
  ylab("Richness")+
  ggtitle("16S Richness")
###export as .png or .tiff: 600x350, PDF: 6x4

#
###
#
##
#
##
#
####### PHILLIPS SAMPLES:
#
####### POST-FIRE comp240 
unique(sample_data(ps16S)$Date)
ps16S.F240pp <- prune_samples(sample_data(ps16S)$Date== "2020-02-17" & 
                                sample_data(ps16S)$Plot== "240" & 
                                sample_data(ps16S)$Who == "Phillip" | 
                                sample_data(ps16S)$Date== "2020-02-17" & 
                                sample_data(ps16S)$Plot== "240" & 
                                sample_data(ps16S)$Who == "PhillipPool", ps16S)
ps16S.F240pp #54 samples, 15026 taxa

# copy phyloseq object for Hellinger Transformation
ps16S.F240pp.hel <- ps16S.F240pp
# Hellinger Transformation
otu_table(ps16S.F240pp.hel) <-otu_table(decostand(otu_table(ps16S.F240pp.hel), method = "hellinger"), taxa_are_rows=FALSE)
otu_table(ps16S.F240pp.hel)[1:5,1:5]

# PCA
ps16S.F240pp.hel.pca <- rda(otu_table(ps16S.F240pp.hel))
summary(ps16S.F240pp.hel.pca)

#### PCA plot by depth
plot_ordination(ps16S.F240pp, ps16S.F240pp.hel.pca, color="depth")+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=depth), alpha=0.1)+
  theme_classic()+
  geom_point(size=2) + 
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  ggtitle("16S Hellinger PCA - February (no burn control) comp240") 
#TIF export: 600x500 or PDF export: 5x4

#PERMANOVA to test for a statistically significant differences between reps:
adonis2(otu_table(ps16S.F240pp.hel) ~ sample_data(ps16S.F240pp)$depth) #p=0.289

#
#
######### POST-FIRE comp321
#
# subset data
unique(sample_data(ps16S)$Date)
ps16S.F321pp <- prune_samples(sample_data(ps16S)$Date== "2020-02-17" & 
                                sample_data(ps16S)$Plot== "321west" & 
                                sample_data(ps16S)$Who == "Phillip" | 
                                sample_data(ps16S)$Date== "2020-02-17" & 
                                sample_data(ps16S)$Plot== "321west" & 
                                sample_data(ps16S)$Who == "PhillipPool", ps16S)

ps16S.F321pp #53 samples, 7622 taxa
# copy phyloseq object for Hellinger Transformation
ps16S.F321pp.hel <- ps16S.F321pp
# Hellinger Transformation
otu_table(ps16S.F321pp.hel) <-otu_table(decostand(otu_table(ps16S.F321pp.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps16S.F321pp.hel.pca <- rda(otu_table(ps16S.F321pp.hel))

# PCA plot by depth:
plot_ordination(ps16S.F321pp.hel, ps16S.F321pp.hel.pca, color="depth")+
  geom_mark_hull(concavity = 5,expand=0,radius=0,alpha=0.1,aes(fill=depth))+
  geom_point(size=2) + 
  theme_classic()+
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  ggtitle("16S Hellinger PCA - February (post-burn) comp321") 
#TIF export: 600x500 or PDF export: 5x4
# PERMANOVA to test for a statistically significant differences between depth:
adonis2(otu_table(ps16S.F321pp.hel) ~ sample_data(ps16S.F321pp)$depth) #p=0.5

#
#
#
###### PRE-FIRE 240
#
# subset data
unique(sample_data(ps16S)$Date)
ps16S.O240pp <- prune_samples(sample_data(ps16S)$Date== "2019-10-08" & 
                                sample_data(ps16S)$Plot== "240" & 
                                sample_data(ps16S)$Who == "Phillip" | 
                                sample_data(ps16S)$Date== "2019-10-08" & 
                                sample_data(ps16S)$Plot== "240" & 
                                sample_data(ps16S)$Who == "PhillipPool", ps16S)

# copy phyloseq object for Hellinger Transformation
ps16S.O240pp.hel <- ps16S.O240pp
# Hellinger Transformation
otu_table(ps16S.O240pp.hel) <-otu_table(decostand(otu_table(ps16S.O240pp.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps16S.O240pp.hell.pca <- rda(otu_table(ps16S.O240pp.hel))

# PCA plot by depth
plot_ordination(ps16S.O240pp.hel, ps16S.O240pp.hell.pca, color="depth")+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=depth))+ 
  geom_point(size=2) + 
  theme_classic()+
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  ggtitle("16S Hellinger PCA - October (no burn control) comp240") 
#TIF export: 600x500 or PDF export: 5x4
#PERMANOVA to test for a statistically significant differences between depths:
adonis2(otu_table(ps16S.O240pp.hel) ~ sample_data(ps16S.O240pp)$depth) 
# p-value=0.183


#
#
#
######### PRE-FIRE comp321
#
# subset data
unique(sample_data(ps16S)$Date)
ps16S.O321pp <- prune_samples(sample_data(ps16S)$Date== "2019-10-08" & 
                                sample_data(ps16S)$Plot== "321west" & 
                                sample_data(ps16S)$Who == "Phillip" | 
                                sample_data(ps16S)$Date== "2019-10-08" & 
                                sample_data(ps16S)$Plot== "321west" & 
                                sample_data(ps16S)$Who == "PhillipPool", ps16S)
ps16S.O321pp #51 samples, 7622 taxa

# copy phyloseq object for Hellinger Transformation
ps16S.O321.hel <- ps16S.O321pp
# Hellinger Transformation
otu_table(ps16S.O321.hel) <-otu_table(decostand(otu_table(ps16S.O321.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps16S.O321.hel.pca <- rda(otu_table(ps16S.O321.hel))
#
# PCA by depth
plot_ordination(ps16S.O321.hel, ps16S.O321.hel.pca, color="depth")+
  #geom_polygon(aes(fill=depth), alpha=0.2) + 
  geom_mark_hull(concavity = 5, expand=0, radius=0, aes(fill=depth))+
  geom_point(size=2) + 
  theme_classic()+
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black"))+
  ggtitle("16S Hellinger PCA - October (pre-burn) comp321") 
#TIF export: 600x500 or PDF export: 5x4
#PERMANOVA to test for a statistically significant differences between rep:
adonis2(otu_table(ps16S.O321.hel) ~ sample_data(ps16S.O321pp)$depth)
# p=0.114


##
##
###### PCAs on POOLED SAMPLES ONLY!
#
#Oct321
ps16S.O321.p <- prune_samples(sample_data(ps16S.O321pp)$depth != "20", ps16S.O321pp)
ps16S.O321.p #12 samples, 15026 taxa
# copy phyloseq object for Hellinger Transformation
ps16S.O321.p.hel <- ps16S.O321.p
# Hellinger Transformation
otu_table(ps16S.O321.p.hel) <-otu_table(decostand(otu_table(ps16S.O321.p.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps16S.O321.p.hel.pca <- rda(otu_table(ps16S.O321.p.hel))

# PCA by depth
plot_ordination(ps16S.O321.p.hel, ps16S.O321.p.hel.pca, color="depth")+
  #geom_polygon(aes(fill=depth), alpha=0.2) + 
  geom_mark_hull(concavity = 5, expand=0, radius=0, aes(fill=depth))+
  geom_point(size=2) + 
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#0077ee", "black"))+
  ggtitle("Hellinger PCA - October (pre-burn) comp321 -- pools only") #600x500
#PERMANOVA to test for a statistically significant differences between rep:
adonis2(otu_table(ps16S.O321.p.hel) ~ sample_data(ps16S.O321.p)$depth)
# p=0.73

#Oct240
ps16S.O240.p <- prune_samples(sample_data(ps16S.O240pp)$depth != "20", ps16S.O240pp)
ps16S.O240.p #12 samples, 7622 taxa
# copy phyloseq object for Hellinger Transformation
ps16S.O240.p.hel <- ps16S.O240.p
# Hellinger Transformation
otu_table(ps16S.O240.p.hel) <-otu_table(decostand(otu_table(ps16S.O240.p.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps16S.O240.p.hel.pca <- rda(otu_table(ps16S.O240.p.hel))

# PCA by depth
plot_ordination(ps16S.O240.p.hel, ps16S.O240.p.hel.pca, color="depth")+
  #geom_polygon(aes(fill=depth), alpha=0.2) + 
  geom_mark_hull(concavity = 5, expand=0, radius=0, aes(fill=depth))+
  geom_point(size=2) + 
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#0077ee", "black"))+
  ggtitle("Hellinger PCA - October (no burn control) comp240 -- pools only") #600x500
#PERMANOVA to test for a statistically significant differences between rep:
adonis2(otu_table(ps16S.O240.p.hel) ~ sample_data(ps16S.O240.p)$depth)
# p=0.92

#Feb321
ps16S.F321.p <- prune_samples(sample_data(ps16S.F321pp)$depth != "20", ps16S.F321pp)
ps16S.F321.p #51 samples, 7622 taxa
# copy phyloseq object for Hellinger Transformation
ps16S.F321.p.hel <- ps16S.F321.p
# Hellinger Transformation
otu_table(ps16S.F321.p.hel) <-otu_table(decostand(otu_table(ps16S.F321.p.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps16S.F321.p.hel.pca <- rda(otu_table(ps16S.F321.p.hel))

# PCA by depth
plot_ordination(ps16S.F321.p.hel, ps16S.F321.p.hel.pca, color="depth")+
  #geom_polygon(aes(fill=depth), alpha=0.2) + 
  geom_mark_hull(concavity = 5, expand=0, radius=0, aes(fill=depth))+
  geom_point(size=2) + 
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#0077ee", "black"))+
  ggtitle("Hellinger PCA - February (post-burn) comp321 -- pools only") #600x500
#PERMANOVA to test for a statistically significant differences between rep:
adonis2(otu_table(ps16S.F321.p.hel) ~ sample_data(ps16S.F321.p)$depth)
# p=0956

#Feb240
ps16S.F240.p <- prune_samples(sample_data(ps16S.F240pp)$depth != "20", ps16S.F240pp)
ps16S.F240.p #51 samples, 7622 taxa
# copy phyloseq object for Hellinger Transformation
ps16S.F240.p.hel <- ps16S.F240.p
# Hellinger Transformation
otu_table(ps16S.F240.p.hel) <-otu_table(decostand(otu_table(ps16S.F240.p.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps16S.F240.p.hel.pca <- rda(otu_table(ps16S.F240.p.hel))

# PCA by depth
plot_ordination(ps16S.F240.p.hel, ps16S.F240.p.hel.pca, color="depth")+
  #stat_ellipse(geom="polygon", aes(fill=depth), alpha=0.2)+ ##requires at least 4 points!
  #geom_polygon(aes(fill=depth), alpha=0.2) +
  #geom_density2d(alpha=0.5)+
  theme_classic()+
  geom_mark_hull(concavity = 5, expand=0, radius=0, aes(fill=depth))+
  geom_point(size=2) + 
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#0077ee", "black"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#0077ee", "black"))+
  ggtitle("Hellinger PCA - February (no burn control) comp240 -- pools only") #600x500
#PERMANOVA to test for a statistically significant differences between rep:
adonis2(otu_table(ps16S.F240.p.hel) ~ sample_data(ps16S.F240.p)$depth)
# p=0.678

##
##
###### PCAs TO COMPARE PRE-FIRE vs. POST-FIRE 
#
# subset data for everything except blanks (data without dates = blanks)

ps2.321prepost.pp <- prune_samples(sample_data(ps2)$Date== "17-Feb-20" & 
                                     sample_data(ps2)$Plot== "321west" & 
                                     sample_data(ps2)$Who == "PhillipPool" &
                                     sample_data(ps2)$depth == "p1-10" | 
                                     sample_data(ps2)$Date== "8-Oct-19" & 
                                     sample_data(ps2)$Plot== "321west" & 
                                     sample_data(ps2)$depth == "p1-10" &
                                     sample_data(ps2)$Who == "PhillipPool", ps2)
ps2.321prepost.pp #234 samples, 15026 taxa
# copy phyloseq object for Hellinger Transformation
ps2.321prepost.pp.hel <- ps2.321prepost.pp
# Hellinger Transformation
otu_table(ps2.321prepost.pp.hel) <-otu_table(decostand(otu_table(ps2.321prepost.pp.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps2.321prepost.pp.hel.pca <- rda(otu_table(ps2.321prepost.pp.hel))
# Broken Stick plot
screeplot(ps2.321prepost.pp.hel.pca, bstick = TRUE, 
          npcs = length(ps2.321prepost.pp.hel.pca$CA$eig),
          main = "Hellinger PCA - All Oct & Feb Phillip's Samples") #1500x600

# PCA by Date
plot_ordination(ps2.321prepost.pp.hel, ps2.321prepost.pp.hel.pca, color="Date")+
  geom_mark_hull(concavity = 5, expand=0, radius=0, aes(fill=Date))+
  geom_point(size=2) + 
  theme_classic() +
  ggtitle("Hellinger PCA: Pre/Post-fire Pooled Phillip Samples 1-10cm") #600x500
#PERMANOVA to test for a statistically significant differences between rep:
adonis(otu_table(ps2.321prepost.pp.hel) ~ sample_data(ps2.321prepost.pp)$Date)

colnames(sample_data(ps16S))


#
#
#





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
ITS16Stax <- rbind(TAXtableITS.sub, TAXtable16S.sub)


# CREATE PHYLOSEQ OBJECT of combined data
#
# Subset metadata for only the samples included in the OTU table:
MetaNeemonika <- MetadataWithReads[Who =="Neemonika" & numdate < 18308 & Plot != "380"]
ITS16S.MetaNeemonika <- MetaNeemonika[MetaNeemonika$SampID %in% rownames(ITS16Sotu), ]
dim(ITS16S.MetaNeemonika) #255 x 26
ITS16S.MetaNeemonika.rn <-column_to_rownames(ITS16S.MetaNeemonika, var="SampID")

psITS16S <- phyloseq(otu_table(ITS16Sotu, taxa_are_rows=FALSE), #matrix (rownames = SampIDs, colnames = ASVs, values = abundance of each ASV in each SeqID)
                     tax_table(ITS16Stax), #matrix (rownames = ASVs, colnames = taxonomic levels, values = taxonomic classification for each ASV)
                     sample_data(ITS16S.MetaNeemonika.rn)) #data.frame, (rownames = SampIDs, colnames & values = additional info for each SeqID)



###############################
### Differential Abundance ####
######## metagenomeSeq  #######
#########  Figure 3 ###########
###############################
library(phyloseq)
library(metagenomeSeq)
library(biomformat)

### below adapted from Thea Whitman's github!
#
# remove samples with less than 10000 reads because most blanks and KC's have <10000reads
psITS16S = prune_samples(sample_sums(psITS16S)>=10000, psITS16S) #nothing removed!

# Copy phyloseq object
psITS16S.biom = psITS16S

# turn phyloseq into a biom table using the biomformat package
biom = make_biom(data = t(otu_table(psITS16S.biom)), 
                 observation_metadata= tax_table(psITS16S.biom), 
                 sample_metadata = sample_data(psITS16S.biom))
#took about 5min to run

# turns our biom file into the type of file needed for this analysis (an MRexperiment object)
biom.MRexp = biom2MRexperiment(biom) #takes about a minute
biom.MRexp

#Preparing object for metegenomeseq analysis
MRexp = biom.MRexp 
MRexp = cumNorm(MRexp, cumNormStat(MRexp, qFlag=TRUE, pFlag = FALSE)) #takes about a minute
ModelData = pData(MRexp)

# Setting factors to be proper values
ModelData$Fire = as.factor(ModelData$Fire)
ModelData$Fire = ordered(ModelData$Fire, levels = c('Pre-fire','Post-fire'))
# Establishing model formula
model = model.matrix(~Fire, data=ModelData)
# Assigning settings
settings = zigControl(tol = 1e-04, maxit = 30, verbose = TRUE, dfMethod = "default", pvalMethod = "default")
# Running fit
fit = fitZig(obj=MRexp, mod=model, control = settings, useCSSoffset = TRUE, zeroMod = NULL, useMixedModel = FALSE)
# useCSSoffset = default, and automatically incudes the CSS scaling normalization factor
# it= 0, nll=366.62, log10(eps+1)=Inf, stillActive=110138
# it= 1, nll=363.20, log10(eps+1)=Inf, stillActive=106870 ...@11:21
# it= 2, nll=358.54, log10(eps+1)=Inf, stillActive=106379 ...@11:22
# it= 3, nll=358.54, log10(eps+1)=Inf, stillActive=106379 ...@11:23
# ....
# it=27, nll=341.30, log10(eps+1)=Inf, stillActive=48316 ....@11:41
# it=28, nll=341.31, log10(eps+1)=Inf, stillActive=48296 ....@11:41
# it=29, nll=341.31, log10(eps+1)=Inf, stillActive=48284 ....@11:42

#
# Effective sample size is calculated, and the average values of this is determined
EffSamp = calculateEffectiveSamples(fit)
MeanEffSamp = mean(EffSamp[!is.na(EffSamp)])
# These are the taxa that that had less than the average number of effective samples
# As recommended in the vignette: https://www.bioconductor.org/packages/devel/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf
rareFeatures = which(rowSums(MRcounts(MRexp) > 0) < MeanEffSamp)
dim(rareFeatures) # NULL
# SKIP THIS STEP, since there are no raraFeaturesTake the data object and remove the rareFeatures (taxa)
# MRexp = MRexp[-rareFeatures, ]

head(fit@fit)

# Extracting the results of interest
modeldesign = fit@fit$design
# Error in fit$fit : $ operator not defined for this S4 class
# fix = changed fit$fit$design to fit@fit$design

modelfit = fit@fit
# Error in fit$fit : $ operator not defined for this S4 class
# fix = changed fit$fit to fit@fit

#convert NAs to zeros (required for treat on teh next line)
modelfit$coefficients[is.na(modelfit$coefficients)] <- 0

# Empirical Bayes statistics for differential expression:
# computes t, F, and log-odds of differential expression
modelfit.treat = treat(modelfit, lfc=0)

dim(modelfit.treat$coefficients) #110138

#extract Differential Abundance Table:
ITS16S_DAtable = topTreat(modelfit.treat, number=110138)
# removed "coef=8" ..."column number or column name specifying..."
# which coefficient or contrast of the linear model is of interest"
# "number" = maximum number of genes to list
# number = 15026 is the complete list for ITS, 110138 for both ITS16S
# reduce this number to only see the top most DA taxa

#
##
#### MA plot (Figure 3)
#
##
ITS16S_DAtable.na <- na.omit(ITS16S_DAtable[,c(1,2,5)])
# add taxonomy
temp <- rownames_to_column(ITS16S_DAtable.na, var="ASV")
colnames(ITS16Stax.df.rn)
ITS16S_DAtable.na.TAX <- merge(temp, ITS16Stax.df.rn, by="ASV")

dim(ITS16S_DAtable.na.TAX[ITS16S_DAtable.na.TAX$Domain=="Bacteria",])
dim(ITS16S_DAtable.na.TAX[ITS16S_DAtable.na.TAX$Domain=="Fungi",])
colnames(ITS16S_DAtable.na.TAX.dt)
ITS16S_DAtable.na.TAX.dt <- data.table(ITS16S_DAtable.na.TAX)

ggplot(ITS16S_DAtable.na.TAX.dt[ITS16S_DAtable.na.TAX.dt$Domain=="Fungi",], 
       aes(x=AveExpr,y=logFC))+ 
  geom_point(aes(color=adj.P.Val<0.01), size=2, alpha=0.25)+
  scale_color_manual(values=setNames(c('#273046','#FAD510'), c(T,F)),
                     labels=c("p.adj>0.01","p.adj<0.01"),) +
  geom_hline(yintercept = c(-2,2), color="#CB2314", linetype=2)+ #turquoise lines at y=2 and y=-2
  geom_hline(yintercept = 0, color="black", alpha=0.5)+ #transparent line at y=0 
  theme_classic()+
  labs(color="")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  scale_x_log10()+  #log-scaled x-axis
  scale_y_continuous(limits=c(-18, 28), breaks=seq(-18, 28, 2))+
  ylab("Log2 Fold Change")+ #label for y-axis
  xlab("Average Adundance")+ #label for x-axis
  ggtitle("FUNGAL Differential Abundance")

ggplot(ITS16S_DAtable.na.TAX.dt[ITS16S_DAtable.na.TAX.dt$Domain=="Bacteria",], 
       aes(x=AveExpr,y=logFC))+ 
  geom_point(aes(color=adj.P.Val<0.01), size=2, alpha=0.15)+
  scale_color_manual(values=setNames(c('#273046','#FAD510'), c(T,F)),
                     labels=c("p.adj>0.01","p.adj<0.01"),) +
  geom_hline(yintercept = c(-2,2), color="#CB2314", linetype=2)+ #turquoise lines at y=2 and y=-2
  geom_hline(yintercept = 0, color="black", alpha=0.5)+ #transparent line at y=0 
  theme_classic()+
  labs(color="")+
  guides(color = guide_legend(override.aes = list(size=3)))+
  scale_x_log10()+  #log-scaled x-axis
  scale_y_continuous(limits=c(-22, 20), breaks=seq(-22, 20, 2))+
  ylab("Log2 Fold Change")+ #label for y-axis
  xlab("Average Adundance")+ #label for x-axis
  ggtitle("BACTERIAL Differential Abundance")
#export 5x3


# Explore how Fungal FunGuilds relate to DA...
ITS16S_DAtable.na.TAXguild <- merge(ITS16S_DAtable.na.TAX.dt, TaxGuild, by="ASV", all.x=TRUE)

unique(TaxGuild$trophicMode)

ITS_DAtaxGuild <- ITS16S_DAtable.na.TAXguild[ITS16S_DAtable.na.TAXguild$Domain=="Fungi" ,]

ggplot(ITS_DAtaxGuild[ITS_DAtaxGuild$trophicMode=="Saprotroph",], 
       aes(x=AveExpr,y=logFC))+ 
  geom_point(aes(color=trophicMode), size=2)+
  geom_hline(yintercept = c(-2,2), color="black")+ #turquoise lines at y=2 and y=-2
  geom_hline(yintercept = 0, color="grey50", alpha=0.5)+ #transparent line at y=0 
  scale_color_manual(values=c("#82ffc1"))+
  #scale_color_manual(values=c("grey80", "#ff9b35", "#ff9b35", "grey80", "#9b35ff",
  #                                    "#82ffc1", "#4ea5ff", "#9b35ff"))+
  theme_classic()+
  labs(color="")+
  #guides(color = guide_legend(override.aes = list(size=3)))+
  scale_x_log10()+  #log-scaled x-axis
  scale_y_continuous(limits = c(min(ITS_DAtaxGuild$logFC), max(ITS_DAtaxGuild$logFC)))+
  ylab("Log2 Fold Change")+ #label for y-axis
  xlab("Average Adundance")+ #label for x-axis
  ggtitle("FUNGAL Differential Abundance")

view(TaxGuild[trophicMode == "Pathotroph"])
view(TaxGuild[trophicMode == "Pathotroph-Saprotroph"]) #plant/animal pathogens
view(TaxGuild[trophicMode == "Pathotroph-Saprotroph-Symbiotroph"]) #catch-all.. pretty meaningless..
view(TaxGuild[trophicMode == "Pathotroph-Symbiotroph"]) # mostly ericoid mycorrhizae, some lichen/ednophyte/epiphyte
view(TaxGuild[trophicMode == "Saprotroph-Symbiotroph"]) #maybe mycorrhizal / maybe just a saprotroph


head(ITS16S_DAtable.na)
#add a column to summarize the if, and the direction of DA for each ASV
ITS16S_DAtable.summary <- ITS16S_DAtable.na %>%
  mutate(DAsummary = case_when(
    adj.P.Val>0.1 ~ "NotDA",
    adj.P.Val<0.1 & logFC > 2 ~ "PosDA", 
    adj.P.Val<0.1 & logFC < -2 ~ "NegDA",
    adj.P.Val<0.1 & logFC > -2 & logFC < 2 ~ "NotDA",
    is.na(adj.P.Val) ~ NA_character_, 
    is.na(logFC) ~ NA_character_
  ))


##################################################
####### TITAN pilot run with ITS data only #######
################# Figure 4 #######################
##################################################
library(TITAN2)

# TITAN2 requires two inputs:
# (1) OTU table -- data.frame: rownames = samples, colnames=OTUs
view(glades.taxa) #example
# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
view(glades.env) #example

# Extract OTU table:
ps2.Neemonika <- prune_samples(sample_data(ps2)$Who == "Neemonika" & 
                                 sample_data(ps2)$Plot != "380" &
                                 sample_data(ps2)$numdate < 18308, ps2)

otu.df <- as.data.frame(otu_table(ps2))
otu.df[1:5,1:5]
dim(otu.df)
# extract Taxonomy table:
tax.df <- as.data.frame(tax_table(ps2))
tax.df[1:5,]
dim(tax.df)
tax.df.rn <- rownames_to_column(tax.df, var="ASV")
#transpose OTU table:
otu.dft <- as.data.frame(t(otu.df))
otu.dft[1:5,1:5]
# Merge taxa info onto OTU table
OTUtax <- as.data.table(merge(otu.dft, tax.df, by=0))
dim(OTUtax)
OTUtax[1:3,665:672]

# Create a NEW AND IMPROVED .taxa table with Genus appened onto ASVs,
# which will make the TITAN output plots a little easier to digest.
# split the "Genus" column to get rid of the "g__" before every genus name
OTUtax[ ,c("garbage", "JustGenus"):=tstrsplit(Genus, "__", fixed=TRUE)]
# create column with GenusOTU info (i.e. "Piloderma_OTU1")
OTUtax <- unite(OTUtax, GenusASV, c(JustGenus, Row.names), sep="_", remove=FALSE)
# Massage table into the format that TITAN is expecting:
#remove excess columns (by number) so that all that remains are the sample columns and the GenusOTU column
OTUtax[,c(2,667:675)] <- list(NULL)
OTUtax <- column_to_rownames(OTUtax, var="GenusASV")
#transpose
OTUtax.t <- as.data.frame(t(OTUtax))
OTUtax.t.rn <- rownames_to_column(OTUtax.t, var="SeqID")
OTUtax.t.rn[1:5, 1:5] ## THIS is the table that should be used as the basis for .taxa arguments in TITAN!
#
# Since TITAN identifies indicator taxa along a gradient,
# I think it makes sense to run TITAN on each site independently to 
# identify indicators along the timeline within each site.
#
#
nrow(Metadata400)

########### comp400 TITAN! 
#
# (1) OTU Table -- Subset OTUtax.t.rn table for only the column names that match the SeqIDs in Metadata400
BS400.taxa <- OTUtax.t.rn[OTUtax.t.rn$SeqID %in% Metadata400$SeqID, ]
BS400.taxa[1:5,1:5]
rownames(BS400.taxa) <- c() #explicitly assign nothing to the rownames, so that column_to_rownames has a blank slate to work with...
BS400.taxa <- column_to_rownames(BS400.taxa, var="SeqID")
# remove OTUs that occur in only 3 or fewer samples (TITAN requirement)
BS400.taxa <- BS400.taxa[colSums(BS400.taxa > 0) >= 3 ]
# "BS400.taxa > 0" creates a TRUE/FALSE table for abundance values >0, 
# colSums() adds up all the values in each column (TRUE = 1, FALSE = 0)
# subset the table for where their are more than 3 TRUE (non-zero) values!
dim(BS400.taxa)
# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
# Subset the Metadata400 table to create the .env argument:
BS400.env <- Metadata400[ ,c("SeqID", "numdate")]
BS400.env <- column_to_rownames(BS400.env, var="SeqID")

#check that everything looks right before launching TITAN
# .env argument should be numeric:
class(BS400.env$numdate)
BS400.env$numdate <- as.numeric(BS400.env$numdate)
class(BS400.env$numdate)
# .env and .taxa tables should have equivallent rownames
head(BS400.env)
BS400.taxa[1:5,1:5]
# both tables should have the same number of rows (first number)
dim(BS400.env)
dim(BS400.taxa)

# Run TITAN2 with default parameters:
# WARNING - it can take hours to run!
#example:: glades.titan <- titan(glades.env, glades.taxa) 
BS400.titan <- titan(BS400.env, BS400.taxa) 
# start @ 8:23am, end @ ~10:00am
#
# "100% occurrence detected 1 times (0.1% of taxa)" 
# "Proportion of pure and reliable taxa = 0.0431" 

##EXPLORE THE RESULTS OF TITAN!
#
# overview of what's contained within out TITAN object:
summary(BS400.titan)
# Each item within the TITAN object can be accessed as if it were a column
# ...for example: BS400.titan$sumz.cp
#
#summary(BS400.titan)
#           Length  Class  Mode   
#sppmax        24096 -none- numeric   ### table of change points by taxa, dim = 1506 x 16
#sumz.cp          24 -none- numeric   ### table of observed change points ("cp"), with selected quantiles, dim = 4 x 6
#env              80 -none- numeric   ### input env table, dim = 80 x 1
#taxa         120480 -none- numeric   ### input taxa table, dim = 80 x 1506
#envcls           71 -none- numeric   ### "A vector of candidate partitions derived from subtracting ’minSplt’ from ’env’"
#srtEnv           80 -none- numeric   ### sorted version of input env table
#ivzScores    427704 -none- numeric   ### "matrix containing group membership, z scores, IndVals, and p values for each taxon at every candidate partition in ’envcls’"
#ivz             142 -none- numeric   ### sum z+ and sum z- scores for each value in encvls
#ivz.f           142 -none- numeric   ### Filtered sum z+ and sum z- scores for each value in encvls
#maxSumz        1000 -none- numeric   ### "matrix of environmental values at sum(z-) and sum(z+) maxima across all bootstrap replicates
#maxFsumz       1000 -none- numeric   ### "matrix of environmental values at filtered sum(z-) and sum(z+) maxima across all bootstrap replicates"
#metricArray 3012000 -none- numeric   ### "An array of group membership, env change points, z scores, and p values for passing to ’plot.IVecdf’"
#arguments        10 -none- numeric   ### arguments used in TITAN function call

# summary table of observed change points ("cp"), with selected quantiles around the cp
BS400.titan$sumz.cp
#output:
#       cp      0.05      0.10    0.50    0.90    0.95
#sumz-  17904   17904.0   17946   18006   18127   18127
#sumz+  18286   18086.0   18100   18206   18305   18305
#fsumz- 18024   17988.0   17988   18024   18057   18057
#fsumz+ 18253   18222.5   18239   18267   18286   18286


# Table of TITAN results for all taxa
# Zscores, IndVal scores, p-values, etc...
head(BS400.titan$sppmax)
# SAVE THIS TABLE!
BS400.titan.dt <- as.data.table(BS400.titan$sppmax)
# Results table doesn't have informative row.names, so re-attach ASVs:
# transpose taxa table and then move rownames to the results table
BS400.taxa <- as.data.frame(t(BS400.taxa)) #BS.taxa was the input to TITAN that can also be seen as BS.titan$taxa
BS400.taxa.rn <- rownames_to_column(BS400.taxa, var="ASV")
# Paste ASVs onto the TITAN results table:
merge1 <- cbind(BS400.taxa.rn$ASV, BS400.titan.dt)
names(merge1)[1] <- "GenusASV"
# split the GenusASV column to make an ASV column that can be matched with the ASV column in tax.df
merge1[ ,c("Genus", "ASV"):=tstrsplit(GenusASV, "_", fixed=TRUE)]

tax.df[1:5, 1:5]
tax.df.rn <- rownames_to_column(tax.df, var="ASV")
tax.df.rn[1:5, 1:5]
# merge taxonomy table with results table:
BS400.titan.AllResults <- merge(merge1, tax.df.rn, by="ASV")

# filter the table for only the results that were ultimately plotted = where the value in the "filter" column is >0
BS400.titan.TopResults <- BS400.titan.AllResults[filter>0]

colnames(BS400.titan.TopResults)
# ienv.cp- environmental change point for each taxon based on IndVal maximum (used if imax = TRUE)
# zenv.cp- environmental change point for each taxon based on z maximum (default, imax = FALSE)
# freq- number of non-zero abundance values per taxon
# maxgrp- 1 if z- (negative response); 2 if z+ (positive response)
# IndVal-Dufrene and Legendre 1997 IndVal statistic, scaled 0-100%
# obsiv.prob- probability of an equal or larger IndVal from random permutation. **squishy! don't trust!
# zscore- IndVal z score
# 5%, 10%, 50%, 90%, 95%- change point quantiles among bootstrap replicates
# purity - proportion of replicates matching observed maxgrp assignment
# reliability - proportion of replicate obsiv.prob values < = 0.05
# z.median- median score magnitude across all bootstrap replicates
# filter- logical (if >0) indicating whether each taxa met purity and reliability criteria, value indicates maxgrp assignment.


#
##
#
#PLOT RESULTS!
#
# Summary plot of total z scores...
plot_sumz_density(BS400.titan)
#illustrates the filtered sumz's
plot_sumz(BS400.titan, xmin = 17800, xmax = 18400, col1 = "blue", col2 = "orange")
#
class(BS400.titan.TopResults$zenv.cp)
#plot z-scores of indicator taxa for the +/- change points! NEW METHOD (colored distributions)
plot_taxa_ridges(BS400.titan, xlabel = "Time")
#plot plot z-scores of indicator taxa for the +/- change points! OLD METHOD (B&W circle+line "boxplots")
plot_taxa(BS400.titan, xlabel = "Time",  xmin = 17800, xmax = 18400)
#generate Z-plot with ggplot2!!
ggplot(BS400.titan.TopResults, aes(x=zenv.cp, y=reorder(GenusASV, -zenv.cp), color=Phylum))+
  geom_point(aes(size=zscore))+
  geom_errorbar(aes(xmax=`95%`, xmin=`5%`))+
  theme_light()+
  scale_x_continuous(breaks=seq(17880, 18320, 20), limits=c(17880, 18320), expand=c(0,0))+
  theme(axis.text.x = element_text(angle=90, color="black", size=12, hjust=0, vjust=0.5))+
  theme(axis.text.y = element_text(color="black", size=12))+
  xlab("Time (days)")+
  ylab("Indicator Taxa")+
  ggtitle("comp400 TITAN results - ASVs")

#plot abundance of each of the indicator taxa over time in a faceted manner:
plot_cps(BS400.titan)
#plot an individual taxon from the above plot
plot_cps(BS400.titan, taxaID = "ASV773", cp.trace = TRUE, xlabel = "Time")
plot_cps(BS400.titan, taxaID = "Pyronema_ASV363", cp.med=FALSE, xlabel = "Time")
#black circles = abundance at each time point (comes directly from input table)
#red circle = observed change point (if "cp.med=TRUE" = median change point across all bootstrap reps)
#blue = histogram of bootstrap replicates by probability density

#compare TITAN plot with ggplot-ed abundance values
plot_cps(BS400.titan, taxaID = "Pyronema_ASV363", cp.med=FALSE, xlabel = "Time")
dat <- merge(BS400otu[,c("SeqID", "ASV363")], MetadataWithReads[,c("SeqID", "numdate")], by="SeqID")
ggplot(dat, aes(x=numdate, y=ASV363))+
  geom_point(size=4, shape=1)+
  theme_minimal()


#generate Z-plot with ggplot2!!
ggplot(BS321Easv.titan.TopResults, aes(x=zenv.cp, y=reorder(GenusASV, -zenv.cp), color=Phylum))+
  geom_point(aes(size=zscore))+
  geom_errorbar(aes(xmax=`95%`, xmin=`5%`))+
  theme_light()+
  xlab("Time (days)")+
  ylab("Indicator Taxa")+
  ggtitle("comp321E TITAN results")



#
##
#
##
#
########### comp240 TITAN!
#
# (1) OTU Table -- Subset OTUtax.t.rn table for only the column names that match the SeqIDs in Metadata240
BS240.taxa <- OTUtax.t.rn[OTUtax.t.rn$SeqID %in% Metadata240$SeqID, ]
rownames(BS240.taxa) <- c() #explicitly assign nothing to the rownames, so that column_to_rownames has a blank slate to work with...
BS240.taxa <- column_to_rownames(BS240.taxa, var="SeqID")
# remove OTUs that occur in only 3 or fewer samples (TITAN requirement)
BS240.taxa <- BS240.taxa[colSums(BS240.taxa > 0) >= 3 ]
# "BS240.taxa > 0" creates a TRUE/FALSE table for abundance values >0, 
# colSums() adds up all the values in each column (TRUE = 1, FALSE = 0)
# subset the table for where their are more than 3 TRUE (non-zero) values!

# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
# Subset the Metadata400 table to create the .env argument:
BS240.env <- Metadata240[ ,c("SeqID", "numdate")]
BS240.env <- column_to_rownames(BS240.env, var="SeqID")

#check that everything looks right before launching TITAN
# .env argument should be numeric:
class(BS240.env$numdate)
BS240.env$numdate <- as.numeric(BS240.env$numdate)
class(BS240.env$numdate)
# .env and .taxa tables should have equivallent rownames
head(BS240.env)
BS240.taxa[1:5,1:5]
# both tables should have the same number of rows (first number)
dim(BS240.env)
dim(BS240.taxa)

# Run TITAN2 with default parameters:
# WARNING - it can take hours to run!
#example:: glades.titan <- titan(glades.env, glades.taxa) 
BS240.titan <- titan(BS240.env, BS240.taxa) 
# start @10:20am, end @ 11:38am
#
#Proportion of pure and reliable taxa = 0.0244
#         cp    0.05    0.10    0.50    0.90    0.95
#sumz-  17855 17835.5 17855 17863.0 17868 17873.0
#sumz+  18239 18158.5 18177 18222.5 18253 18267.0
#fsumz- 17862 17858.5 17862 17878.0 18006 18023.4
#fsumz+ 18239 18176.8 18177 18206.0 18239 18253.0


# Table of TITAN results for all taxa
# Zscores, IndVal scores, p-values, etc...
head(BS240.titan$sppmax)
# SAVE THIS TABLE!
BS240.titan.dt <- as.data.table(BS240.titan$sppmax)
# Results table doesn't have informative row.names, so re-attach ASVs:
# transpose taxa table and then move rownames to the results table
BS240.taxa <- as.data.frame(t(BS240.taxa)) #BS.taxa was the input to TITAN that can also be seen as BS.titan$taxa
BS240.taxa.rn <- rownames_to_column(BS240.taxa, var="ASV")

# Paste ASVs onto the TITAN results table:
merge1 <- cbind(BS240.taxa.rn$ASV, BS240.titan.dt)
names(merge1)[1] <- "GenusASV"
# split the GenusASV column to make an ASV column that can be matched with the ASV column in tax.df
merge1[ ,c("Genus", "ASV"):=tstrsplit(GenusASV, "_", fixed=TRUE)]
# merge taxonomy table with results table:
BS240.titan.AllResults <- merge(merge1, tax.df.rn, by="ASV")
head(BS240.titan.AllResults)

# filter the table for only the results that were ultimately plotted = where the value in the "filter" column is >0
BS240.titan.TopResults <- BS240.titan.AllResults[filter>0]

#
#
#PLOT RESULTS!
#
#plot z scores...
plot_sumz_density(BS240.titan)
#illustrates the filtered sumz's
#

#plot z-scores of indicator taxa for the +/- change points! NEW METHOD (colored distributions)
plot_taxa_ridges(BS240.titan, xlabel = "Time")
#plot plot z-scores of indicator taxa for the +/- change points! OLD METHOD (B&W circle+line "boxplots")
plot_taxa(BS240.titan, xlabel = "Time", xmin = 17800, xmax = 18400)
#if it says "figure margins too large", clear your plot history!!
#generate Z-plot with ggplot2!!
ggplot(BS240.titan.TopResults, aes(x=zenv.cp, y=reorder(GenusASV, -zenv.cp), color=Phylum))+
  geom_point(aes(size=zscore))+
  geom_errorbar(aes(xmax=`95%`, xmin=`5%`))+
  theme_light()+
  scale_x_continuous(breaks=seq(17820, 18300, 20), limits=c(17820, 18300), expand=c(0,0))+
  theme(axis.text.x = element_text(angle=90, color="black", size=12, hjust=0, vjust=0.5))+
  theme(axis.text.y = element_text(color="black", size=12))+
  xlab("Time (days)")+
  ylab("Indicator Taxa")+
  ggtitle("comp240 TITAN results")


#plot abundance of each of the indicator taxa over time in a faceted manner:
plot_cps(BS240.titan)
#plot an individual taxon from the above plot
plot_cps(BS240.titan, taxaID = "ASV773", cp.trace = TRUE, xlabel = "Time")
plot_cps(BS240.titan, taxaID = "ASV2770", xlabel = "Time")
#abundance at each time point = black circles


#
#
#
#
########### comp321east TITAN!
#
# Subset OTUtax table for only the column names that match the SeqIDs in Metadata400
BS321e.taxa <- OTUtax.t.rn[OTUtax.t.rn$SeqID %in% Metadata321E$SeqID, ]
rownames(BS321e.taxa) <- c() #explicitly assign nothing to the rownames, so that column_to_rownames has a blank slate to work with...
BS321e.taxa <- column_to_rownames(BS321e.taxa, var="SeqID")
# remove OTUs that occur in only 3 or fewer samples (TITAN requirement)
BS321e.taxa <- BS321e.taxa[colSums(BS321e.taxa > 0) >= 3 ]
# "BS321e.taxa > 0" creates a TRUE/FALSE table for abundance values >0, 
# colSums() adds up all the values in each column (TRUE = 1, FALSE = 0)
# subset the table for where their are more than 3 TRUE (non-zero) values!

# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
# Subset the Metadata321E table to create the .env argument:
BS321e.env <- Metadata321E[ ,c("SeqID", "numdate")]
BS321e.env <- column_to_rownames(BS321e.env, var="SeqID")

#check that everything looks right before launching TITAN
# .env argument should be numeric:
class(BS321e.env$numdate)
BS321e.env$numdate <- as.numeric(BS321e.env$numdate)
class(BS321e.env$numdate)
# .env and .taxa tables should have equivallent rownames
head(BS321e.env)
BS321e.taxa[1:5,1:5]
# both tables should have the same number of rows (first number)
dim(BS321e.env)
dim(BS321e.taxa)

# Run TITAN2 with default parameters:
# WARNING - it can take hours to run!
#example:: glades.titan <- titan(glades.env, glades.taxa) 
BS321e.titan <- titan(BS321e.env, BS321e.taxa) 
# start @ 12:41... end @ ~7:20ish
# 100% occurrence detected 24 times (1.2% of taxa)
# "Proportion of pure and reliable taxa = 0.0504
#       cp    0.05    0.10  0.50  0.90    0.95
#sumz-  17855 17833.0 17844 17855 17868 17888
#sumz+  18267 18071.5 18086 18239 18267 18286
#fsumz- 17862 17844.0 17855 17862 17865 17868
#fsumz+ 18086 17888.0 18057 18100 18114 18140


# Table of TITAN results for all taxa
# Zscores, IndVal scores, p-values, etc...
head(BS321e.titan$sppmax)
# SAVE THIS TABLE!
BS321e.titan.dt <- as.data.table(BS321e.titan$sppmax)
# Results table doesn't have informative row.names, so re-attach ASVs:
# transpose taxa table and then move rownames to the results table
BS321e.taxa.t <- as.data.frame(t(BS321e.taxa)) #BS.taxa was the input to TITAN that can also be seen as BS.titan$taxa
BS321e.taxa.t.rn <- rownames_to_column(BS321e.taxa.t, var="ASV")
# Paste ASVs onto the TITAN results table:
merge1 <- cbind(BS321e.taxa.t.rn$ASV, BS321e.titan.dt)
names(merge1)[1] <- "GenusASV"
# split the GenusASV column to make an ASV column that can be matched with the ASV column in tax.df
merge1[ ,c("Genus", "ASV"):=tstrsplit(GenusASV, "_", fixed=TRUE)]
# merge taxonomy table with results table:
BS321e.titan.AllResults <- merge(merge1, tax.df.rn, by="ASV")

# filter the table for only the results that were ultimately plotted = where the value in the "filter" column is >0
BS321e.titan.TopResults <- BS321e.titan.AllResults[filter>0]
#
#
#PLOT RESULTS!
#
#plot z scores...
plot_sumz_density(BS321e.titan)
#illustrates the filtered sumz's
#
#plot z-scores of indicator taxa for the +/- change points! NEW METHOD (colored distributions)
plot_taxa_ridges(BS321e.titan, xlabel = "Time")
#plot plot z-scores of indicator taxa for the +/- change points! OLD METHOD (B&W circle+line "boxplots")
plot_taxa(BS321e.titan, xlabel = "Time", xmin = 17800, xmax = 18400)
#if it says "figure margins too large", clear your plot history!!

#generate Z-plot with ggplot2!!
ggplot(BS321e.titan.TopResults, aes(x=zenv.cp, y=reorder(GenusASV, -zenv.cp), color=Phylum))+
  geom_point(aes(size=zscore))+
  geom_errorbar(aes(xmax=`95%`, xmin=`5%`))+
  theme_light()+
  scale_x_continuous(breaks=seq(17800, 18320, 40), limits=c(17800, 18320), expand=c(0,0))+
  theme(axis.text.x = element_text(angle=90, color="black", size=12, hjust=0, vjust=0.5))+
  theme(axis.text.y = element_text(color="black", size=12))+
  xlab("Time (days)")+
  ylab("Indicator Taxa")+
  ggtitle("comp321E TITAN results - ASVs")

#plot abundance of each of the indicator taxa over time in a faceted manner:
plot_cps(BS321e.titan)
# "figure margins too large" -- try clicking the broom icon in the plot viewer, 
# and/or try changing the physical dimensions of the plot viewer in R Studio
#plot an individual taxon from the above plot
plot_cps(BS321e.titan, taxaID = "ASV2517", cp.trace = TRUE, xlabel = "Time")
plot_cps(BS321e.titan, taxaID = "ASV2517", xlabel = "Time")
#abundance at each time point = black circles

# Zplot with ggplot
ggplot(BS321e.titan.TopResults, aes(x=zenv.cp, y=reorder(GenusOTU, -zenv.cp), color=Phylum))+
  geom_point(aes(size=zscore))+
  geom_errorbar(aes(xmax=`95%`, xmin=`5%`))+
  theme_light()+
  xlab("Time (days)")+
  ylab("Indicator Taxa")+
  ggtitle("comp321E TITAN results - ASVs up until 321W burn")
#
#
#
#
#
########### comp321west TITAN!
#
dim(OTUtax.t.rn)
OTUtax.t.rn[1:5, 15020:15027]
# Subset OTUtax table for only the column names that match the SeqIDs in Metadata321w
BS321w.taxa <- OTUtax.t.rn[OTUtax.t.rn$SeqID %in% Metadata321W$SeqID, ]
dim(BS321w.taxa)
BS321w.taxa[1:5,1:5]
rownames(BS321w.taxa) <- c() #explicitly assign nothing to the rownames, so that column_to_rownames has a blank slate to work with...
BS321w.taxa <- column_to_rownames(BS321w.taxa, var="SeqID")
# remove OTUs that occur in only 3 or fewer samples (TITAN requirement)
BS321w.taxa <- BS321w.taxa[colSums(BS321w.taxa > 0) >= 3 ]
# "BS321w.taxa > 0" creates a TRUE/FALSE table for abundance values >0, 
# colSums() adds up all the values in each column (TRUE = 1, FALSE = 0)
# subset the table for where their are more than 3 TRUE (non-zero) values!

# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
# Subset the Metadata321E table to create the .env argument:
BS321w.env <- Metadata321W[ ,c("SeqID", "numdate")]
BS321w.env <- column_to_rownames(BS321w.env, var="SeqID")

#check that everything looks right before launching TITAN
# .env argument should be numeric:
class(BS321w.env$numdate)
BS321w.env$numdate <- as.numeric(BS321w.env$numdate)
class(BS321w.env$numdate)
# .env and .taxa tables should have equivallent rownames
head(BS321w.env)
BS321w.taxa[1:6,1:5]
# both tables should have the same number of rows (first number)
dim(BS321w.env)
dim(BS321w.taxa)

# Run TITAN2 with default parameters:
# WARNING - it can take hours to run!
#example:: glades.titan <- titan(glades.env, glades.taxa) 
BS321w.titan <- titan(BS321w.env, BS321w.taxa) 
# start @ 9:18am, end @ 10:09am 
# 100% occurrence detected 11 times (0.8% of taxa)
# Proportion of pure and reliable taxa = 0.0182
#       cp    0.05    0.10    0.50    0.90    0.95
#sumz-   17855.0 17837.5 17837.5 17865.0 18158.5 18158.5
#sumz+  18330.0 18177.0 18251.6 18309.0 18330.0 18330.0
#fsumz- 18191.5 17868.0 17873.0 18158.5 18206.0 18253.0
#fsumz+ 18253.0 18158.5 18177.0 18253.0 18307.0 18307.0


# Table of TITAN results for all taxa
# Zscores, IndVal scores, p-values, etc...
head(BS321w.titan$sppmax)
# SAVE THIS TABLE!
BS321w.titan.dt <- as.data.table(BS321w.titan$sppmax)
# Results table doesn't have informative row.names, so re-attach ASVs:
# transpose taxa table and then move rownames to the results table
BS321w.taxa.t <- as.data.frame(t(BS321w.taxa)) #BS.taxa was the input to TITAN that can also be seen as BS.titan$taxa
BS321w.taxa.t.rn <- rownames_to_column(BS321w.taxa.t, var="ASV")
# Paste ASVs onto the TITAN results table:
merge1 <- cbind(BS321w.taxa.t.rn$ASV, BS321w.titan.dt)
names(merge1)[1] <- "GenusASV"
# split the GenusOTU column to make an OTU column that can be matched with the OTU column in tax.df
merge1[ ,c("Genus", "ASV"):=tstrsplit(GenusASV, "_", fixed=TRUE)]
# merge taxonomy table with results table:
BS321w.titan.AllResults <- merge(merge1, tax.df.rn, by="ASV")

# filter the table for only the results that were ultimately plotted = where the value in the "filter" column is >0
BS321w.titan.TopResults <- BS321w.titan.AllResults[filter>0]
#

#
#PLOT RESULTS!
#
#plot z scores...
plot_sumz_density(BS321w.titan)
#illustrates the filtered sumz's
# if you get an error -- just clear your plotting window!!

#plot z-scores of indicator taxa for the +/- change points! NEW METHOD (colored distributions)
plot_taxa_ridges(BS321w.titan, xlabel = "Time")
#plot plot z-scores of indicator taxa for the +/- change points! OLD METHOD (B&W circle+line "boxplots")
plot_taxa(BS321w.titan, xlabel = "Time", xmin = 17800, xmax = 18400)
#if it gives an error (i.e. plot size is too large).. clear your plotting window!!
#generate Z-plot with ggplot2!!
ggplot(BS321w.titan.TopResults, aes(x=zenv.cp, y=reorder(GenusASV, -zenv.cp), color=Phylum))+
  geom_point(aes(size=zscore))+
  geom_errorbar(aes(xmax=`95%`, xmin=`5%`))+
  theme_light()+
  scale_x_continuous(breaks=seq(17860, 18300, 20), limits=c(17860, 18300), expand=c(0,0))+
  theme(axis.text.x = element_text(angle=90, color="black", size=12, hjust=0, vjust=0.5))+
  theme(axis.text.y = element_text(color="black", size=12))+
  xlab("Time (days)")+
  ylab("Indicator Taxa")+
  ggtitle("comp321W TITAN results - ASVs up until 321W burn")


#plot abundance of each of the indicator taxa over time in a faceted manner:
plot_cps(BS321w.titan)
#plot an individual taxon from the above plot
plot_cps(BS321w.titan, taxaID = "ASV1550", cp.trace = TRUE, xlabel = "Time")
plot_cps(BS321w.titan, taxaID = "ASV2517", xlabel = "Time")
#abundance at each time point = black circles
#
##
#
#########################################################################
######  Subset tables to create FastSpar inputs, and run FastSpar #######
############################## Figure 5 #################################
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
##################### Figure 5 #########################
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
################### Figure 5 #####################
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

##################################################
######  Summarizing Network & DA patterns #######
################ Figure 5A&B ####################
#################################################

colnames(nodesALL.DA)
dim(nodesALL.DA)
mean(nodesALL.DA$logFC, na.rm=TRUE) # -0.399
sd(nodesALL.DA$logFC, na.rm=TRUE)   # 3.1368
hist(nodesALL.DA$logFC, breaks=100)

#add a column to denote in which condition a taxon is DA
nodesALL.DA <- nodesALL.DA %>%
  mutate(DAsummary = case_when(
    adj.P.Val>0.5 ~ "NotDA",
    adj.P.Val<0.5 & logFC > 2 ~ "BurnDA", 
    adj.P.Val<0.5 & logFC < -2 ~ "ControlDA",
    adj.P.Val<0.5 & logFC > -2 & logFC < 2 ~ "NotDA",
    is.na(adj.P.Val) ~ NA_character_, 
    is.na(logFC) ~ NA_character_
  ))

unique(nodesALL.DA$DAsummary)
# for some reason there's still a few NA's... change them to NotDA:
nodesALL.DA[is.na(nodesALL.DA$DAsummary)] <- "NotDA"
unique(nodesALL.DA$DAsummary)

#add a column summarizing TITAN columns
nodesALL.DA.TITAN <- nodesALL.DA %>%
  mutate(TITANsummary = case_when(
    TITAN400 == "TITAN400" ~ "TITANind",
    TITAN321E == "TITAN321E" ~ "TITANind", 
    TITAN321W == "TITAN321W" ~ "TITANind",
    TITAN240 == "TITAN240" ~ "TITANind",
  ))
unique(nodesALL.DA.TITAN$TITANsummary)


nrow(nodesALL.DA[nodesALL.DA$DAsummary == "BurnDA",]) #1285
nrow(nodesALL.DA[nodesALL.DA$DAsummary == "ControlDA",]) #1757
nrow(nodesALL.DA[nodesALL.DA$DAsummary == "NotDA",]) #2767

colnames(nodesALL.DA)
unique(nodesALL.DA$NetworkRole)
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Peripheral",]) #5225
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Module Hub",]) #182
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Connector",]) #392
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Network Hub",]) #13
nrow(nodesALL.DA)

#How many hub taxa are also BurnDA?
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Module Hub" & 
                   nodesALL.DA$DAsummary == "BurnDA",]) #39
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Module Hub" & 
                   nodesALL.DA$DAsummary == "ControlDA",]) #55
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Module Hub" & 
                   nodesALL.DA$DAsummary == "NotDA",]) #88
ModHubs = c(39,55,88)

nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Network Hub" & 
                   nodesALL.DA$DAsummary == "BurnDA",]) #7
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Network Hub" & 
                   nodesALL.DA$DAsummary == "ControlDA",]) #2
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Network Hub" & 
                   nodesALL.DA$DAsummary == "NotDA",]) #4

#How many connector taxa are also BurnDA?
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Connector" & 
                   nodesALL.DA$DAsummary == "BurnDA",]) #150
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Connector" & 
                   nodesALL.DA$DAsummary == "ControlDA",]) #105
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Connector" & 
                   nodesALL.DA$DAsummary == "NotDA",]) #137

#How many peripheral taxa are also BurnDA?
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Peripheral" & 
                   nodesALL.DA$DAsummary == "BurnDA",]) #1089
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Peripheral" & 
                   nodesALL.DA$DAsummary == "ControlDA",]) #1595
nrow(nodesALL.DA[nodesALL.DA$NetworkRole == "Peripheral" & 
                   nodesALL.DA$DAsummary == "NotDA",]) #2538

#view
temp <- nodesALL.DA[nodesALL.DA$NetworkRole == "Network Hub" & 
                      nodesALL.DA$DAsummary == "BurnDA",]



ZiPiDAsum <- data.table(Group = c("BurnDA", "ControlDA", "NotDA"),
                        ModHubs = c(39,55,88),
                        NetHubs = c(7,2,4),
                        Connect = c(150,105, 137),
                        Periph = c(1089,1595,2538))
# convert values to percentages so they're comparable:
ZiPiDAsum.norm <- ZiPiDAsum[, lapply(.SD, function(x) 100*(x/sum(x))), .SDcols = c(2:5)]
ZiPiDAsum.norm <- data.table(Group = ZiPiDAsum$Group, ZiPiDAsum.norm)
ZiPiDAsum.norm.m <- melt(ZiPiDAsum.norm, id.var="Group")
head(ZiPiDAsum.norm.m)

ggplot(ZiPiDAsum.norm.m)+
  geom_bar(aes(y=value, x=variable, fill=Group), stat="identity")+
  theme_minimal()+
  scale_fill_manual(values=c("#ff0102", "#0670ff", "grey"))+
  xlab("Network Role")+
  ylab("Proportion of Taxa (%)")


##
#
## COUNT EDGES BY ASV
#
colnames(edges)

edges[666,]
edges[Taxon1 == "OTU65" & Taxon2 == "ASV100",]
edges[Taxon1 == "ASV100" & Taxon2 == "OTU65",]
#one row per couple... stack these to answer questions about edges/taxon:
edges.temp1 <- edges[,c(2:5)]
edges.temp2 <- edges[,c(1,3:5)]
names(edges.temp1)[1] <- "ASV"
names(edges.temp2)[1] <- "ASV"
#stack tables:
edges.stacked <- data.table(rbind(edges.temp1, edges.temp2))
#add relevant columns:
colnames(nodesALL.DA.TITAN)
edges.stacked.DA <- merge(edges.stacked, nodesALL.DA.TITAN[,c(1,7,46,47)], by="ASV")
head(edges.stacked.DA)
edges.stacked.DA$DAsummary2 <- ifelse(edges.stacked.DA$DAsummary == "NotDA", "NotDA", "DA")

# Does BurnDA = more edges?
nrow(edges.stacked.DA[edges.stacked.DA$DAsummary == "BurnDA", ])    #311466
nrow(edges.stacked.DA[edges.stacked.DA$DAsummary == "ControlDA", ]) #149676
nrow(edges.stacked.DA[edges.stacked.DA$DAsummary == "NotDA", ])     #380623
length(which(is.na(edges.stacked.DA$DAsummary))) #317
nrow(edges.stacked.DA) #842082
311466+149676+380623+317 #842082

# TOTAL NETWORK
nrow(edges.stacked.DA[edges.stacked.DA$Sign == "pos",]) #621397
nrow(edges.stacked.DA[edges.stacked.DA$Sign == "neg",]) #220368
621397/842082 #73.8% positive
220368/842082 #26.2% negative

#BURN DA
nrow(edges.stacked.DA[edges.stacked.DA$DAsummary == "BurnDA" & 
                        edges.stacked.DA$Sign == "pos",]) #208962
nrow(edges.stacked.DA[edges.stacked.DA$DAsummary == "BurnDA" & 
                        edges.stacked.DA$Sign == "neg",]) #102504
208962/(208962+102504) #67.1% positive
102504/(208962+102504) #32.9% negative

#CONTROL DA
nrow(edges.stacked.DA[edges.stacked.DA$DAsummary == "ControlDA" & 
                        edges.stacked.DA$Sign == "pos",]) #126374
nrow(edges.stacked.DA[edges.stacked.DA$DAsummary == "ControlDA" & 
                        edges.stacked.DA$Sign == "neg",]) #23302
126374/(126374+23302) #84.4% positive
23302/(126374+23302) #15.6% negative

# NOT DA
nrow(edges.stacked.DA[edges.stacked.DA$DAsummary == "NotDA" & 
                        edges.stacked.DA$Sign == "pos",]) #286061
nrow(edges.stacked.DA[edges.stacked.DA$DAsummary == "NotDA" & 
                        edges.stacked.DA$Sign == "neg",]) #94562
286061/(286061+94562) #75.2% positive
94562/(286061+94562) #24.8% negative

# TITAN
nrow(edges.stacked.DA[edges.stacked.DA$TITANsummary == "TITANind" & 
                        edges.stacked.DA$Sign == "pos",]) #58863
nrow(edges.stacked.DA[edges.stacked.DA$TITANsummary == "TITANind" & 
                        edges.stacked.DA$Sign == "neg",]) #26101
58863/(58863+26101) #69.3% positive
26101/(58863+26101) #30.7% negative

nrow(edges.stacked.DA) #842082
sum(208962, 102504, 126374, 23302, 286061, 94562) #841765 + 317 NAs = 842082

nrow(edges.stacked.DA[edges.stacked.DA$Sign == "neg",]) #220448
nrow(edges.stacked.DA[edges.stacked.DA$Sign == "pos",]) #621634
220448+621634 #842082

# TITAN + DA
nrow(edges.stacked.DA[edges.stacked.DA$TITANsummary == "TITANind" & 
                        edges.stacked.DA$DAsummary2 =="DA" &
                        edges.stacked.DA$Sign == "pos",]) #36053
nrow(edges.stacked.DA[edges.stacked.DA$TITANsummary == "TITANind" & 
                        edges.stacked.DA$DAsummary2 =="DA" &
                        edges.stacked.DA$Sign == "neg",]) #18215
36053/(36053+18215) #66.4.3% positive
18215/(36053+18215) #33.6% negative


# Pie chart summary
slices <- c(208962, 102504, 126374, 23302, 286061, 94562)
LABELS <- c("BurnDA-pos", "BurnDA-neg", "ControlDA-pos", "ControlDA-neg", "NotDA-pos", "NotDA-neg")
pie(slices, 
    labels = LABELS,
    main = "Network Edges Summary",
    border = FALSE,
    radius = 1,
    col = c("#ff0102", "#ff9b9b", "#009acd", "#b5ecff", "grey30", "grey80"))


# TOTAL NETWORK
621397/842082 #73.8% positive
220368/842082 #26.2% negative
#BURN DA
208962/(208962+102504) #67.1% positive
102504/(208962+102504) #32.9% negative
#CONTROL DA
126374/(126374+23302) #84.4% positive
23302/(126374+23302) #15.6% negative
# NOT DA
286061/(286061+94562) #75.2% positive
94562/(286061+94562) #24.8% negative
# TITAN
58863/(58863+26101) #69.3% positive
26101/(58863+26101) #30.7% negative
# TITAN + DA
36053/(36053+18215) #66.4.3% positive
18215/(36053+18215) #33.6% negative

PosNegEdgeSum <- data.table(Group = c("Total", "BurnDA", "ControlDA", "NotDA", "TITANind", "TITANda"),
                            Positive = c(73.8, 67.1, 84.4, 75.2, 69.3, 66.4),
                            Negative = c(26.2, 32.9, 15.6, 24.8, 30.7, 33.6))

PosNegEdgeSum.m <- melt(PosNegEdgeSum)
head(PosNegEdgeSum.m)

ggplot(PosNegEdgeSum.m)+
  geom_bar(aes(y=value, x=Group, fill=variable), stat="identity")+
  theme_minimal()+
  scale_x_discrete(limits=c("Total", "NotDA", "ControlDA", "BurnDA", "TITANind", "TITANda"),
                   labels=c("Total", "NotDA", "NegDA", "PosDA", "TITANind", "TITANda"))+ 
  scale_fill_manual(values=c("#b77eb7", "grey30"))+
  labs(fill="Correlation")+
  xlab("")+
  ylab("Proportion of Taxa (%)")




##
#
## COUNT EDGES BY INTERACTION TYPE: bacteria-fungi, fungi-fungi, bacteria-bacteria
#
colnames(edges)
colnames(nodesALL.DA.TITAN)

# merge Domain for Taxon1
edges.temp0 <- copy(edges)
head(edges.temp0)
names(edges.temp0)[1] <- "ASV" #change Taxon1 to "ASV"
head(edges.temp0)
edges.temp1 <- merge(edges.temp0, nodesALL.DA.TITAN[,c(1,7,46,47)], by="ASV") #Merge DA and TITAN for Taxon1
head(edges.temp1)
#re-name columns to facilitate the next merge:
names(edges.temp1)[6] <- "Domain1"
names(edges.temp1)[1] <- "Taxon1"
names(edges.temp1)[2] <- "ASV"
head(edges.temp1)
edges.temp2 <- merge(edges.temp1, nodesALL.DA.TITAN[,c(1,7,46,47)], by="ASV") #Merge DA and TITAN for Taxon2
head(edges.temp2)
#re-name columns back to what they should be:
names(edges.temp2)[9] <- "Domain2"
names(edges.temp2)[1] <- "Taxon2"
names(edges.temp2)[7] <- "DAsummary1"
names(edges.temp2)[8] <- "TITANsummary1"
names(edges.temp2)[10] <- "DAsummary2"
names(edges.temp2)[11] <- "TITANsummary2"
head(edges.temp2)
edges.dom <- edges.temp2[,c(2,1,3:11)] 
head(edges.dom)
unique(edges.dom$Domain2)
#remove Archaea and NA's
edges.dom <- edges.dom[Domain1 != 'Archaea']
edges.dom <- edges.dom[Domain2 != 'Archaea']
edges.dom <- edges.dom[is.na(edges.dom$Domain1) == FALSE]
edges.dom <- edges.dom[is.na(edges.dom$Domain2) == FALSE]

nrow(edges.dom[edges.dom$Domain1 == "Bacteria" & 
                 edges.dom$Domain2 == "Bacteria",]) #178656
nrow(edges.dom[edges.dom$Domain1 == "Fungi" & 
                 edges.dom$Domain2 == "Bacteria",]) #135922
nrow(edges.dom[edges.dom$Domain1 == "Bacteria" & 
                 edges.dom$Domain2 == "Fungi",]) #0
nrow(edges.dom[edges.dom$Domain1 == "Fungi" & 
                 edges.dom$Domain2 == "Fungi",]) #97111


178656+135922+97111 #411689
nrow(edges.dom) #411689

# %edges in the whole network:
135922/411689 #33% bacteria-fungi
178656/411689 #43.5% bacteria-bacteria
97111/411689 #23.5% fungi-fungi

#subset for only BurnDA taxa
head(edges.dom)
edges.dom.BurnDA <- edges.dom[DAsummary1 == "BurnDA" | DAsummary2 == "BurnDA"] #237279 rows with BurnDA in either column
edges.dom.BurnDA <- edges.dom[DAsummary1 == "BurnDA" & DAsummary2 == "BurnDA"] #69933 rows with BurnDA in both columns!
dim(edges.dom.BurnDA)

nrow(edges.dom.BurnDA[edges.dom.BurnDA$Domain1 == "Bacteria" & 
                        edges.dom.BurnDA$Domain2 == "Bacteria",]) #2971
nrow(edges.dom.BurnDA[edges.dom.BurnDA$Domain1 == "Fungi" & 
                        edges.dom.BurnDA$Domain2 == "Bacteria",]) #13694
nrow(edges.dom.BurnDA[edges.dom.BurnDA$Domain1 == "Bacteria" & 
                        edges.dom.BurnDA$Domain2 == "Fungi",]) #0
nrow(edges.dom.BurnDA[edges.dom.BurnDA$Domain1 == "Fungi" & 
                        edges.dom.BurnDA$Domain2 == "Fungi",]) #53268

# %edges in the BurnDA-only network:
13694/69933 #19.6% bacteria-fungi
2971/69933 #4.2% bacteria-bacteria
53268/69933 #76.2% fungi-fungi

edges.dom.ControlDA <- edges.dom[DAsummary1 == "ControlDA" | DAsummary2 == "ControlDA"] #123961 rows with ControlDA in either column
edges.dom.ControlDA <- edges.dom[DAsummary1 == "ControlDA" & DAsummary2 == "ControlDA"] #20133 rows with ControlDA in both columns!
dim(edges.dom.ControlDA)

nrow(edges.dom.ControlDA[edges.dom.ControlDA$Domain1 == "Bacteria" & 
                           edges.dom.ControlDA$Domain2 == "Bacteria",]) #19890
nrow(edges.dom.ControlDA[edges.dom.ControlDA$Domain1 == "Fungi" & 
                           edges.dom.ControlDA$Domain2 == "Bacteria",]) #234
nrow(edges.dom.ControlDA[edges.dom.ControlDA$Domain1 == "Bacteria" & 
                           edges.dom.ControlDA$Domain2 == "Fungi",]) #0
nrow(edges.dom.ControlDA[edges.dom.ControlDA$Domain1 == "Fungi" & 
                           edges.dom.ControlDA$Domain2 == "Fungi",]) #9

# %edges in the ControlDA-only network:
234/20133 #1.2% bacteria-fungi
19890/20133 #98.8% bacteria-bacteria
9/20133 #0.04% fungi-fungi

edges.dom.NotDA <- edges.dom[DAsummary1 == "NotDA" | DAsummary2 == "NotDA"] #283084 rows with NotDA in either column
edges.dom.NotDA <- edges.dom[DAsummary1 == "NotDA" & DAsummary2 == "NotDA"] #88988 rows with NotDA in both columns!
dim(edges.dom.NotDA)

nrow(edges.dom.NotDA[edges.dom.NotDA$Domain1 == "Bacteria" & 
                       edges.dom.NotDA$Domain2 == "Bacteria",]) #60924
nrow(edges.dom.NotDA[edges.dom.NotDA$Domain1 == "Fungi" & 
                       edges.dom.NotDA$Domain2 == "Bacteria",]) #20405
nrow(edges.dom.NotDA[edges.dom.NotDA$Domain1 == "Bacteria" & 
                       edges.dom.NotDA$Domain2 == "Fungi",]) #0
nrow(edges.dom.NotDA[edges.dom.NotDA$Domain1 == "Fungi" & 
                       edges.dom.NotDA$Domain2 == "Fungi",]) #7659

# %edges in the NotDA-only network:
60924/88988 #68.5% bacteria-fungi
20405/88988 #22.9% bacteria-bacteria
7659/88988 #8.6% fungi-fungi


#subset for only TITANinds taxa
head(edges.dom)
edges.dom.Ind <- edges.dom[TITANsummary1 == "TITANind" | TITANsummary2 == "TITANind"] #77911 rows with TITANind in either column
edges.dom.Ind <- edges.dom[TITANsummary1 == "TITANind" & TITANsummary2 == "TITANind"] #4883 rows with TITANind in both columns!
dim(edges.dom.Ind)

nrow(edges.dom.Ind[edges.dom.Ind$Domain1 == "Bacteria" & 
                     edges.dom.Ind$Domain2 == "Bacteria",]) #1233
nrow(edges.dom.Ind[edges.dom.Ind$Domain1 == "Fungi" & 
                     edges.dom.Ind$Domain2 == "Bacteria",]) #1633
nrow(edges.dom.Ind[edges.dom.Ind$Domain1 == "Bacteria" & 
                     edges.dom.Ind$Domain2 == "Fungi",]) #0
nrow(edges.dom.Ind[edges.dom.Ind$Domain1 == "Fungi" & 
                     edges.dom.Ind$Domain2 == "Fungi",]) #2017

# %edges in the TITANinds-only network:
1633/4883 #33.4% bacteria-fungi
1233/4883 #25.3% bacteria-bacteria
2017/4883 #41.3% fungi-fungi


head(edges.dom)
edges.dom$DAsummary1.2 <- ifelse(edges.dom$DAsummary1 == "NotDA", "NotDA", "DA")
edges.dom$DAsummary2.2 <- ifelse(edges.dom$DAsummary2 == "NotDA", "NotDA", "DA")

#subset for only TITANinds & DA taxa
head(edges.dom)
edges.dom.Ind.DA <- edges.dom[TITANsummary1 == "TITANind" & TITANsummary2 == "TITANind" &
                                DAsummary1.2 == "DA" & DAsummary2.2 == "DA"] #2085 rows with TITANind and DA in both columns!
dim(edges.dom.Ind.DA)

nrow(edges.dom.Ind.DA[edges.dom.Ind.DA$Domain1 == "Bacteria" & 
                        edges.dom.Ind.DA$Domain2 == "Bacteria",]) #194
nrow(edges.dom.Ind.DA[edges.dom.Ind.DA$Domain1 == "Fungi" & 
                        edges.dom.Ind.DA$Domain2 == "Bacteria",]) #488
nrow(edges.dom.Ind.DA[edges.dom.Ind.DA$Domain1 == "Bacteria" & 
                        edges.dom.Ind.DA$Domain2 == "Fungi",]) #0
nrow(edges.dom.Ind.DA[edges.dom.Ind.DA$Domain1 == "Fungi" & 
                        edges.dom.Ind.DA$Domain2 == "Fungi",]) #1403
# %edges in the TITANda-only network:
488/2085 #23.4% bacteria-fungi
194/2085 #9.3% bacteria-bacteria
1403/2085 #67.3% fungi-fungi



#
## TURN THE BELOW INTO A TABLE TO GRAPH IT!
## Stacked bar: x = (total, ControlDA, BurnDA, NotDA), y = (%b-b,%b-f, %f-f)
#
# %edges in the whole network:
178656/411689 #43.5% bacteria-bacteria
135922/411689 #33% bacteria-fungi
97111/411689 #23.5% fungi-fungi

# %edges in the ControlDA-only network:
19890/20133 #98.793% bacteria-bacteria
234/20133 #1.162% bacteria-fungi
9/20133 #0.045% fungi-fungi

# %edges in the BurnDA-only network:
2971/69933 #4.2% bacteria-bacteria
13694/69933 #19.6% bacteria-fungi
53268/69933 #76.2% fungi-fungi

# %edges in the NotDA-only network:
20405/88988 #22.9% bacteria-bacteria
60924/88988 #68.5% bacteria-fungi
7659/88988 #8.6% fungi-fungi

# %edges in the TITANinds-only network:
1233/4883 #25.3% bacteria-bacteria
1633/4883 #33.4% bacteria-fungi
2017/4883 #41.3% fungi-fungi

# %edges in the TITANda-only network:
488/2085 #23.4% bacteria-fungi
194/2085 #9.3% bacteria-bacteria
1403/2085 #67.3% fungi-fungi

CorEdgeSum <- data.table(Group = c("Total", "ControlDA", "BurnDA", "NotDA", "TITANind", "TITANda"),
                         BB = c(43.5, 98.793, 4.2, 22.9, 25.3, 9.3),
                         BF = c(33, 1.162, 19.6, 68.5, 33.4, 23.4), 
                         FF = c(23.5, 0.045, 76.2, 8.6, 41.3, 67.3))
CorEdgeSum.m <- melt(CorEdgeSum)
head(CorEdgeSum.m)

ggplot(CorEdgeSum.m)+
  geom_bar(aes(y=value, x=Group, fill=variable), stat="identity")+
  theme_minimal()+
  scale_x_discrete(limits=c("Total", "NotDA", "ControlDA", "BurnDA", "TITANind", "TITANda"),
                   labels=c("Total", "NotDA", "NegDA", "PosDA", "TITANind", "TITANda"))+ 
  scale_fill_manual(values=c("#798E87", "#C27D38", "#CBC590"), 
                    labels=c("Bacteria-Bacteria", "Bacteria-Fungi", "Fungi-Fungi"))+
  labs(fill="Interaction")+
  xlab("")+
  ylab("Proportion of Edges (%)")

### RE-MAKE THIS AS A BARPLOT
## x = total, burnDA, ControlDA, NotDA
# y = %pos and %neg (stacked bar)
# Pie chart summary
slices <- c(208962, 102504, 126374, 23302, 286061, 94562, 317)
LABELS <- c("BurnDA-pos", "BurnDA-neg", "ControlDA-pos", "ControlDA-neg", "NotDA-pos", "NotDA-neg")
pie(slices, 
    labels = LABELS,
    main = "Network Edges Summary",
    border = FALSE,
    radius = 1,
    col = c("#ff0102", "#ff9b9b", "#009acd", "#b5ecff", "grey30", "grey80"))

nrow(edges.stacked.DA) #842082
sum(208962, 102504, 126374, 23302, 286061, 94562) #841765 + 317 NAs = 842082

#remove NAs and then calculate network totals
edges.stacked.DA.na <- edges.stacked.DA[is.na(edges.stacked.DA$DAsummary) == FALSE]
nrow(edges.stacked.DA.na[edges.stacked.DA.na$Sign == "neg",]) #220368 
nrow(edges.stacked.DA.na[edges.stacked.DA.na$Sign == "pos",]) #621397
220368+621397 #841765


PosNegEdgeSum <- data.table(Group = c("Total", "ControlDA", "BurnDA", "NotDA"),
                            Pos = c(100*621397/841765, 100*126374/(126374+23302), 100*208962/(208962+102504), 100*286061/(286061+94562)),
                            Neg = c(100*220368/841765, 100*23302/(126374+23302), 100*102504/(208962+102504), 100*94562/(286061+94562)))

PosNegEdgeSum.m <- melt(PosNegEdgeSum)
head(PosNegEdgeSum.m)

ggplot(PosNegEdgeSum.m)+
  geom_bar(aes(y=value, x=Group, fill=variable), stat="identity")+
  theme_minimal()+
  scale_x_discrete(limits=c("Total", "NotDA", "ControlDA", "BurnDA"),
                   labels=c("Total", "NotDA", "Control", "Burn"))+ 
  scale_fill_manual(values=c("grey80", "grey20"), 
                    labels=c("Positive", "Negative"))+
  labs(fill="Correlation")+
  xlab("")+
  ylab("Proportion of Taxa (%)")


############################################
### sncm.ft() function from Burns, et al ###
############## Figure 6 ####################
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

############################
#### SNCM Neutral Model ####
### stats & scatterplots ###
######### Figure 6 #########
############################

# Step 1: assign the scnm.fit funtion in the previous section
#
## RUN SCNM.FIT ON COMBINED ITS & 16S DATA!
# fit the neutral model to *all* the data in each plot: 
# do not split data by time! only by plot!
# 
# add Plot information to combined ITS16S otu_table:
colnames(MetadataWithReads)
dim(ITS16Sotu) # 255 x 110139 (110138 taxa + SampID column)
dim(TAXtableITS16S) # 110138 x 8
class(ITS16Sotu)
ITS16Sotu.df <- as.data.frame(ITS16Sotu)
ITS16Sotu.df[1:5,1:5] #ASVs=Fungi
ITS16Sotu.df <- rownames_to_column(ITS16Sotu.df, var="SampID")
# ITS16Sotu table comes from generating the input table for FastSpar!
# Merge Plot column onto the OTU table by matching the SampID columns:
ITS16Sotu.sites <- merge(MetadataWithReads[,c(6,8)], ITS16Sotu.df, by="SampID")
ITS16Sotu.sites[1:5, 1:5]
#convert ASV/OTU columns to numeric (somehow they became integers...)
ITS16Sotu.sites.num <- data.frame(lapply(ITS16Sotu.sites[,3:110140], as.numeric))
ITS16Sotu.sites.num[1:5, 1:5]
# add back plot column:
ITS16Sotu.Plot <- as.data.table(cbind(ITS16Sotu.sites$Plot, ITS16Sotu.sites.num))
ITS16Sotu.Plot[1:5, 1:5]
names(ITS16Sotu.Plot)[1] <- "Plot"


# Filter Master otu_table for ASVs that we want to highlight on SCNM scatterplots
# (see final section in this script for the Master OTU table!)
PyrophileVennSection <- ITS16S_ResultsMASTER[TITANsummary == "TITANburn" & 
                                               DA_logFC > 2 &
                                               DA_adj.P.Val < 0.01,] 
colnames(PyrophileVennSection)
Dat4SCNM <- PyrophileVennSection[,c(1:3,6,71,72,79,80,81,88,89,90,97,98,99,106,111:117, 127:134)]
Dat4SCNM.400 <- Dat4SCNM[SCNMhiNonNeu == "NonNeutral" & SCNMsummary == "SCNMburn",] #13
Dat4SCNM.321E <- Dat4SCNM[SCNMloNonNeu == "NonNeutral" & SCNMsummary == "SCNMburn",] #24
Dat4SCNM.240 <- Dat4SCNM[SCNM1cNonNeu == "NonNeutral" & SCNMsummary == "SCNMburn",] #0
Dat4SCNM.321W <- Dat4SCNM[SCNM2cNonNeu == "NonNeutral" & SCNMsummary == "SCNMburn",] #0
#
##
#
## comp400! (Hi)
#
#
# subset for data only from our comp400 site:
ITS16Sotu400 <- ITS16Sotu.Plot[Plot == "400",] #67 samples
class(ITS16Sotu400$OTU25) #numeric
ITS16Sotu400[1:5, 1:5]
ITS16Sotu400[1:5, 110130:110138]
ITS16Sotu400$Plot <- NULL # remove Plot column that was only needed for subsetting
ITS16Sotu400.m <- as.matrix(ITS16Sotu400) #convert to matrix
ITS16Sotu400.m[1:5, 1:5]
# keep only columns that sum to a value greater than zero:
ITS16Sotu400.nonzero <- ITS16Sotu400.m[ ,colSums(ITS16Sotu400.m) > 0] 
ITS16Sotu400.nonzero[1:5, 1:5]
dim(ITS16Sotu400.nonzero) #67 x 24599
dim(ITS16Sotu400.m) #67 x 110138

# Fit the Neutral Model
BS400fit <- sncm.fit(spp=ITS16Sotu400.nonzero, stats=FALSE) #outputs model fit values
BS400fit.stats <- sncm.fit(spp=ITS16Sotu400.nonzero, stats=TRUE) #outputs general stats like R-squared, m, and AIC!

view(BS400fit.stats)
head(BS400fit)
tail(BS400fit)
# move rownames to a column
BS400fit <- tibble::rownames_to_column(BS400fit, "ASV")

# add a "group" column that can be used ot assign colors during plotting
# "if freq > pred.upr put "Above" in the "group" column, 
# otherwise, if freq < pred.lwr, put "Below" in the "group" column... for all else, "Neutral"
BS400fit_grouped <- BS400fit %>%
  mutate(group = ifelse(freq > pred.upr, "Above", ifelse(freq < pred.lwr,"Below","Neutral")))

colnames(BS400fit_grouped)
BS400fit_grouped.tax <- merge(BS400fit_grouped, ITS16Stax.df.rn, by="ASV", all.x=TRUE)

# There are two 16S ASVs that are simply classified as "Eukaryota"... remove them:
BS400fit_grouped.tax.sub <- BS400fit_grouped.tax[BS400fit_grouped.tax$Domain !="Eukaryota",]

# Create tables to highlight the SCNM fit of taxa-of-interest
class(ITS16S_DAtable.summary) # differential abundance restuls
ITS16S_DAtable.summary <- rownames_to_column(ITS16S_DAtable.summary, var="ASV")
ITS16S_DAtable.summary.dt <- data.table(ITS16S_DAtable.summary)
# merge SCNM results with TITAN and/or DA tables:
BS400fit_TITAN <- merge(nodesALL.DA.TITAN[TITAN400 == "TITAN400",], BS400fit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID", all.x=TRUE)

# Plot:
ggplot(BS400fit_grouped.tax.sub)+
  geom_point(aes(x = log(p), y = freq, color=group, shape=Domain), size=1, stroke=0.7)+
  #circle ASVs-of-interest:
  geom_point(data=Dat4SCNM, shape=1, stroke=1, size=1.5, aes(x = log(SCNMhi_p), y = SCNMhi_freq), color="#20f038")+ 
  guides(color = guide_legend(override.aes = list(size=2)))+
  guides(shape = guide_legend(override.aes = list(size=2)))+
  geom_line(aes(x = log(p), y = freq.pred), color="grey30")+
  geom_line(aes(x = log(p), y = pred.lwr), color="grey30", linetype=2)+
  geom_line(aes(x = log(p), y = pred.upr), color="grey30", linetype=2)+
  theme_classic()+
  scale_color_manual(values= c("orange", "purple", "grey"))+
  scale_shape_manual(values=c(17,15,4))+
  ggtitle("Hi-Burn Neutral Model fit, all time-points, ITS & 16S")+
  xlab("log(Mean Relative Abundance)")+
  ylab("Frequency")+
  annotate(geom="text", x=min(log(BS400fit_grouped.tax.sub$p)), y=max(BS400fit_grouped.tax.sub$freq), 
           label=paste0("R-squared = ", formatC(BS400fit.stats$Rsqr, digits=4)), 
           size=4, vjust = "inward", hjust = "inward")+
  annotate(geom="text", x=min(log(BS400fit_grouped.tax.sub$p)), y=0.9, 
           label=paste0("m = ", formatC(BS400fit.stats$m, digits=2)), # "digits" are non-zero numbers... digits=2: 0.0031, digits=4: 0.003192
           size=4, vjust = "inward", hjust = "inward")
#export PDF: 7x3 

# 
##
#
## comp321E!
#
# Subset otu_table:
ITS16Sotu321E <- ITS16Sotu.Plot[Plot == "321east"] #94 samples
ITS16Sotu321E$Plot <- NULL
ITS16Sotu321E[1:5, 1:5]
ITS16Sotu321E[1:5, 110130:110138]
ITS16Sotu321E.m <- as.matrix(ITS16Sotu321E)
ITS16Sotu321E.m[1:5, 1:5]
# keep only columns that sum to a value greater than zero:
ITS16Sotu321E.nonzero <- ITS16Sotu321E.m[ ,colSums(ITS16Sotu321E.m) > 0]
dim(ITS16Sotu321E.nonzero) #30071 taxa x 94 samples
dim(ITS16Sotu321E.m) #110138 taxa
ITS16Sotu321E.nonzero[1:5, 1:5]
ITS16Sotu321E.nonzero[1:5, 30065:30071]

# Fit the Neutral Model
BS321Efit <- sncm.fit(spp=ITS16Sotu321E.nonzero, stats=FALSE) #outputs model fit values
BS321Efit.stats <- sncm.fit(spp=ITS16Sotu321E.nonzero, stats=TRUE) #outputs general stats like R-squared, m, and AIC!

view(BS321Efit.stats)
head(BS321Efit)

# move rownames to a column
BS321Efit <- tibble::rownames_to_column(BS321Efit, "ASV")
# add a "group" column that can be used ot assign colors during plotting
# "if freq > pred.upr color the dot red, otherwise, if freq < pred.lwr, color the point yellow, otherwise black"
BS321Efit_grouped <- BS321Efit %>%
  mutate(group = ifelse(freq > pred.upr, "Above", ifelse(freq<pred.lwr,"Below","Neutral")))

BS321Efit_grouped.tax <- merge(BS321Efit_grouped, ITS16Stax.df.rn, by="ASV", all.x=TRUE)

#plot:
ggplot(BS321Efit_grouped.tax)+
  geom_point(aes(x = log(p), y = freq, color=group, shape=Domain), size=1, stroke=0.7)+
  #circle ASVs-of-interest:
  geom_point(data=Dat4SCNM, shape=1, stroke=1, size=1.5, aes(x = log(SCNMlo_p), y = SCNMlo_freq), color="#5aff00")+  
  guides(color = guide_legend(override.aes = list(size=2)))+
  guides(shape = guide_legend(override.aes = list(size=2)))+
  geom_line(aes(x = log(p), y = freq.pred), color="grey30")+
  geom_line(aes(x = log(p), y = pred.lwr), color="grey30", linetype=2)+
  geom_line(aes(x = log(p), y = pred.upr), color="grey30", linetype=2)+
  theme_classic()+
  scale_color_manual(values= c("orange", "purple", "grey"))+
  scale_shape_manual(values=c(17,15,4))+
  ggtitle("Lo-Burn Neutral Model fit, all time-points, ITS & 16S")+
  xlab("log(Mean Relative Abundance)")+
  ylab("Frequency")+
  annotate(geom="text", x=min(log(BS321Efit_grouped.tax$p)), y=max(BS321Efit_grouped.tax$freq), 
           label=paste0("R-squared = ", formatC(BS321Efit.stats$Rsqr, digits=4)), 
           size=4, vjust = "inward", hjust = "inward")+
  annotate(geom="text", x=min(log(BS321Efit_grouped.tax$p)), y=0.9, 
           label=paste0("m = ", formatC(BS321Efit.stats$m, digits=2)), 
           size=4, vjust = "inward", hjust = "inward")
#export PDF: 7x3 


#
##
#
## comp321W!
#
# subset otu_table:
ITS16Sotu321W <- ITS16Sotu.Plot[Plot == "321west"] #46 samples
ITS16Sotu321W$Plot <- NULL
ITS16Sotu321W[1:5, 1:5]
ITS16Sotu321W[1:5, 110130:110138]
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

view(BS321Wfit.stats)
head(BS321Wfit)

# move rownames to a column
BS321Wfit <- tibble::rownames_to_column(BS321Wfit, "ASV")
# add a "group" column that can be used ot assign colors during plotting
# "if freq > pred.upr color the dot red, otherwise, if freq < pred.lwr, color the point yellow, otherwise black"
BS321Wfit_grouped <- BS321Wfit %>%
  mutate(group = ifelse(freq > pred.upr, "Above", ifelse(freq<pred.lwr,"Below","Neutral")))

BS321Wfit_grouped.tax <- merge(BS321Wfit_grouped, ITS16Stax.df.rn, by="ASV", all.x=TRUE)

#plot
ggplot(BS321Wfit_grouped.tax)+
  geom_point(aes(x = log(p), y = freq, color=group, shape=Domain), size=1, stroke=0.7)+
  #circle ASVs-of-interest:
  geom_point(data=Dat4SCNM, shape=1, stroke=1, size=1.5, aes(x = log(SCNM2c_p), y = SCNM2c_freq), color="#20f038")+ 
  guides(color = guide_legend(override.aes = list(size=2)))+
  guides(shape = guide_legend(override.aes = list(size=2)))+
  geom_line(aes(x = log(p), y = freq.pred), color="grey30")+
  geom_line(aes(x = log(p), y = pred.lwr), color="grey30", linetype=2)+
  geom_line(aes(x = log(p), y = pred.upr), color="grey30", linetype=2)+
  theme_classic()+
  scale_color_manual(values= c("orange", "purple", "grey"))+
  scale_shape_manual(values=c(17,15,4))+
  ggtitle("Control2 Neutral Model fit, all time-points, ITS & 16S")+
  xlab("log(Mean Relative Abundance)")+
  ylab("Frequency")+
  annotate(geom="text", x=min(log(BS321Wfit_grouped.tax.sub$p)), y=max(BS321Wfit_grouped.tax$freq), 
           label=paste0("R-squared = ", formatC(BS321Wfit.stats$Rsqr, digits=4)), 
           size=4, vjust = "inward", hjust = "inward")+
  annotate(geom="text", x=min(log(BS321Wfit_grouped.tax.sub$p)), y=0.9, 
           label=paste0("m = ", formatC(BS321Wfit.stats$m, digits=2)), 
           size=4, vjust = "inward", hjust = "inward")
#export PDF: 7x3 


#
##
#
## comp240!
#
# subset otu_table:
ITS16Sotu240 <- ITS16Sotu.Plot[Plot == "240"] #48 samples
ITS16Sotu240$Plot <- NULL
ITS16Sotu240[1:5, 1:5]
ITS16Sotu240[1:5, 110130:110138]
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

view(BS240fit.stats)
head(BS240fit)

# move rownames to a column
BS240fit <- tibble::rownames_to_column(BS240fit, "ASV")
# add a "group" column that can be used ot assign colors during plotting
# "if freq > pred.upr color the dot red, otherwise, if freq < pred.lwr, color the point yellow, otherwise black"
BS240fit_grouped <- BS240fit %>%
  mutate(group = ifelse(freq > pred.upr, "Above", ifelse(freq<pred.lwr,"Below","Neutral")))

BS240fit_grouped.tax <- merge(BS240fit_grouped, ITS16Stax.df.rn, by="ASV", all.x=TRUE)

#plot:
ggplot(BS240fit_grouped.tax)+
  geom_point(aes(x = log(p), y = freq, color=group, shape=Domain), size=1, stroke=0.7)+
  #circle ASVs-of-interest:
  geom_point(data=Dat4SCNM, shape=1, stroke=1, size=1.5, aes(x = log(SCNM1c_p), y = SCNM1c_freq), color="#20f038")+  
  guides(color = guide_legend(override.aes = list(size=2)))+
  guides(shape = guide_legend(override.aes = list(size=2)))+
  geom_line(aes(x = log(p), y = freq.pred), color="grey30")+
  geom_line(aes(x = log(p), y = pred.lwr), color="grey30", linetype=2)+
  geom_line(aes(x = log(p), y = pred.upr), color="grey30", linetype=2)+
  theme_classic()+
  scale_color_manual(values= c("orange", "purple", "grey"))+
  scale_shape_manual(values=c(17,15,4))+
  ggtitle("Control1 Neutral Model fit, all time-points, ITS & 16S")+
  xlab("log(Mean Relative Abundance)")+
  ylab("Frequency")+
  annotate(geom="text", x=min(log(BS240fit_grouped.tax$p)), y=max(BS240fit_grouped.tax$freq), 
           label=paste0("R-squared = ", formatC(BS240fit.stats$Rsqr, digits=4)), 
           size=4, vjust = "inward", hjust = "inward")+
  annotate(geom="text", x=min(log(BS240fit_grouped.tax$p)), y=0.9, 
           label=paste0("m = ", formatC(BS240fit.stats$m, digits=2)), 
           size=4, vjust = "inward", hjust = "inward")
#export PDF: 7x3 

#
#
####
#
###
#
##
#

######################################################
##### Stacked Bar Plots of SCNCM fits (Figure 6A) ####
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


## Barplots connecting Network Modules and SCNM!
#
# use the table() function to count how many of each varibale is in the group column (above, below, or neutral)
ModuleTaxa <- fread(file.choose())
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


#
##
#
##
#
##
#



## Barplots connecting DIFFERENTIAL ABUNDANCE and SCNM!
head(ITS16S_DAtable.summary)
#
# use the table() function to count how many of each varibale is in the group column (above, below, or neutral)
colnames(nodesALL.DA.TITAN)
colnames(BS400fit_grouped.tax.sub)

## POSITIVE DA
# merge DAsummary column onto SCNM results for each plot:
BS400fit_Bda <- merge(nodesALL.DA.TITAN[DAsummary == "PosDA", c(1,46)], BS400fit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")
BS240fit_Bda <- merge(nodesALL.DA.TITAN[DAsummary == "PosDA", c(1,46)], BS240fit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")
BS321Wfit_Bda <- merge(nodesALL.DA.TITAN[DAsummary == "PosDA", c(1,46)], BS321Wfit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")
BS321Efit_Bda <- merge(nodesALL.DA.TITAN[DAsummary == "PosDA", c(1,46)], BS321Efit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")

# sum up the number of rows (taxa) by group (above, below, or neutral):
BS400fit_Bda.sub <- as.data.frame(table(BS400fit_Bda$group)) 
BS400fit_Bda.sub$Plot <- "hi" #column designating the plot
BS400fit_Bda.sub.prop <- mutate_if(BS400fit_Bda.sub, is.numeric, funs(./sum(.)*100)) #convert row count to percent

BS240fit_Bda.sub <- as.data.frame(table(BS240fit_Bda$group)) 
BS240fit_Bda.sub$Plot <- "1c"
BS240fit_Bda.sub.prop <- mutate_if(BS240fit_Bda.sub, is.numeric, funs(./sum(.)*100)) 

BS321Wfit_Bda.sub <- as.data.frame(table(BS321Wfit_Bda$group)) 
BS321Wfit_Bda.sub$Plot <- "2c"
BS321Wfit_Bda.sub.prop <- mutate_if(BS321Wfit_Bda.sub, is.numeric, funs(./sum(.)*100)) 

BS321Efit_Bda.sub <- as.data.frame(table(BS321Efit_Bda$group)) 
BS321Efit_Bda.sub$Plot <- "lo"
BS321Efit_Bda.sub.prop <- mutate_if(BS321Efit_Bda.sub, is.numeric, funs(./sum(.)*100)) 

fitBDAprop <- as.data.frame(rbind(BS321Efit_Bda.sub.prop, BS321Wfit_Bda.sub.prop, BS240fit_Bda.sub.prop, BS400fit_Bda.sub.prop))
fitBDAprop$Plot <- factor(fitBDAprop$Plot, levels=c("1c", "2c", "lo", "hi"))

ggplot(fitBDAprop)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=30))+
  theme(axis.title = element_text(color="black", size=30))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("PosDA")+
  ylab("Frequency (%)")+
  labs(fill="Group")

#
##
#
## NEGATIVE DA
# merge DAsummary column onto SCNM results for each plot:
BS400fit_Cda <- merge(nodesALL.DA.TITAN[DAsummary == "NegDA", c(1,46)], BS400fit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")
BS240fit_Cda <- merge(nodesALL.DA.TITAN[DAsummary == "NegDA", c(1,46)], BS240fit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")
BS321Wfit_Cda <- merge(nodesALL.DA.TITAN[DAsummary == "NegDA", c(1,46)], BS321Wfit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")
BS321Efit_Cda <- merge(nodesALL.DA.TITAN[DAsummary == "NegDA", c(1,46)], BS321Efit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")

# sum up the number of rows (taxa) by group (above, below, or neutral):
BS400fit_Cda.sub <- as.data.frame(table(BS400fit_Cda$group)) 
BS400fit_Cda.sub$Plot <- "hi" #column designating the plot
BS400fit_Cda.sub.prop <- mutate_if(BS400fit_Cda.sub, is.numeric, funs(./sum(.)*100)) #convert row count to percent

BS240fit_Cda.sub <- as.data.frame(table(BS240fit_Cda$group)) 
BS240fit_Cda.sub$Plot <- "1c"
BS240fit_Cda.sub.prop <- mutate_if(BS240fit_Cda.sub, is.numeric, funs(./sum(.)*100)) 

BS321Wfit_Cda.sub <- as.data.frame(table(BS321Wfit_Cda$group)) 
BS321Wfit_Cda.sub$Plot <- "2c"
BS321Wfit_Cda.sub.prop <- mutate_if(BS321Wfit_Cda.sub, is.numeric, funs(./sum(.)*100)) 

BS321Efit_Cda.sub <- as.data.frame(table(BS321Efit_Cda$group)) 
BS321Efit_Cda.sub$Plot <- "lo"
BS321Efit_Cda.sub.prop <- mutate_if(BS321Efit_Cda.sub, is.numeric, funs(./sum(.)*100)) 

fitCDAprop <- as.data.frame(rbind(BS321Efit_Cda.sub.prop, BS321Wfit_Cda.sub.prop, BS240fit_Cda.sub.prop, BS400fit_Cda.sub.prop))
fitCDAprop$Plot <- factor(fitCDAprop$Plot, levels=c("1c", "2c", "lo", "hi"))

ggplot(fitCDAprop)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=30))+
  theme(axis.title = element_text(color="black", size=30))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("NegDA")+
  ylab("Frequency (%)")+
  labs(fill="Group")

#
##
#
## NOT DA
# merge DAsummary column onto SCNM results for each plot:
BS400fit_Nda <- merge(nodesALL.DA.TITAN[DAsummary == "NotDA", c(1,46)], BS400fit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")
BS240fit_Nda <- merge(nodesALL.DA.TITAN[DAsummary == "NotDA", c(1,46)], BS240fit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")
BS321Wfit_Nda <- merge(nodesALL.DA.TITAN[DAsummary == "NotDA", c(1,46)], BS321Wfit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")
BS321Efit_Nda <- merge(nodesALL.DA.TITAN[DAsummary == "NotDA", c(1,46)], BS321Efit_grouped.tax.sub, by.x="ASV", by.y="OTU_ID")

# sum up the number of rows (taxa) by group (above, below, or neutral):
BS400fit_Nda.sub <- as.data.frame(table(BS400fit_Nda$group)) 
BS400fit_Nda.sub$Plot <- "hi" #column designating the plot
BS400fit_Nda.sub.prop <- mutate_if(BS400fit_Nda.sub, is.numeric, funs(./sum(.)*100)) #convert row count to percent

BS240fit_Nda.sub <- as.data.frame(table(BS240fit_Nda$group)) 
BS240fit_Nda.sub$Plot <- "1c"
BS240fit_Nda.sub.prop <- mutate_if(BS240fit_Nda.sub, is.numeric, funs(./sum(.)*100)) 

BS321Wfit_Nda.sub <- as.data.frame(table(BS321Wfit_Nda$group)) 
BS321Wfit_Nda.sub$Plot <- "2c"
BS321Wfit_Nda.sub.prop <- mutate_if(BS321Wfit_Nda.sub, is.numeric, funs(./sum(.)*100)) 

BS321Efit_Nda.sub <- as.data.frame(table(BS321Efit_Nda$group)) 
BS321Efit_Nda.sub$Plot <- "lo"
BS321Efit_Nda.sub.prop <- mutate_if(BS321Efit_Nda.sub, is.numeric, funs(./sum(.)*100)) 

fitNDAprop <- as.data.frame(rbind(BS321Efit_Nda.sub.prop, BS321Wfit_Nda.sub.prop, BS240fit_Nda.sub.prop, BS400fit_Nda.sub.prop))
fitNDAprop$Plot <- factor(fitNDAprop$Plot, levels=c("1c", "2c", "lo", "hi"))

ggplot(fitNDAprop)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=30))+
  theme(axis.title = element_text(color="black", size=30))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("NotDA")+
  ylab("Frequency (%)")+
  labs(fill="Group")



Dat4SCNMHi.sub <- as.data.frame(table(Dat4SCNM$SCNMhi_group)) 
Dat4SCNMHi.sub$Plot <- "Hi"
Dat4SCNMHi.sub.prop <- mutate_if(Dat4SCNMHi.sub, is.numeric, funs(./sum(.)*100)) 

Dat4SCNMlo.sub <- as.data.frame(table(Dat4SCNM$SCNMlo_group)) 
Dat4SCNMlo.sub$Plot <- "Lo"
Dat4SCNMlo.sub.prop <- mutate_if(Dat4SCNMlo.sub, is.numeric, funs(./sum(.)*100)) 

Dat4SCNM1c.sub <- as.data.frame(table(Dat4SCNM$SCNM1c_group)) 
Dat4SCNM1c.sub$Plot <- "1c"
Dat4SCNM1c.sub.prop <- mutate_if(Dat4SCNM1c.sub, is.numeric, funs(./sum(.)*100)) 

Dat4SCNM2c.sub <- as.data.frame(table(Dat4SCNM$SCNM2c_group)) 
Dat4SCNM2c.sub$Plot <- "2c"
Dat4SCNM2c.sub.prop <- mutate_if(Dat4SCNM2c.sub, is.numeric, funs(./sum(.)*100)) 


fitTDAprop <- as.data.frame(rbind(Dat4SCNMHi.sub.prop, Dat4SCNMlo.sub.prop, Dat4SCNM1c.sub.prop, Dat4SCNM2c.sub.prop))
fitTDAprop$Plot <- factor(fitTDAprop$Plot, levels=c("1c", "2c", "Lo", "Hi"))

ggplot(fitTDAprop)+
  geom_bar(aes(x=Plot, y=Freq, fill=Var1), stat="identity")+
  theme_classic()+
  theme(axis.text = element_text(color="black", size=30))+
  theme(axis.title = element_text(color="black", size=30))+
  theme(legend.text = element_text(size=14))+
  theme(legend.title = element_text(size=14))+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= c("orange", "purple", "grey"))+
  ggtitle("TITAN+DA")+
  ylab("Frequency (%)")+
  labs(fill="Group")




#################################################
######### Pyrophile Abundance Over Time #########
################## Figure 7 #####################
#################################################
#
PyrophileVennSection <- ITS16S_ResultsMASTER[TITANsummary == "TITANburn" & 
                                               SCNMsummary == "SCNMburn" &
                                               DA_logFC > 2 &
                                               DA_adj.P.Val < 0.01,] 

## plot average abundance of ASVs that passed TITANburn+SCNMburn+PosDA filter
## 29 ASVs total --> 19 classified to genus, representing 15 unique genera
## the remaining 11 ASVs are all Fungi = 
## 3 Basidiomycota (Agaricomycetes), 7 Ascomycota (all different), and 1 just "Fungi"


#Subset OTU table for just the ASVs that passed the filter
Pyrophile_OTUtable <- ITS16Sotu.df[ , PyrophileVennSection[!is.na(PyrophileVennSection$Genus), ASV] ]
Pyrophile_OTUtable$SampID <- ITS16Sotu.df$SampID
#transpose and add taxonomy:
Pyrophile_OTUtable.t <- as.data.table(t(Pyrophile_OTUtable))
names(Pyrophile_OTUtable.t) <- as.character(Pyrophile_OTUtable.t[19, ]) #copy row19 into the colnames spot
Pyrophile_OTUtable.t <- Pyrophile_OTUtable.t[1:18, ASV:=colnames(Pyrophile_OTUtable[,1:18])] #remove redundant row19 and add ASV column
Pyrophile_OTUtable.t <- Pyrophile_OTUtable.t[,c(256,1:255)] #move ASV column up front
Pyrophile_OTUtable.tax <- as.data.frame(merge(Pyrophile_OTUtable.t, ITS16Stax.df.rn[,c(1,7)], by="ASV", all.x=TRUE)) #add Genus and convert to data.table for next step
Pyrophile_OTUtable.tax$Genus

#Subset table for comp240
#subset for OTUtable columns that match comp240 SampIDs, keep Genus columnm but not ASV column...
BS240.PyroTaxOTUtable <- Pyrophile_OTUtable.tax[ , c("Genus", ITS16S.MetaNeemonika[Plot==240, SampID]) ]
dim(BS240.PyroTaxOTUtable)

#melt table to make it easier to feed ggplot
BS240.PyroTaxOTUtable.m <- reshape2::melt(BS240.PyroTaxOTUtable, id.vars="Genus") 
head(BS240.PyroTaxOTUtable.m)
# merge numdates onto table! and convert to data.table for next step
BS240.PyroTaxOTUtable.m <- as.data.table(merge(BS240.PyroTaxOTUtable.m, ITS16S.MetaNeemonika[Plot==240, c("SampID", "numdate")], 
                                               by.x="variable", by.y="SampID", all=TRUE))

BS240.PyroTaxOTUtable.m$value <- as.numeric(BS240.PyroTaxOTUtable.m$value)
head(BS240.PyroTaxOTUtable.m)

#calculate means
BS240.means <- BS240.PyroTaxOTUtable.m[ ,mean(value), by=c("Genus", "numdate"), ]
colnames(BS240.means)[3] <- "mean"


# print the max value for each genus in the burned plots
# (each genus reached it's maximum abundnace in one of these two plots)
BS400.means[ , max(mean), by = Genus]
BS321E.means[ , max(mean), by = Genus]
#normalize each ASV to the max value for it's genus:
BS240.means[grepl("Hebeloma", Genus), NormVal:=mean/13.6666667]
BS240.means[grepl("Calyptrozyma", Genus), NormVal:=mean/557.166667]
BS240.means[grepl("Botryobasidium", Genus), NormVal:=mean/16.166667]
BS240.means[grepl("Pyronema", Genus), NormVal:=mean/371.6666667]
BS240.means[grepl("Geoglossum", Genus), NormVal:=mean/181.666667]
BS240.means[grepl("Geopyxis", Genus), NormVal:=mean/422.400000]
BS240.means[grepl("Geminibasidium", Genus), NormVal:=mean/176.250000]
BS240.means[grepl("Xerocomellus", Genus), NormVal:=mean/285.583333]
BS240.means[grepl("Myxomphalia", Genus), NormVal:=mean/508.800000]
BS240.means[grepl("Solicoccozyma", Genus), NormVal:=mean/2793.200000]
BS240.means[grepl("Spizellomyces", Genus), NormVal:=mean/140.800000]
BS240.means[grepl("Hyaloscypha", Genus), NormVal:=mean/58.400000] 
BS240.means[grepl("Flavobacterium", Genus), NormVal:=mean/54.166667] 
BS240.means[grepl("Gemmatimonas", Genus), NormVal:=mean/170.833333] 
BS240.means[grepl("Massilia", Genus), NormVal:=mean/49.333333] 

# line plot of Blodgett data:
ggplot(BS240.means, aes(x=numdate, y=NormVal, color=Genus))+
  geom_point(shape=16, size=2)+
  geom_line(size=0.5)+
  scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0, 1, 2))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#1F77B4", "#3E6735", "#9467BD", "#2CA02C",  "#F0027F",
                                         "#FDB863", "#202020", "#FF7F0E", "#909090", "#8C564B",
                                         "#B2182B", "#17BECF","#98DF8A", "#AEC7E8", "#DED42F"
  ))+
  theme(panel.background = element_blank())+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90,  hjust=0, vjust=0.5))+
  xlab("Time")+
  ylab("Normalized Average Abundance")+
  ggtitle("Control #1 (ITS & 16S)")
#export PDF = 8x4

#
##
#
##
#Subset table for comp400
#subset for OTUtable columns that match comp240 SampIDs, keep Genus columnm but not ASV column...
BS400.PyroTaxOTUtable <- Pyrophile_OTUtable.tax[ , c("Genus", ITS16S.MetaNeemonika[Plot==400, SampID]) ]
dim(BS400.PyroTaxOTUtable)

#melt table to make it easier to feed ggplot
BS400.PyroTaxOTUtable.m <- reshape2::melt(BS400.PyroTaxOTUtable, id.vars="Genus") 
head(BS400.PyroTaxOTUtable.m)
# merge numdates onto table! and convert to data.table for next step
BS400.PyroTaxOTUtable.m <- as.data.table(merge(BS400.PyroTaxOTUtable.m, ITS16S.MetaNeemonika[Plot==400, c("SampID", "numdate")], 
                                               by.x="variable", by.y="SampID", all=TRUE))

BS400.PyroTaxOTUtable.m$value <- as.numeric(BS400.PyroTaxOTUtable.m$value)

#calculate means
BS400.means <- BS400.PyroTaxOTUtable.m[ ,mean(value), by=c("Genus", "numdate"), ]
colnames(BS400.means)[3] <- "mean"
head(BS400.means)

#normalize each ASV to the max value for it's genus:
# print the max value for each genus:
BS400.means[ , max(mean), by = Genus]
BS321E.means[ , max(mean), by = Genus]
BS400.means[grepl("Hebeloma", Genus), NormVal:=mean/13.6666667]
BS400.means[grepl("Calyptrozyma", Genus), NormVal:=mean/557.166667]
BS400.means[grepl("Botryobasidium", Genus), NormVal:=mean/16.166667]
BS400.means[grepl("Pyronema", Genus), NormVal:=mean/371.6666667]
BS400.means[grepl("Geoglossum", Genus), NormVal:=mean/181.666667]
BS400.means[grepl("Geopyxis", Genus), NormVal:=mean/422.400000]
BS400.means[grepl("Geminibasidium", Genus), NormVal:=mean/176.250000]
BS400.means[grepl("Xerocomellus", Genus), NormVal:=mean/285.583333]
BS400.means[grepl("Myxomphalia", Genus), NormVal:=mean/508.800000]
BS400.means[grepl("Solicoccozyma", Genus), NormVal:=mean/2793.200000]
BS400.means[grepl("Spizellomyces", Genus), NormVal:=mean/140.800000]
BS400.means[grepl("Hyaloscypha", Genus), NormVal:=mean/58.400000] 
BS400.means[grepl("Flavobacterium", Genus), NormVal:=mean/54.166667] 
BS400.means[grepl("Gemmatimonas", Genus), NormVal:=mean/170.833333] 
BS400.means[grepl("Massilia", Genus), NormVal:=mean/49.333333] 

# line plot of Blodgett data:
ggplot(BS400.means, aes(x=numdate, y=NormVal, color=Genus))+
  geom_point(shape=16, size=2)+
  geom_line(size=0.5)+
  scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0, 1, 2))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#1F77B4", "#3E6735", "#9467BD", "#2CA02C",  "#F0027F",
                                         "#FDB863", "#202020", "#FF7F0E", "#909090", "#8C564B",
                                         "#B2182B", "#17BECF","#98DF8A", "#AEC7E8", "#DED42F"
  ))+
  theme(panel.background = element_blank())+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90,  hjust=0, vjust=0.5))+
  xlab("Time")+
  ylab("Normalized Average Abundance")+
  ggtitle("Hi Plot (ITS & 16S)")
#export PDF = 8x4




#
##
#
##
#Subset table for comp321e
#subset for OTUtable columns that match comp240 SampIDs, keep Genus columnm but not ASV column...
BS321E.PyroTaxOTUtable <- Pyrophile_OTUtable.tax[ , c("Genus", ITS16S.MetaNeemonika[Plot=="321east", SampID]) ]
dim(BS321E.PyroTaxOTUtable)

#melt table to make it easier to feed ggplot
BS321E.PyroTaxOTUtable.m <- reshape2::melt(BS321E.PyroTaxOTUtable, id.vars="Genus") 
head(BS321E.PyroTaxOTUtable.m)
# merge numdates onto table! and convert to data.table for next step
BS321E.PyroTaxOTUtable.m <- as.data.table(merge(BS321E.PyroTaxOTUtable.m, ITS16S.MetaNeemonika[Plot=="321east", c("SampID", "numdate")], 
                                                by.x="variable", by.y="SampID", all=TRUE))

BS321E.PyroTaxOTUtable.m$value <- as.numeric(BS321E.PyroTaxOTUtable.m$value)

#calculate means
BS321E.means <- BS321E.PyroTaxOTUtable.m[ ,mean(value), by=c("Genus", "numdate"), ]
colnames(BS321E.means)[3] <- "mean"
head(BS321E.means)

#normalize each ASV to the max value for it's genus:
# print the max value for each genus:
BS400.means[ , max(mean), by = Genus]
BS321E.means[ , max(mean), by = Genus]
BS321E.means[grepl("Hebeloma", Genus), NormVal:=mean/13.6666667]
BS321E.means[grepl("Calyptrozyma", Genus), NormVal:=mean/557.166667]
BS321E.means[grepl("Botryobasidium", Genus), NormVal:=mean/16.166667]
BS321E.means[grepl("Pyronema", Genus), NormVal:=mean/371.6666667]
BS321E.means[grepl("Geoglossum", Genus), NormVal:=mean/181.666667]
BS321E.means[grepl("Geopyxis", Genus), NormVal:=mean/422.400000]
BS321E.means[grepl("Geminibasidium", Genus), NormVal:=mean/176.250000]
BS321E.means[grepl("Xerocomellus", Genus), NormVal:=mean/285.583333]
BS321E.means[grepl("Myxomphalia", Genus), NormVal:=mean/508.800000]
BS321E.means[grepl("Solicoccozyma", Genus), NormVal:=mean/2793.200000]
BS321E.means[grepl("Spizellomyces", Genus), NormVal:=mean/140.800000]
BS321E.means[grepl("Hyaloscypha", Genus), NormVal:=mean/58.400000] 
BS321E.means[grepl("Flavobacterium", Genus), NormVal:=mean/54.166667] 
BS321E.means[grepl("Gemmatimonas", Genus), NormVal:=mean/170.833333] 
BS321E.means[grepl("Massilia", Genus), NormVal:=mean/49.333333] 

# line plot of Blodgett data:
ggplot(BS321E.means, aes(x=numdate, y=NormVal, color=Genus))+
  geom_point(shape=16, size=2)+
  geom_line(size=0.5)+
  scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0, 1, 2))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#1F77B4", "#3E6735", "#9467BD", "#2CA02C",  "#F0027F",
                                         "#FDB863", "#202020", "#FF7F0E", "#909090", "#8C564B",
                                         "#B2182B", "#17BECF","#98DF8A", "#AEC7E8", "#DED42F"
  ))+
  theme(panel.background = element_blank())+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90,  hjust=0, vjust=0.5))+
  xlab("Time")+
  ylab("Normalized Average Abundance")+
  ggtitle("Lo plot (ITS & 16S)")
#export PDF = 8x4

#
##
#
##
#Subset table for comp321W
#subset for OTUtable columns that match comp240 SampIDs, keep Genus columnm but not ASV column...
BS321W.PyroTaxOTUtable <- Pyrophile_OTUtable.tax[ , c("Genus", ITS16S.MetaNeemonika[Plot=="321west", SampID]) ]
dim(BS321W.PyroTaxOTUtable)

#melt table to make it easier to feed ggplot
BS321W.PyroTaxOTUtable.m <- reshape2::melt(BS321W.PyroTaxOTUtable, id.vars="Genus") 
head(BS321W.PyroTaxOTUtable.m)
# merge numdates onto table! and convert to data.table for next step
BS321W.PyroTaxOTUtable.m <- as.data.table(merge(BS321W.PyroTaxOTUtable.m, ITS16S.MetaNeemonika[Plot=="321west", c("SampID", "numdate")], 
                                                by.x="variable", by.y="SampID", all=TRUE))

BS321W.PyroTaxOTUtable.m$value <- as.numeric(BS321W.PyroTaxOTUtable.m$value)

#calculate means
BS321W.means <- BS321W.PyroTaxOTUtable.m[ ,mean(value), by=c("Genus", "numdate"), ]
colnames(BS321W.means)[3] <- "mean"
head(BS321W.means)

#normalize each ASV to the max value for it's genus:
# print the max value for each genus:
BS400.means[ , max(mean), by = Genus]
BS321E.means[ , max(mean), by = Genus]
BS321W.means[grepl("Hebeloma", Genus), NormVal:=mean/13.6666667]
BS321W.means[grepl("Calyptrozyma", Genus), NormVal:=mean/557.166667]
BS321W.means[grepl("Botryobasidium", Genus), NormVal:=mean/16.166667]
BS321W.means[grepl("Pyronema", Genus), NormVal:=mean/371.6666667]
BS321W.means[grepl("Geoglossum", Genus), NormVal:=mean/181.666667]
BS321W.means[grepl("Geopyxis", Genus), NormVal:=mean/422.400000]
BS321W.means[grepl("Geminibasidium", Genus), NormVal:=mean/176.250000]
BS321W.means[grepl("Xerocomellus", Genus), NormVal:=mean/285.583333]
BS321W.means[grepl("Myxomphalia", Genus), NormVal:=mean/508.800000]
BS321W.means[grepl("Solicoccozyma", Genus), NormVal:=mean/2793.200000]
BS321W.means[grepl("Spizellomyces", Genus), NormVal:=mean/140.800000]
BS321W.means[grepl("Hyaloscypha", Genus), NormVal:=mean/58.400000] 
BS321W.means[grepl("Flavobacterium", Genus), NormVal:=mean/54.166667] 
BS321W.means[grepl("Gemmatimonas", Genus), NormVal:=mean/170.833333] 
BS321W.means[grepl("Massilia", Genus), NormVal:=mean/49.333333] 

# line plot of Blodgett data:
ggplot(BS321W.means, aes(x=numdate, y=NormVal, color=Genus))+
  geom_point(shape=16, size=2)+
  geom_line(size=0.5)+
  scale_shape_manual(values=c(15, 16, 17, 6, 8, 11, 0, 1, 2))+
  scale_y_continuous(breaks=seq(0, 1, 0.25), limits=c(0, 1))+
  scale_x_continuous(breaks=seq(17800, 18330, 30), limits=c(17800, 18330))+
  scale_color_manual(values = c("#1F77B4", "#3E6735", "#9467BD", "#2CA02C",  "#F0027F",
                                         "#FDB863", "#202020", "#FF7F0E", "#909090", "#8C564B",
                                         "#B2182B", "#17BECF","#98DF8A", "#AEC7E8", "#DED42F"
  ))+
  theme(panel.background = element_blank())+
  theme(panel.grid.major = element_line(color="grey95"))+
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 0.3))+ #boarder around plot
  theme(legend.key = element_rect(fill = NA))+
  theme(axis.text.y = element_text(size = 12, color = "black"))+
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 90,  hjust=0, vjust=0.5))+
  xlab("Time")+
  ylab("Normalized Average Abundance")+
  ggtitle("Control #2 plot (ITS & 16S)")
#export PDF = 8x4

#
##
#
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
################################################
#### create one MEGA table with ALL RESULTS ####
################################################
#
## DA, TITAN, SCNM, Correlation Network, Taxonomy, FunGuilds...
#
#output table from DA analysis contains *all* ITS and 16S ASVs,
# use this table as a base to build one mega-table with all results
class(ITS16S_DAtable.rn)
head(ITS16S_DAtable.rn)

#
## add TITAN results
TITANhi <- fread(file.choose())
TITANlo <- fread(file.choose())
TITAN1c <- fread(file.choose())
TITAN2c <- fread(file.choose())
colnames(TITANhi)
# Add TITAN indicators to OTUtable with DA values
ITS16S_DA_TITANhi <- merge(ITS16S_DAtable.rn, TITANhi[,c(1,3:18)], by="ASV", all.x=TRUE)
colnames(ITS16S_DA_TITANhi)[2:6] <- paste("DA", colnames(ITS16S_DA_TITANhi)[2:6], sep = "_")
colnames(ITS16S_DA_TITANhi)[7:22] <- paste("TITANhi", colnames(ITS16S_DA_TITANhi)[7:22], sep = "_")
colnames(ITS16S_DA_TITANhi)

ITS16S_DA_TITANlo <- merge(ITS16S_DA_TITANhi, TITANlo[,c(1,3:18)], by="ASV", all.x=TRUE)
colnames(ITS16S_DA_TITANlo)[23:38] <- paste("TITANlo", colnames(ITS16S_DA_TITANlo)[23:38], sep = "_")
colnames(ITS16S_DA_TITANlo)

ITS16S_DA_TITAN1c <- merge(ITS16S_DA_TITANlo, TITAN1c[,c(1,3:18)], by="ASV", all.x=TRUE)
colnames(ITS16S_DA_TITAN1c)[39:54] <- paste("TITAN1c", colnames(ITS16S_DA_TITAN1c)[39:54], sep = "_")
colnames(ITS16S_DA_TITAN1c)

ITS16S_DA_TITAN2c <- merge(ITS16S_DA_TITAN1c, TITAN2c[,c(1,3:18)], by="ASV", all.x=TRUE)
colnames(ITS16S_DA_TITAN2c)[55:70] <- paste("TITAN2c", colnames(ITS16S_DA_TITAN2c)[55:70], sep = "_")
colnames(ITS16S_DA_TITAN2c)

#
## Add SCNM results:
head(BS321Efit_grouped)
ITS16S_DA_TITAN_SCNMhi <- merge(ITS16S_DA_TITAN2c, BS400fit_grouped[,c(1:10)], by="ASV", all.x=TRUE)
colnames(ITS16S_DA_TITAN_SCNMhi)[71:79] <- paste("SCNMhi", colnames(ITS16S_DA_TITAN_SCNMhi)[71:79], sep = "_")
colnames(ITS16S_DA_TITAN_SCNMhi)

ITS16S_DA_TITAN_SCNMlo <- merge(ITS16S_DA_TITAN_SCNMhi, BS321Efit_grouped[,c(1:10)], by="ASV",  all.x=TRUE)
colnames(ITS16S_DA_TITAN_SCNMlo)[80:88] <- paste("SCNMlo", colnames(ITS16S_DA_TITAN_SCNMlo)[80:88], sep = "_")
colnames(ITS16S_DA_TITAN_SCNMlo)

ITS16S_DA_TITAN_SCNM1c <- merge(ITS16S_DA_TITAN_SCNMlo, BS240fit_grouped[,c(1:10)],by="ASV",  all.x=TRUE)
colnames(ITS16S_DA_TITAN_SCNM1c)[89:97] <- paste("SCNM1c", colnames(ITS16S_DA_TITAN_SCNM1c)[89:97], sep = "_")
colnames(ITS16S_DA_TITAN_SCNM1c)

ITS16S_DA_TITAN_SCNM2c <- merge(ITS16S_DA_TITAN_SCNM1c, BS321Wfit_grouped[,c(1:10)], by="ASV",  all.x=TRUE)
colnames(ITS16S_DA_TITAN_SCNM2c)[98:106] <- paste("SCNM2c", colnames(ITS16S_DA_TITAN_SCNM2c)[98:106], sep = "_")
colnames(ITS16S_DA_TITAN_SCNM2c)

#
## Add Correlation Network, nodes, modules, Zi & Pi, etc.
colnames(nodesALL)
ITS16S_DA_TITAN_SCNM_nodes <- merge(ITS16S_DA_TITAN_SCNM2c, nodesALL[,c(1,3:6)], by="ASV", all.x=TRUE)
colnames(ITS16S_DA_TITAN_SCNM_nodes)

#
## Add Taxonomy
colnames(ITS16Stax.df.rn)
ITS16S_DA_TITAN_SCNM_nodes_tax <- merge(ITS16S_DA_TITAN_SCNM_nodes, ITS16Stax.df.rn, by="ASV", all.x=TRUE)

#
## Add FunGuilds
colnames(TaxGuild)
ITS16S_DA_TITAN_SCNM_nodes_taxGuilds <- merge(ITS16S_DA_TITAN_SCNM_nodes_tax, TaxGuild[,c(1,9:17)], by="ASV", all.x=TRUE)
colnames(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds)

#
## Create summary columns...
# add a column to denote in which condition a taxon is DA
ITS16S_DA_TITAN_SCNM_nodes_taxGuilds <- ITS16S_DA_TITAN_SCNM_nodes_taxGuilds %>%
  mutate(DAsummary = case_when(
    DA_adj.P.Val>0.1 ~ "NotDA",
    DA_adj.P.Val<0.1 & DA_logFC > 2 ~ "PosDA", 
    DA_adj.P.Val<0.1 & DA_logFC < -2 ~ "NegDA",
    DA_adj.P.Val<0.1 & DA_logFC > -2 & DA_logFC < 2 ~ "NotDA",
    is.na(DA_adj.P.Val) ~ NA_character_, 
    is.na(DA_logFC) ~ NA_character_
  ))

colnames(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds[, c(22,38,54,70)])

# add a column to summarize TITAN indicators for each ASV:
# Use the "filter" column in the TITAN outputs!
# if filter>0, that means the indicator met purity and reliability criteria
class(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$TITAN2c_filter) #integer
length(which(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$TITANhi_filter > 0))
73+121+161+32 #387 total TITAN indicators (duplicated ASVs)

#change zero's to NAs...
ITS16S_DA_TITAN_SCNM_nodes_taxGuilds[, c(22,38,54,70)][ITS16S_DA_TITAN_SCNM_nodes_taxGuilds[, c(22,38,54,70)] == 0] <- NA

ITS16S_DA_TITAN_SCNM_nodes_taxGuilds <- ITS16S_DA_TITAN_SCNM_nodes_taxGuilds %>%
  mutate(TITANsummary = case_when(
    # ASV is an indicator in one plot only:
    TITANhi_filter > 0 & is.na(TITANlo_filter) & is.na(TITAN1c_filter) & is.na(TITAN2c_filter) ~ "TITANburn",
    is.na(TITANhi_filter) & TITANlo_filter > 0 & is.na(TITAN1c_filter) & is.na(TITAN2c_filter) ~ "TITANburn",
    is.na(TITANhi_filter) & is.na(TITANlo_filter) & TITAN1c_filter > 0 & is.na(TITAN2c_filter) ~ "TITANcontrol",
    is.na(TITANhi_filter) & is.na(TITANlo_filter) & is.na(TITAN1c_filter) & TITAN2c_filter > 0 ~ "TITANcontrol",
    # ASV is an indicator in two plots:
    TITANhi_filter > 0 & TITANlo_filter > 0 & is.na(TITAN1c_filter) & is.na(TITAN2c_filter) ~ "TITANburn",
    TITANhi_filter > 0 & is.na(TITANlo_filter) & TITAN1c_filter > 0 & is.na(TITAN2c_filter) ~ "TITAN_NA",
    TITANhi_filter > 0 & is.na(TITANlo_filter) & is.na(TITAN1c_filter) & TITAN2c_filter > 0 ~ "TITAN_NA",
    is.na(TITANhi_filter) & TITANlo_filter > 0 & TITAN1c_filter > 0 & is.na(TITAN2c_filter) ~ "TITAN_NA",
    is.na(TITANhi_filter) & TITANlo_filter > 0 & is.na(TITAN1c_filter) & TITAN2c_filter > 0 ~ "TITAN_NA",
    is.na(TITANhi_filter) & is.na(TITANlo_filter) & TITAN1c_filter > 0 & TITAN2c_filter > 0 ~ "TITANcontrol",
    # ASV is an indicator in three plots:
    TITANhi_filter > 0 & TITANlo_filter > 0 & TITAN1c_filter > 0 & is.na(TITAN2c_filter) ~ "TITAN_NA",
    TITANhi_filter > 0 & is.na(TITANlo_filter) & TITAN1c_filter > 0 & TITAN2c_filter > 0 ~ "TITAN_NA",
    TITANhi_filter > 0 & TITANlo_filter > 0 & is.na(TITAN1c_filter) & TITAN2c_filter > 0 ~ "TITAN_NA",
    is.na(TITANhi_filter) & TITANlo_filter > 0 & TITAN1c_filter > 0 & TITAN2c_filter > 0 ~ "TITAN_NA",
    # ASV is an indicator in all plots:
    TITANhi_filter > 0 & TITANlo_filter > 0 & TITAN1c_filter > 0 & TITAN2c_filter > 0 ~ "TITAN_NA",
    # ASV is not an indicator in any plot:
    is.na(TITANhi_filter) & is.na(TITANlo_filter) & is.na(TITAN1c_filter) & is.na(TITAN2c_filter) ~ "NA"
  ))

length(which(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$TITANsummary == "TITANburn")) #271
length(which(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$TITANsummary == "TITANcontrol")) #98
length(which(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$TITANsummary == "TITAN_NA")) #5
length(which(is.na(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$TITANsummary))) #0
271+98+5 #374

#
## Add a Network Summary column
ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$NetworkSummary <- ifelse(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$Module==c(1:3), 
                                                              "Network", NA)

#
## Add SCNM Summary columns:
ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNMhiNonNeu <- ifelse(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNMhi_group=="Below", "NonNeutral",
                                                            ifelse(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNMhi_group=="Above", "NonNeutral",
                                                                   NA))
ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNMloNonNeu <- ifelse(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNMlo_group=="Below", "NonNeutral",
                                                            ifelse(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNMlo_group=="Above", "NonNeutral",
                                                                   NA))
ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNM1cNonNeu <- ifelse(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNM1c_group=="Below", "NonNeutral",
                                                            ifelse(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNM1c_group=="Above", "NonNeutral",
                                                                   NA))
ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNM2cNonNeu <- ifelse(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNM2c_group=="Below", "NonNeutral",
                                                            ifelse(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNM2c_group=="Above", "NonNeutral",
                                                                   NA))

ITS16S_DA_TITAN_SCNM_nodes_taxGuilds <- ITS16S_DA_TITAN_SCNM_nodes_taxGuilds %>%
  mutate(SCNMsummary = case_when(
    # ASV is non-neutral in one plot only:
    SCNMhiNonNeu == "NonNeutral" & is.na(SCNMloNonNeu) & is.na(SCNM1cNonNeu) & is.na(SCNM2cNonNeu) ~ "SCNMburn",
    is.na(SCNMhiNonNeu) & SCNMloNonNeu == "NonNeutral" & is.na(SCNM1cNonNeu) & is.na(SCNM2cNonNeu) ~ "SCNMburn",
    is.na(SCNMhiNonNeu) & is.na(SCNMloNonNeu) & SCNM1cNonNeu == "NonNeutral" & is.na(SCNM2cNonNeu) ~ "SCNMcontrol",
    is.na(SCNMhiNonNeu) & is.na(SCNMloNonNeu) & is.na(SCNM1cNonNeu) & SCNM2cNonNeu == "NonNeutral" ~ "SCNMcontrol",
    # ASV is non-neutral in two plots:
    SCNMhiNonNeu == "NonNeutral" & SCNMloNonNeu == "NonNeutral" & is.na(SCNM1cNonNeu) & is.na(SCNM2cNonNeu) ~ "SCNMburn",
    SCNMhiNonNeu == "NonNeutral" & is.na(SCNMloNonNeu) & SCNM1cNonNeu == "NonNeutral" & is.na(SCNM2cNonNeu) ~ "SCNM_NA",
    SCNMhiNonNeu == "NonNeutral" & is.na(SCNMloNonNeu) & is.na(SCNM1cNonNeu) & SCNM2cNonNeu == "NonNeutral"~ "SCNM_NA",
    is.na(SCNMhiNonNeu) & SCNMloNonNeu == "NonNeutral" & SCNM1cNonNeu == "NonNeutral" & is.na(SCNM2cNonNeu) ~ "SCNM_NA",
    is.na(SCNMhiNonNeu) & SCNMloNonNeu == "NonNeutral" & is.na(SCNM1cNonNeu) & SCNM2cNonNeu == "NonNeutral" ~ "SCNM_NA",
    is.na(SCNMhiNonNeu) & is.na(SCNMloNonNeu) & SCNM1cNonNeu == "NonNeutral" & SCNM2cNonNeu == "NonNeutral" ~ "SCNMcontrol",
    # ASV is non-neutral in three plots:
    SCNMhiNonNeu == "NonNeutral" & SCNMloNonNeu == "NonNeutral" & SCNM1cNonNeu == "NonNeutral" & is.na(SCNM2cNonNeu) ~ "SCNM_NA",
    SCNMhiNonNeu == "NonNeutral" & is.na(SCNMloNonNeu) & SCNM1cNonNeu == "NonNeutral" & SCNM2cNonNeu == "NonNeutral" ~ "SCNM_NA",
    SCNMhiNonNeu == "NonNeutral" & SCNMloNonNeu == "NonNeutral"& is.na(SCNM1cNonNeu) & SCNM2cNonNeu == "NonNeutral" ~ "SCNM_NA",
    is.na(SCNMhiNonNeu) & SCNMloNonNeu == "NonNeutral" & SCNM1cNonNeu == "NonNeutral" & SCNM2cNonNeu == "NonNeutral" ~ "SCNM_NA",
    # ASV is non-neutral in all plots:
    SCNMhiNonNeu == "NonNeutral" & SCNMloNonNeu== "NonNeutral" & SCNM1cNonNeu == "NonNeutral" & SCNM2cNonNeu == "NonNeutral" ~ "SCNM_NA",
    # ASV is not non-neutral in any plot:
    is.na(SCNMhiNonNeu) & is.na(SCNMloNonNeu) & is.na(SCNM1cNonNeu) & is.na(SCNM2cNonNeu) ~ "NA"
  ))


colnames(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds)
unique(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$TITANsummary) #TITANburn, TITANcontrol, TITAN_NA, and NA
unique(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$SCNMsummary) #SCNMburn, SCNMcontrol, SCNM_NA, and NA
unique(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds$DAsummary) #PosDA, NegDA, NotDA, and NA (p<0.01, -2>logFC>2)

#
##
###
ITS16S_ResultsMASTER <- data.table(ITS16S_DA_TITAN_SCNM_nodes_taxGuilds)
###
##
#
colnames(ITS16S_ResultsMASTER)
# 1 = ASV
# 2-6 = DA
# 7-70 = complete TITAN outputs
# 71-106 = complete SCNM outputs
# 107-110 = Network (Pi, Zi, Roles, Modules)
# 111-117 = Taxonomy
# 118-126 = FUNGuild
# 127 = DAsummary
# 128 = TITANsummary
# 129 = NetworkSummary
# 130-134 = SCNM summaries

unique(ITS16S_ResultsMASTER$TITANsummary)

ITS16S_ResultsMASTER <- ITS16S_ResultsMASTER %>% mutate_all(na_if,"") #fill blanks with NAs

# use this new table to re-create venn diagram from old Figure 6... 
PyrophileVennSection <- ITS16S_ResultsMASTER[TITANsummary == "TITANburn" & 
                                               SCNMsummary == "SCNMburn" &
                                               DAsummary == "PosDA"] 

nrow(PyrophileVennSection) 
unique(PyrophileVennSection$Genus)
#29 ASVs: TITAN + SCNM + DA_logFC > 2 + DA_adj.P.Val < 0.01   
#15 unique genera:
# Hebeloma, Pyronema, Geopyxis, Geoglossum, Myxomphalia, Botryobasidium, Calyptrozyma, 
# Hyaloscypha, Hebeloma, Xerocomellus, Solicoccozyma, Geminibasidium
# Flavobacterium, Massilia, Gemmatimonas

#
#
ITS16S_ResultsMASTER$Genus[which(ITS16S_ResultsMASTER$ASV == "ASV3949")] <- "Tricharina"
ITS16S_ResultsMASTER$Genus[which(ITS16S_ResultsMASTER$ASV == "ASV1588")] <- "Tricharina"
ITS16S_ResultsMASTER$Genus[which(ITS16S_ResultsMASTER$ASV == "ASV1681")] <- "Lyophyllum"
ITS16S_ResultsMASTER$Genus[which(ITS16S_ResultsMASTER$ASV == "ASV879")] <- "Lyophyllum"
ITS16S_ResultsMASTER$Genus[which(ITS16S_ResultsMASTER$ASV == "ASV921")] <- "Lyophyllum"

#
#
### NUMDATES!
# start = 17816
# Autumn2018 = 17816-17886 **321E burned on 17820
# Winter2019 = 17887-17976 **400 burned on 17898 = start date for comp400
# Spring2019 = 17977-18068
# Summer2019 = 18069-18160
# Autumn2019 = 18161-18251
# Winter2020 = 18252-18342 **321W burned on 18308
# end = 18308

