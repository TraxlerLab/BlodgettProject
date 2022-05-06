library(decontam)
library(dada2)
library(phyloseq)
library(ggplot2)
library(data.table)
library(tidyverse)
library(vegan)
library(reshape2)
####
##
#
### SOME USEFUL/IMPORTANT DATA OBJECTS
#
# phyloseq object before decontam:
ps
# phyloseq object output from decontam:
ps2
# complete taxonomy table extracted from decontam output:
head(tax.df)
# complete sequence table extracted from decontam output:
head(seq.df, 2)

# complete OTU table extracted from decontam output (rownames = seqID):
otu.df[1:5,1:5]
# complete OTU table extracted from decontam output (column1 = seqID):
otu.df.rownames[1:5,1:5]
# transposed OTU table:
otu.dft[1:5,1:5]
# OTU table with taxonomic classes appended to the end of the table
OTUtax.df[1:5, 1:5]
OTUtax.df[1:5, 662:672]
# OTUtax table, but with obscure phyla removed
OTUtax.df.sub[1:5, 662:672]

# Metadata
head(MetadataWithReads)

#######################################################################################
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

#################################################
#### CREATE PHYLOSEQ OBJECT and RUN DECONTAM ####
#################################################
##
##
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html#identifying-contaminants-in-marker-gene-and-metagenomics-data
##
##

library(phyloseq)
library(ggplot2)
library(decontam)
library(data.table)
library(tidyverse)
library(vegan)

## MASTER METADATA TABLE
SampleMetadata <- fread("~/Desktop/Box Sync/All16S_Processed/SampleMetadata_Blodgett_UPDATED3.csv")
head(SampleMetadata)
## WARNING -- the sequencing facility changed some of our sample names...
## here's how to fix it so that our Metadata Table matches the seq facility output:
SampleMetadata$SeqID <- gsub("rp", "REP", SampleMetadata$SeqID) #find and replace "rep" with "rp" to match sequencing facility output...
SampleMetadata$SeqID <- gsub("p", "P", SampleMetadata$SeqID) #find and replace "p" with "nothing"P" to match sequencing facility output...
## CORRECT DATE SNAFU! ("replace any '8-Oct-20' with '8-Oct-19'")
SampleMetadata$Date <- gsub("8-Oct-20", "8-Oct-19", SampleMetadata$Date) 
## CORRECT POOLED LABEL SNAFU! ("replace any 'p1-5' with 'p1-3'")
SampleMetadata$depth <- gsub("p1-5", "p1-3", SampleMetadata$depth) 

write.csv(SampleMetadata, "~/Desktop/Blodgett Data & Manuscript/Blodgett_16S_Comm/SampleMetadata_Blodgett_UPDATED.csv")
#NOTE that there were a handful more edits made in Excel
#USE THIS TABLE:
SampleMetadata <- fread("~/Desktop/Box Sync/All16S_Processed/SampleMetadata_Blodgett_UPDATED3.csv")
head(SampleMetadata)
##16S Table, Fixed
SampleMetadata <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/16S_SampleMetadata.csv")
SampleMetadata <- SampleMetadata[-c(651:nrow(SampleMetadata)),]
SampleMetadata <- column_to_rownames(SampleMetadata, var = "SeqID")
head(SampleMetadata)
class(SampleMetadata)

## decontam uses dada2 outputs converted into phyloseq objects, and phyloseq requires OTU and Taxa tables as matrices
#
# Load OTU table:
seqtab.nochim <- as.matrix(read.delim("~/Desktop/Box Sync/All16S_Processed/QIIMEobjects/exported-feature-table/feature-table.tsv", sep = "\t", header = TRUE, skip = 1, row.names = 1))

seqtab.nochim[1:5,1:2] #since there's a million columns, this is better than head()

seqtab.nochim.t <- t(seqtab.nochim)
seqtab.nochim.t[1:5,1:2]

# Load Taxonomy Table:
taxa <- as.matrix(fread("~/Desktop/Box Sync/All16S_Processed/QIIMEobjects/TaxonomySep.csv",
                        header = TRUE), rownames=1)
head(taxa)


# Extract the SeqIDs from the OTUtable
samplenames <- data.table(SeqID = row.names(seqtab.nochim.t))

# add read counts per sample to the Sample Metadata Table.. (reads output by DADA2 in the post-chimera-removal sanity check)
reads <- fread("~/Desktop/Box Sync/All16S_Processed/QIIMEobjects/stats-metadata_3C.csv") 
head(reads)
colnames(reads) <- c("SeqID", "dada2_reads", "nochim_reads")
reads <- mutate(reads, SeqID = str_replace_all(SeqID, "-", "."))


class(reads$nochim_reads)
class(reads$dada2_reads)
reads$nochim_reads <- as.numeric(reads$nochim_reads)
max(reads$nochim_reads)
min(reads$nochim_reads)

# subset metadata table, add reads counts, and create table for phyloseq
dat1 <- merge(samplenames, SampleMetadata, by="SeqID", all.x=TRUE)
dat2 <- merge(dat1, reads, by="SeqID", all.x=TRUE)
samp16S <- column_to_rownames(dat2, var="SeqID") #rownames = SeqID will be required for phyloseq

test.ps <- phyloseq(sample_data(samp16S))


#plot the number of reads by sample to get an general idea of what the data look like
ggplot(data=dat2, aes(x=SeqID, y=nochim_reads, color=Who)) +
  geom_point() + 
  scale_color_manual(values = c("#00ff95", "#ff9500", "#9500ff", "#ff006a")) +
  scale_y_continuous(expand= c(0,0), breaks=seq(0, 130000, 10000), limits=c(0, 130000)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))

ggplot(data=dat2, aes(x=SeqID, y=nochim_reads, color=SamplePool)) +
  geom_point() + 
  scale_color_manual(values = c("#00ff95", "#ff9500", "#9500ff")) +
  scale_y_continuous(expand= c(0,0), breaks=seq(0, 130000, 10000), limits=c(0, 130000)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1))

##
# CREATE PHYLOSEQ OBJECT
ps <- phyloseq(otu_table(seqtab.nochim.t, taxa_are_rows=FALSE), #matrix (rownames = SeqIDs, colnames = ASVs, values = abundance of each ASV in each SeqID)
               tax_table(taxa), #matrix (rownames = ASVs, colnames = taxonomic levels, values = taxonomic classification for each ASV)
               sample_data(samp16S)) #data.frame, (rownames = SeqIDs, colnames & values = additional info for each SeqID)

# check that nothing catastrophic happened when creating the phyloseq object, table dimensions should match:
dim(tax_table(ps)) #phyloseq object
dim(taxa) #original
(otu_table(ps))[1:5, 1:5] #phyloseq object  
dim(seqtab.nochim) #original         
dim(sample_data(ps)) #phyloseq object
head(samp16S) #original
head(tax_table(ps))

otu <- otu_table(ps)
class(otu)
otu.t <- t(otu)
otu.t[1:5, 1:5]

## ASVs are currently identified by their DNA sequence, which is cumbersome
## Move this DNA sequence to the refseq slot of the phyloseq object
## and give each ASV a simpler name like ASV1, ASV2, etc.
dna <- Biostrings::DNAStringSet(taxa_names(ps2copy))
names(dna) <- taxa_names(ps2copy)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("OTU", seq(ntaxa(ps)))
ps

# To run decontam, we will need a column ("is.neg") in the Metadata Table;
# control samples are indicated as "TRUE" and all other samples are "FALSE"
# Create the is.neg column:
sample_data(ps)$is.neg <- sample_data(ps)$Who == "Control"
head(sample_data(ps)) #check that it looks correct
view(sample_data(ps))

#
## RUN DECONTAM!
# outputs a table of stats for all taxa and decontam's decision about whether or not each taxon is a contaminant
# play around with changing the "threshold" value..
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold = 0.1)
# Summary table: how many true contaminants were identified?
table(contamdf.prev$contaminant) 
#tresh0.1 = 190contaminants, tresh0.05 = 93contaminants
#for 16S, tresh0.5 = 103contaminants, tresh0.1 = 361contaminants, tresh0.05 = 34contaminants 

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
## Each ASV is assigned a P-value/score
## In general, the threshold P value should be less-than-or-equal-to 0.5
## higher scores = likely not a contaminant

unique(sample_data(ps)$Who)
sum(sample_data(ps)$Who == "Control") #46 #for 16S = 42
sum(sample_data(ps)$Who == "Neemonika") #325 #for 16S = 316
sum(sample_data(ps)$Who == "Neemonika", sample_data(ps)$Who =="Phillip", sample_data(ps)$Who =="PhillipPool") #618 ...618+46=664 #for 16S = 608

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
  #geom_abline(intercept = 0, slope = 1)+ #add a line showing the threshold for equal representation in controls and samples
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")+
  ggtitle("Threshold = 0.1")


# REMOVE CONTAMINANTS FROM DATASET
contam <- rownames_to_column(contamdf.prev, var = "rowname") #move rownames (ASVs) to a column that data.table can work with
contam.dt <- data.table(contam) #convert to a data.table
badTaxa <- contam.dt[contaminant == "TRUE", rowname] #use data.table to generate a list of contaminant ASVs
goodTaxa <- setdiff(taxa_names(ps), badTaxa) #list of all ASVs minus the contaminant ASVs
ps2 <- prune_taxa(goodTaxa, ps) #create a new phyloseq object that only contains the good (non-contaminant) ASVs!
#check out the difference:
ps #15216 taxa #16S = 95178
ps2 #15026 taxa #16S = 94817
95178 - 94817 # =361, which is equivalent to the number of taxa that decontam identified as bad! this worked! hooray!

########################################
##### ORDINATION PLOTS & PERMANOVA  #####
#########################################
library(vegan)
library(phyloseq)
library(ggplot2)
library(data.table)
## NMDS represents as well as possible the ordering relationships among objects in a small number of axes
# NMDS does not preserve the exact dissimilarities among objects
# The user sets the number of axes,
# the the NMDS algorithm fits the data into the specified number or axes
# Ouput stress value describes the fit of the data to the axes
# The goal is to minimize stress (different number of axes will change stress)
# a random R ecology blog says stress >0.3 is bad, stress <0.01 great
#...it seems the default for phyloseq is 2 axes?
##
# PCA assumes data are normally distributed, but is more infomative than NMDS because the axes are meaningful
# The axes are meaningful in PCA because the exact dissimilarities among objects is maintained,
# and there can be as many (or as few) different axes as are needed to account for all the variation in the data
# RDA is specific type of PCA, in which a regression is preformed prior to the PCA:
# "In RDA, one can truly say that the axes explain or model (in the statistical sense) the variation of the..
# ..dependent matrix. Furthermore, a global hypothesis (H0) of absence of linear relationship between Y and X can..
# ..be tested in RDA; this is not the case in PCA"
# Running the rda() function without an environmental variable matrix is just a plain old PCA
# Hellinger Transformation can be used to normalize data prior to doing a PCA
# PERMANOVA can then be used to identify whether or not there is a significant difference between treatments
# View the entire summary of a PCA with the summary() command... like this:
pca <- rda(data)
summary(pca)

#### HELLLINGER-TRANSFORMATION + PCA

# subset data
ps.F240pp <- prune_samples(sample_data(ps2)$Date== "2020-02-17" &
                             sample_data(ps2)$Plot== "240" &
                             sample_data(ps2)$Who== "Phillip", ps2)


#ps.F240pp <- prune_samples(ps2, Who== "Phillip", Plot== "240", Date== "2/17/20")
ps.F240pp

ps2.F240pp <- prune_samples(sample_data(ps2)$Who== "Phillip", ps2)
unique (sample_data(ps.F240pp)$depth)
unique (sample_data(ps.F240pp)$Who)

##ps2.F240 <- prune_samples(sample_data(ps2)$Date== "17-Feb-20" & sample_data(ps2)$Plot== "240", ps2)
##ps2.F240 #54 samples, 7622 taxa

# copy phyloseq object for Hellinger Transformation
ps.F240pp.hel <- ps.F240pp
# Hellinger Transformation
otu_table(ps.F240pp.hel) <-otu_table(decostand(otu_table(ps.F240pp.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps.F240pp.hel.pca <- rda(otu_table(ps.F240pp.hel))
# Broken Stick Model
screeplot(ps.F240pp.hel.pca, bstick = TRUE,
          npcs = length(ps.F240pp.hel.pca$CA$eig),
          main = "Hellinger PCA - February (no burn control) comp240") #1500x600
# PCA plot by Rep
plot_ordination(ps.F240pp.hel, ps.F240pp.hel.pca, color="Rep")+
  stat_ellipse(geom="polygon", aes(fill=Rep), alpha=0.2) +
  geom_point(size=2) +
  geom_text(label=sample_data(ps.F240pp.hel)$depth, color="black", nudge_x = 0.02, nudge_y = 0.02)+
  ggtitle("Hellinger PCA - February (no burn control) comp240") #600x500
#PERMANOVA to test for a statistically significant differences between reps:
adonis(otu_table(ps.F240pp.hel) ~ sample_data(ps.F240pp)$Rep)
# PCA plot by depth
plot_ordination(ps.F240pp.hel, ps.F240pp.hel.pca, color="depth")+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=depth))+ 
  geom_point(size=2) + 
  theme_classic()+
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black", "pink"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black", "pink"))+
  ggtitle("Hellinger PCA - February (no burn control) comp240") #600x500

#PERMANOVA to test for a statistically significant differences between reps:
adonis(otu_table(ps.F240pp.hel) ~ sample_data(ps.F240pp)$depth)


# subset data
ps.F321pp <- prune_samples(sample_data(ps2)$Date== "2020-02-17" & 
                             sample_data(ps2)$Plot== "321west", ps2)
ps.F321pp #53 samples, 7622 taxa
unique (sample_data(ps.F321pp)$depth)
# copy phyloseq object for Hellinger Transformation
ps.F321pp.hel <- ps.F321pp
# Hellinger Transformation
otu_table(ps.F321pp.hel) <-otu_table(decostand(otu_table(ps.F321pp.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps.F321pp.hel.pca <- rda(otu_table(ps.F321pp.hel))
# Broken Stick Modle:
screeplot(ps.F321pp.hel.pca, bstick = TRUE,
          npcs = length(ps.F321pp.hel.pca$CA$eig),
          main = "Hellinger PCA - February (post-burn) comp321") #1500x600
# PCA plot by Rep:
plot_ordination(ps.F321pp.hel, ps.F321pp.hel.pca, color="Rep")+
  stat_ellipse(geom="polygon", aes(fill=Rep), alpha=0.2) +
  geom_point(size=2) +
  geom_text(label=sample_data(ps.F321pp.hel)$depth, color="black", nudge_x = 0.02, nudge_y = 0.02)+
  ggtitle("Hellinger PCA - February (post-burn) comp321") #600x500
# PERMANOVA to test for a statistically significant differences between reps:
adonis(otu_table(ps.F321pp.hel) ~ sample_data(ps.F321pp)$Rep)
# PCA plot by depth:
plot_ordination(ps.F321pp.hel, ps.F321pp.hel.pca, color="depth")+
  geom_mark_hull(concavity = 5,expand=0,radius=0,aes(fill=depth))+ 
  geom_point(size=2) + 
  theme_classic()+
  scale_color_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                     values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black", "pink", "purple"))+
  scale_fill_manual(breaks=c("1","2","3","4", "5", "10", "20", "p1-3", "p1-10"),
                    values=c("#ff3b9d", "#ff3b6c", "#ff3b3b", "#ff6c3b", "#ff9d3b", "#ffce3b", "#ceff3b", "#0077ee", "black", "pink", "purple"))+
  ggtitle("Hellinger PCA - February (post-burn) comp321") #600x500

### For Blodgett Data ###

ps.BS <- prune_samples(sample_data(ps2)$Plot== "240" &
                         sample_data(ps2)$Who== "Neemonika" |
                         sample_data(ps2)$Plot== "400" &
                         sample_data(ps2)$Who== "Neemonika" |
                         sample_data(ps2)$Plot== "321west" &
                         sample_data(ps2)$Who== "Neemonika" |
                         sample_data(ps2)$Plot== "321east" &
                         sample_data(ps2)$Who== "Neemonika", ps2)
ps.BS
unique (sample_data(ps.BS)$Who)


# copy phyloseq object for Hellinger Transformation
ps.BS.hel <- ps.BS
# Hellinger Transformation
otu_table(ps.BS.hel) <-otu_table(decostand(otu_table(ps.BS.hel), method = "hellinger"), taxa_are_rows=FALSE)
# PCA
ps.BS.hel.pca <- rda(otu_table(ps.BS.hel))
# Broken Stick Model
screeplot(ps.BS.hel.pca, bstick = TRUE,
          npcs = length(ps.BS.hel.pca$CA$eig),
          main = "Hellinger PCA - All Blodgett Samples") #1500x600

### Define the order in which the Plot variable should be plotted:
sample_data(ps.BS.hel)$Plot <- factor(sample_data(ps.BS.hel)$Plot, levels=c("240", "321west", "321east", "400"))
### PCA plot with ellipses:
library(ggplot2)
plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="Plot")+
  stat_ellipse(geom="polygon", aes(fill=Plot), alpha=0.05)+
  #scale_color_gradient2(midpoint=6, low="green", mid="orange", high="blue")+ ##for plotting continuous variables!
  scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Hellinger PCA - 16S") 
### Export as .png or .tiff: 600x500
### Export as .pdf TBD
#PERMANOVA to test for a statistically significant differences between reps:
adonis(otu_table(ps.BS.hel) ~ sample_data(ps.BS)$Plot)

# PCA plot by depth
plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="depth")+
  stat_ellipse(geom="polygon", aes(fill=Plot), alpha=0.05)+
  geom_point(size=2) +
  scale_color_manual(breaks=c("0-3cm","3-6cm","10cm"),
                     values=c("#ff0102", "#0077ee", "#00cd9a","#ffce3b"))+
  scale_fill_manual(breaks=c("0-3cm","3-6cm","10cm"),
                    values=c("#ff0102", "#0077ee", "#00cd9a", "#ffce3b"))+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Depth") #600x500
#PERMANOVA to test for a statistically significant differences between reps:
adonis(otu_table(ps.BS.hel) ~ sample_data(ps.BS)$Rep)

# PCA plot by Fire
plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="Fire")+
  stat_ellipse(geom="polygon", aes(fill=Fire), alpha=0.05)+
  #scale_color_gradient2(midpoint=6, low="green", mid="orange", high="blue")+ ##for plotting continuous variables!
  scale_fill_manual(values=c("#00cd9a", "#ff0102"),
                    labels=c("Pre-Fire", "Post-Fire"))+
  scale_color_manual(values=c("#00cd9a", "#ff0102"),
                     labels=c("Pre-Fire", "Post-Fire"))+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Burned vs. Unburned") 
#PERMANOVA
adonis(otu_table(ps.BS.hel) ~ sample_data(ps.BS)$Fire)

## Oval plots
plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="Season")+
  stat_ellipse(geom="polygon", aes(fill=Plot), alpha=0.05)+
  #scale_color_gradient2(midpoint=6, low="green", mid="orange", high="blue")+ ##for plotting continuous variables!
  scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                    labels=c("1c", "2c", "lo", "hi"))+
  scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
                     labels=c("1c", "2c", "lo", "hi"))+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Hellinger PCA - All Neemonika (3-6cm samples omitted)") 

plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="Season")+
  stat_ellipse(geom="polygon", aes(fill=Season), alpha=0.2)+
  geom_point(size=2) +
  scale_color_manual(breaks=c("Autumn","Winter","Spring","Summer"),
                     values=c("#ff6c3b", "#0077ee", "#00cd9a", "#ff3b9d"))+
  scale_fill_manual(breaks=c("Autumn","Winter","Spring","Summer"),
                    values=c("#ff6c3b", "#0077ee", "#00cd9a", "#ff3b9d"))+
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Season") 


##Depth
plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="depth")+
  stat_ellipse(geom="polygon", aes(fill=depth), alpha=0.2)+
  geom_point(size=2) +
  theme_classic()+
  scale_color_manual(breaks=c("0-3cm","3-6cm","10cm"),
                     values=c("#ff0102", "#0077ee", "#00cd9a","#ffce3b"))+
  scale_fill_manual(breaks=c("0-3cm","3-6cm","10cm"),
                    values=c("#ff0102", "#0077ee", "#00cd9a", "#ffce3b"))+
  ggtitle("Depth") 

unique (sample_data(ps2.BS)$depth)

### PCA plot with gradient colors
plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="CtempDailyAvg")+
  #stat_ellipse(geom="polygon", aes(fill=pH), alpha=0.2)+
  scale_color_gradient2(midpoint=10, low="green", mid="orange", high="blue")+
  geom_point(size=2) +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Avg. Daily Air Temp (C)") #600x500

plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="Precip4wkSum_in")+
  #stat_ellipse(geom="polygon", aes(fill=pH), alpha=0.2)+
  scale_color_gradient2(midpoint=10, low="green", mid="orange", high="blue")+
  geom_point(size=2) +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Avg. Precip - 4wkSum") #600x500

plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="pH")+
  #stat_ellipse(geom="polygon", aes(fill=pH), alpha=0.2)+
  #stat_ellipse(geom="polygon", aes(fill=Fire), alpha=0.2)+
  scale_color_gradient2(midpoint=6, low="green", mid="orange", high="blue")+
  geom_point(size=2) +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("pH") #600x500
#permanova
adonis(otu_table(ps.BS.hel) ~ sample_data(ps2.Neemonika)$pH)

plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="Precip4wkSum_in")+
  #stat_ellipse(geom="polygon", aes(fill=Fire), alpha=0.2)+
  scale_color_gradient2(midpoint=15, low="green", mid="orange", high="blue")+
  geom_point(size=2) +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Avg. Weekly Preciptation") #600x500
#permanova
adonis(otu_table(ps.BS.hel) ~ sample_data(ps2.Neemonika)$Precip4wkSum_in)

plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="numdate")+
  #stat_ellipse(geom="polygon", aes(fill=pH), alpha=0.2)+
  scale_color_gradient2(midpoint=18000, low="green", mid="orange", high="blue")+
  geom_point(size=2) +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Blodgett Samples by Time") #600x500
#permanova
adonis(otu_table(ps.BS.hel) ~ sample_data(ps2.Neemonika)$numdate)

plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="DaysSinceBurn")+
  scale_color_gradient2(midpoint=100, low="green", mid="orange", high="blue")+
  geom_point(size=2) +
  theme_classic()+
  theme(axis.text = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  theme(legend.text = element_text(size=14, color="black"))+
  theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Time after Treatment") #600x500
#permanova
adonis(otu_table(ps.BS.hel) ~ sample_data(ps2.Neemonika)$DaysSinceBurn)

##
plot_ordination(ps.BS.hel, ps.BS.hel.pca, color="pH")+
  #stat_ellipse(geom="polygon", aes(fill=Plot), alpha=0.05)+
  scale_color_gradient2(midpoint=6, low="green", mid="orange", high="blue")+ ##for plotting continuous variables!
  geom_point(size=4) +
  #scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
  #                  labels=c("1c", "2c", "lo", "hi"))+
  #scale_color_manual(values=c("#009acd", "#00cd9a", "#ff7f01", "#ff0102"),
  #                   labels=c("1c", "2c", "lo", "hi"))+
  theme_classic()+
  #theme(axis.text = element_text(color = "black", size=14))+
  #theme(axis.title = element_text(color = "black", size=14))+
  #theme(legend.text = element_text(size=14, color="black"))+
  #theme(legend.title = element_text(size=14, face="bold", color="black"))+
  ggtitle("Hellinger PCA - pH") 
### Export as .png or .tiff: 600x500
### Export as .pdf TBD

### NUMDATES!
#
# Autumn2018 = 17816-17886 **321E burned on 17820
# Winter2019 = 17887-17976 **400 burned on 17898
# Spring2019 = 17977-18068
# Summer2019 = 18069-18160
# Autumn2019 = 18161-18251
# Winter2020 = 18252-18342 **321W burned on 18308
#

#
###################################
## Diversity Metrics across time ##
###################################
#
# subset ps2 (decontam output) for only Neemonika data,
# up until the day that 321W burned, and exclude comp380:
ps2.Neemonika <- prune_samples(sample_data(ps2)$Who == "Neemonika" &     #only Neemonika data (no Phillip or Controls)
                                 sample_data(ps2)$Plot != "380" &  
                                 sample_data(ps2)$numdate < 18308, ps2)  #exclude dates on/after 321W burn

#ps2.Neemonika <- prune_samples(sample_data(ps2)$Who == "Neemonika" &     #only Neemonika data (no Phillip or Controls)
#                                 sample_data(ps2)$Plot != "380" &  
#                                 sample_data(ps2)$depth != "3-6cm" & #exclude comp380 data
#                                 sample_data(ps2)$numdate < 18308, ps2)  #exclude dates on/after 321W burn


# Extract OTU table from phyloseq into a data.frame:
otu.df <- as.data.frame(otu_table(ps2.Neemonika))
otu.df[1:5,1:5] # rownames=SeqIDs, colnames=ASVs
dim(otu.df) ##255 samples.

SampleMetadata <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/16S_SampleData_decontam_620.csv")
# Create Metadata tables for each Plot, which will be useful for subsetting the otu.df table!
Metadata321W <- SampleMetadata[Who =="Neemonika" & Plot == "321west" & numdate < 18308]
Metadata321E <- SampleMetadata[Who =="Neemonika" & Plot == "321east" & numdate < 18308]
Metadata400 <- SampleMetadata[Who =="Neemonika" & Plot == "400" & numdate < 18308]
Metadata240 <- SampleMetadata[Who =="Neemonika" & Plot == "240" & numdate < 18308]
MetadataBS <- SampleMetadata[Who =="Neemonika" & Plot != "380" & numdate < 18308]
head(MetadataBS)
dim(MetadataBS)
#
# Subset otu_table for only the rownames that match the SeqIDs in MetadataBS:
otu_tableBSasv <- otu.df[rownames(otu.df) %in% MetadataBS$SeqID, ]
otu_tableBSasv[1:10, 1:5] #all SeqIDs should be from comp400!
dim(otu_tableBSasv) #255 samples.

# Calculate diversity metric with the vegan package
ShanDivBS <- as.data.table(diversity(otu_tableBSasv, index="shannon")) #index options = shannon, simpson, or invsimpson
ShanDivBS[ ,Date:=MetadataBS$numdate]
ShanDivBS[ ,Plot:=MetadataBS$Plot]
ggplot(ShanDivBS, aes(x=numdate, y=V1, color=Plot))+
  geom_point()+
  geom_smooth(aes(fill=Plot))+
  theme_minimal()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black"))+
  scale_x_continuous(breaks=seq(17890, 18310, 30), limits=c(17890, 18310))+
  xlab("Time")+
  ylab("Shannon Diversity Index")+
  ggtitle("Blodgett 16S Community Shannon Diversity")

SimpDivBS <- as.data.table(diversity(otu_tableBSasv, index="simpson"))
SimpDivBS[ ,Date:=MetadataBS$numdate]
SimpDivBS[ ,Plot:=MetadataBS$Plot]
ggplot(SimpDivBS, aes(x=Date, y=V1, color=Plot))+
  geom_point()+
  geom_smooth(aes(fill=Plot))+
  theme_minimal()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black"))+
  scale_x_continuous(breaks=seq(17890, 18310, 30), limits=c(17890, 18310))+
  xlab("Time")+
  ylab("Simpson Diversity Index")+
  ggtitle("Blodgett 16S Community Simpson Diversity")


#calculate species RICHNESS
# Where ASV abundance value is greater than 0, give it a 1, then sum the 1's by sample:
otu_tableBSasv.rich <- apply(otu_tableBSasv > 0,1,sum)
otu_tableBSasv.rich.dt <- as.data.table(otu_tableBSasv.rich)
otu_tableBSasv.rich.dt[ ,Plot:=MetadataBS$Plot]
otu_tableBSasv.rich.dt[ ,Date:=MetadataBS$numdate]
names(otu_tableBSasv.rich.dt)[1] <- "Richness"
ggplot(otu_tableBSasv.rich.dt, aes(x=Date, y=Richness, color=Plot))+
  geom_point()+
  geom_smooth(aes(fill=Plot))+
  theme_minimal()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black"))+
  scale_x_continuous(breaks=seq(17890, 18310, 30), limits=c(17890, 18310))+
  xlab("Time")+
  ylab("Richness")+
  ggtitle("Blodgett 16S Community Richness")


#calculate species EVENNESS
otu_tableBSasv.even <- as.data.table(diversity(otu_tableBSasv, index="simpson"))/log(otu_tableBSasv.rich)
otu_tableBSasv.even[ ,Date:=MetadataBS$numdate]
otu_tableBSasv.even[ ,Plot:=MetadataBS$Plot]
names(otu_tableBSasv.even)[1] <- "Evenness"
ggplot(otu_tableBSasv.even, aes(x=Date, y=Evenness, color=Plot))+
  geom_point()+
  geom_smooth(aes(fill=Plot))+
  theme_minimal()+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black"))+
  scale_x_continuous(breaks=seq(17890, 18310, 30), limits=c(17890, 18310))+
  xlab("Time")+
  ylab("Evenness")+
  ggtitle("Blodgett 16S Community Evenness")


##############################
####### TITAN Analysis #######
##############################
library(TITAN2)

# TITAN2 requires two inputs:
# (1) OTU table -- data.frame: rownames = samples, colnames=OTUs
View(glades.taxa) #example
# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
View(glades.env) #example


########################################
####### COMBINED TITAN Analysis #######
#######################################
library(TITAN2)

# TITAN2 requires two inputs:
# (1) OTU table -- data.frame: rownames = samples, colnames=OTUs
View(glades.taxa) #example
# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
View(glades.env) #example

SampleMetadata <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/16S_SampleData_decontam_620.csv")
SampleMetadataITS <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/ITS_Data/Blodgett_CommData_Decontam-selected/ITS_ps2_SampleData_decontamTITAN.csv")

SampleMetadata$SeqID <-gsub("\\..*","",SampleMetadata$SeqID)
SampleMetadataITS$SeqID <-gsub("\\-.*","",SampleMetadataITS$SeqID)

view(SampleMetadata)
view(SampleMetadataITS)

dim(SampleMetadata)
dim(SampleMetadataITS)

AllSampleMetadata <- cbind(SampleMetadata, SampleMetadataITS, by="SeqID")

#Change Variant IDs to OTUs
#taxa_names(ps2) <- paste0("OTU", seq(ntaxa(ps2)))

# Subset phyloseq object: 
#ps2.Neemonika <- prune_samples(sample_data(ps2)$Who == "Neemonika" & 
#                                 sample_data(ps2)$Plot != "380" &
#                                 sample_data(ps2)$numdate < 18308, ps2)
# Extract OTU table:
#otu.df <- as.data.frame(otu_table(ps2.Neemonika))
#otu.df[1:5,1:5]

#otu2.df <- as.data.frame(otu_table(ps3.Neemonika))
#otu2.df[1:5,1:5]

#Stringtable <- colSums(otu.df)
#OTUtable <- colSums(otu2.df)

#Stringtable.df <- as.data.frame(Stringtable)
#OTUtable.df <- as.data.frame(OTUtable)

#Stringtable.df <- rownames_to_column(Stringtable.df, var= "String")
#OTUtable.df <- rownames_to_column(OTUtable.df, var= "OTU")

#COMBO <- cbind(Stringtable.df, OTUtable.df)

#if(COMBO$Stringtable == COMBO$OTUtable)
#{
#  print("Column A and B are identical")
#}

#identical(COMBO[['Stringtable']],COMBO[['OTUtable']])

#write.csv(COMBO, file = "~/Desktop/Box Sync/Blodgett_16S_Comm/16S_OTU_String_Combo.csv")

#Read ITS/16S OTU table:
otu <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/otu_table_16SITS.csv")
otu.df <- as.data.frame(otu)
otu.df$V1 <- NULL
colnames(otu.df)[1] <- "OTU"
otu.df[1:5,1:5] 
dim(otu.df) #255 rows, 110139 columns
otu.df <- column_to_rownames(otu.df, var = "OTU")

#Read ITS/16S Taxa table
tax <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/tax_table_16SITS.csv", header = TRUE)
tax.df <- as.data.frame(tax)
tax.df[2:5, 1:8]
tax.df$V1 <- NULL
dim(tax.df) #110138 rows, 8 columns
colnames(tax.df)[1] <- "OTU"
rownames(tax.df) <- tax.df[,1]
tax.df$OTU <- NULL
#tax.df <- tax.df[,-1]

head(SampDat.dt)
rownames(SampDat.dt) <- SampDat.dt[,2]

ps <- phyloseq(otu_table(otu.dft, taxa_are_rows=TRUE), #matrix (rownames = SeqIDs, colnames = ASVs, values = abundance of each ASV in each SeqID)
               tax_table(tax.df), #matrix (rownames = ASVs, colnames = taxonomic levels, values = taxonomic classification for each ASV)
               sample_data(AllSampleMetadata)) #data.frame, (rownames = SeqIDs, colnames & values = additional info for each SeqID)


# Extract Taxonomy table:
#tax.df <- as.data.frame(tax_table(ps2.Neemonika))
#tax.df[1:5,1:5]
#transpose OTU table:
otu.dft <- as.data.frame(t(otu.df))
otu.dft[1:5,1:5]
# Merge taxa info onto OTU table by matching rownames
OTUtax <- as.data.table(merge(otu.dft, tax.df, by=0))
dim(OTUtax) #110138 rows, 264 columns
OTUtax[1:5,257:264] #taxonomy should be at the end of the table

# Create the .taxa table for TITAN; with Genus appened onto ASVs,
# which will make the TITAN output plots a little easier to look at...
# So the y-axis will be labeled like "Pyronema_ASV257" instead of just "ASV257"
# Note that you can choose whatever taxonomic level you want for this! 
# I feel Genus is most informative for the Fungi, 
# but maybe Phylum or a different level is more informative for bacteria... ?!
#
# Split the "Genus" column to get rid of the "g__" before every genus name
OTUtax[ ,c("JustGenus")] <- copy(OTUtax$Genus)
OTUtax$JustGenus <- gsub("_", "-", OTUtax$JustGenus) # replace all "_" with "-" in the JustGenus column

dim(OTUtax) #110138, 265
head(OTUtax)
OTUtax[1:5,257:265]
OTUtax[1:5, 1:5]
# Create column with GenusASV info (i.e. "Pyronema_ASV257")
OTUtax <- unite(OTUtax, GenusOTU, c(JustGenus, Row.names), sep="_", remove=FALSE)
# Massage table into the format that TITAN is expecting:
#remove excess columns (by number) so that all that remains are the sample columns and the GenusASV column
dim(OTUtax) #110138,266
OTUtax[1:5,258:266] #note which columns are non-numeric, taxonomy-related columns...
OTUtax[1:5,1:5] #note that GenusASV is column#1
OTUtax[,c(2,258:266)] <- list(NULL) #remove all taxonomy columns
dim(OTUtax) #110138, 256
OTUtax[1:5,250:256] #check to make sure taxonomy columns were successfully removed
OTUtax[1:5, 1:5] #check to make sure GenusASV column is still there
OTUtax <- column_to_rownames(OTUtax, var="GenusOTU") #move GenusASV to the rownames
#transpose the table
OTUtax.t <- as.data.frame(t(OTUtax))
OTUtax.t[1:5, 1:5] #check that rownames=SeqIDs and colnames=GenusASV
OTUtax.t.rn <- rownames_to_column(OTUtax.t, var="SeqID") #move rownames to a column
OTUtax.t.rn[1:5, 1:5] ## THIS is the table that should be used as the basis for .taxa argumens in TITAN!
#
#
#
# Run TITAN on each Plot separately!
#
Metadata321W <- AllSampleMetadata[Who =="Neemonika" & Plot == "321west" & numdate < 18308]
Metadata321E <- AllSampleMetadata[Who =="Neemonika" & Plot == "321east" & numdate < 18308]
Metadata400 <- AllSampleMetadata[Who =="Neemonika" & Plot == "400" & numdate < 18308]
Metadata240 <- AllSampleMetadata[Who =="Neemonika" & Plot == "240" & numdate < 18308]
MetadataBS <- AllSampleMetadata[Who =="Neemonika" & Plot != "380" & numdate < 18308]
########### comp400 TITAN! 
#
# (1) OTU Table -- Subset OTUtax.t.rn table for only the column names that match the SeqIDs in Metadata400
BS400.taxa <- OTUtax.t.rn[OTUtax.t.rn$SeqID %in% Metadata400$SeqID,] #67 samples, 110139
BS400.taxa[1:5,1:5]
rownames(BS400.taxa) <- c() #explicitly assign nothing to the rownames, so that column_to_rownames has a blank slate to work with...
BS400.taxa <- column_to_rownames(BS400.taxa, var="SeqID")
dim(BS400.taxa) #67, 110138 before filtering / 31, 110138 omitting 3-6cm
# remove OTUs that occur in only 3 or fewer samples (TITAN requirement)
BS400.taxa <- BS400.taxa[colSums(BS400.taxa > 0) >= 3 ] # 6449 after filtering / 3151 filtered for omitted samples
# "BS400.taxa > 0" creates a TRUE/FALSE table for abundance values >0, 
# colSums() adds up all the values in each column (TRUE = 1, FALSE = 0)
# subset the table for where their are more than 3 TRUE (non-zero) values!
dim(BS400.taxa)
# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
# Subset the Metadata400 table to create the .env argument:
BS400.env <- Metadata400[ ,c("SeqID", "numdate")]
head(BS400.env)
BS400.env <- unique(BS400.env)
dim(BS400.env)
BS400.env <- column_to_rownames(BS400.env, var="SeqID")
###check that everything looks right before launching TITAN
# .env argument should be numeric:
class(BS400.env$numdate)
# if it's not numeric, tell it to be numeric:
BS400.env$numdate <- as.numeric(BS400.env$numdate)
# .env and .taxa tables should have equivallent rownames:
head(BS400.env)
BS400.taxa[1:5,1:5]
# both tables should have the same number of rows (first number)
dim(BS400.env)
dim(BS400.taxa)

# Run TITAN2 with default parameters:
# WARNING - it can take hours to run!
#example:: glades.titan <- titan(glades.env, glades.taxa) 
BS400.titan <- titan(BS400.env, BS400.taxa) 
# start @ 1110, end @ 
# start @ 2:59pm, end @ 4:28! after 500 permutations
# MAKE NOTE THE OUTPUTS! ...(so you don't have to run it again!)
#
# "100% occurrence detected 1 times (0.1% of taxa)"
# 
# "Proportion of pure and reliable taxa = 0.0576"

##EXPLORE THE RESULTS OF TITAN!
#
# overview of what's contained within out TITAN object:
summary(BS400.titan)
# Each item within the TITAN object can be accessed as if it were a column
# ...for example: BS400.titan$sumz.cp
#
#             Length  Class  Mode   
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
#output for ASVs:
#       cp      0.05      0.10      0.50    0.90      0.95
#sumz-  17904   17904.0   17904.0   18006   18127.0   18127 ##negatively responding taxa, "17904" = Jan2019 = the first post-fire sample!
#sumz+  18330   18267.0   18286.0   18305   18317.5   18330 ##positively responding taxa, "18330" = March2020
#fsumz- 18024   17988.0   17988.0   18024   18040.5   18057 ##filtered sumz-, "18024" = May2019
#fsumz+ 18305   18191.5   18191.5   18267   18286.0   18305 ##filtered sumz+, "18304" = Feb2020

#output for OTUs:
#       cp      0.05      0.10    0.50    0.90    0.95
#sumz-  18006   17904.00  17904   18006   18127   18127.0
#sumz+  18286   18071.50  18086   18206   18305   18305.0
#fsumz- 18006   17988.00  17988   18006   18024   18040.5
#fsumz+ 18267   18084.55  18127   18239   18267   18286.0


# Table of TITAN results for all taxa
# Zscores, IndVal scores, p-values, etc...
head(BS400.titan$sppmax)
dim(BS400.titan$sppmax) #6449x16 for combined data, 3151x16 for omitted
# SAVE THIS TABLE!
BS400.titan.dt <- as.data.table(BS400.titan$sppmax)
# Results table doesn't have informative rownames, so re-attach ASVs:
# transpose taxa table and then move rownames to the results table
BS400.taxa <- as.data.frame(t(BS400.taxa)) #BS.taxa was the input to TITAN that can also be seen as BS.titan$taxa
head(BS400.taxa)
BS400.taxa.rn <- rownames_to_column(BS400.taxa, var="OTU")
BS400.taxa.rn[1:5, 1:2]
# Paste OTUs onto the TITAN results table:
merge1 <- cbind(BS400.taxa.rn$OTU, BS400.titan.dt)
head(merge1)
names(merge1)[1] <- "GenusOTU"
# split the GenusASV column to make an ASV column that can be matched with the ASV column in tax.df
merge1[ , c("Genus", "OTU"):=tstrsplit(GenusOTU, "_", fixed=TRUE)]
#merge1[ ,c("OTU")] <- grep(pattern = "ASV*", merge1$GenusOTU)


tax.df[1:5, 1:5]
tax.df.rn <- rownames_to_column(tax.df, var="OTU")
tax.df.rn[1:5, 1:5]
# merge taxonomy table with results table:
BS400.titan.AllResults <- merge(merge1, tax.df.rn, by="OTU")
head(BS400.titan.AllResults)
# EXPORT RESULTS to your computer!
write.csv(BS400.titan.AllResults, "~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS_0-3cmONLY/TITAN400_ITS+16S.csv")
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
#
#
#PLOT RESULTS!
#
# Summary plot of total z scores...
plot_sumz_density(BS400.titan)

# illustrate the filtered sumz's...
plot_sumz(BS400.titan, xmin = 17800, xmax = 18400, col1 = "blue", col2 = "orange")
#

#plot z-scores of indicator taxa for the +/- change points! NEW METHOD (colored distributions)
plot_taxa_ridges(BS400.titan, xlabel = "Time")
#plot plot z-scores of indicator taxa for the +/- change points! OLD METHOD (B&W circle+line "boxplots")
plot_taxa(BS400.titan, prob95 = TRUE, xlabel = "Time",  xmin = 17800, xmax = 18400, leg.x = 0.1, leg.y = 10)
#plot Z-score with ggplot2!!
ggplot(BS400.titan.TopResults, aes(x=zenv.cp, y=reorder(GenusOTU, -zenv.cp), color=Phylum))+
  geom_point(aes(size=zscore))+
  geom_errorbar(aes(xmax=`95%`, xmin=`5%`))+
  theme_light()+
  xlab("Time (days)")+
  ylab("Indicator Taxa")+
  ggtitle("High-Intensity TITAN results")


#plot abundance of each of the indicator taxa over time in a faceted manner:
plot_cps(BS400.titan)
#plot an individual taxon from the above plot
plot_cps(BS400.titan, taxaID = "Streptomyces_OTU208", cp.med=FALSE, xlabel = "Time")
#black circles = abundance at each time point (comes directly from input table)
#red circle = observed change point (if "cp.med=TRUE" = median change point across all bootstrap reps)
#blue = histogram of bootstrap replicates by probability density
# IF YOU GET AN ERROR THAT SAYS "FIGURE MARGINS TOO LARGE"... 
## ...clear your plotting window by clicking the broom!
#
#
########### comp321W TITAN! 
#
#
# (1) OTU Table -- Subset OTUtax.t.rn table for only the column names that match the SeqIDs in Metadata400
BS321W.taxa <- OTUtax.t.rn[OTUtax.t.rn$SeqID %in% Metadata321W$SeqID,]
BS321W.taxa[1:5,1:5]
rownames(BS321W.taxa) <- c() #explicitly assign nothing to the rownames, so that column_to_rownames has a blank slate to work with...
BS321W.taxa <- column_to_rownames(BS321W.taxa, var="SeqID")
dim(BS321W.taxa) #46, 110138
# remove OTUs that occur in only 3 or fewer samples (TITAN requirement)
BS321W.taxa <- BS321W.taxa[colSums(BS321W.taxa > 0) >= 3 ]
dim(BS321W.taxa) # 46, 5302
# "BS321W.taxa > 0" creates a TRUE/FALSE table for abundance values >0, 
# colSums() adds up all the values in each column (TRUE = 1, FALSE = 0)
# subset the table for where their are more than 3 TRUE (non-zero) values!

# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
# Subset the Metadata400 table to create the .env argument:
BS321W.env <- Metadata321W[ ,c("SeqID", "numdate")]
head(BS321W.env)
BS321W.env <- unique(BS321W.env)
dim(BS321W.env)
BS321W.env <- column_to_rownames(BS321W.env, var="SeqID")
###check that everything looks right before launching TITAN
# .env argument should be numeric:
class(BS321W.env$numdate)
# if it's not numeric, tell it to be numeric:
BS321W.env$numdate <- as.numeric(BS321W.env$numdate)
# .env and .taxa tables should have equivallent rownames:
head(BS321W.env)
BS321W.taxa[1:5,1:5]
# both tables should have the same number of rows (first number)
dim(BS321W.env)
dim(BS321W.taxa)

# Run TITAN2 with default parameters:
# WARNING - it can take hours to run!
#example:: glades.titan <- titan(glades.env, glades.taxa) 
BS321W.titan <- titan(BS321W.env, BS321W.taxa) 
# start @ 0918
# start @ 2:59pm, end @ 4:28! after 500 permutations
# MAKE NOTE THE OUTPUTS! ...(so you don't have to run it again!)
#
# "100% occurrence detected 1 times (0.1% of taxa)"
# 
# "Proportion of pure and reliable taxa = 0.0576"

##EXPLORE THE RESULTS OF TITAN!
#
# overview of what's contained within out TITAN object:
summary(BS321W.titan)
# Each item within the TITAN object can be accessed as if it were a column
# ...for example: BS321W.titan$sumz.cp
#
#             Length  Class  Mode   
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
BS321W.titan$sumz.cp
#output for ASVs:
#       cp      0.05      0.10      0.50    0.90      0.95
#sumz-  17904   17904.0   17904.0   18006   18127.0   18127 ##negatively responding taxa, "17904" = Jan2019 = the first post-fire sample!
#sumz+  18330   18267.0   18286.0   18305   18317.5   18330 ##positively responding taxa, "18330" = March2020
#fsumz- 18024   17988.0   17988.0   18024   18040.5   18057 ##filtered sumz-, "18024" = May2019
#fsumz+ 18305   18191.5   18191.5   18267   18286.0   18305 ##filtered sumz+, "18304" = Feb2020

#output for OTUs:
#       cp      0.05      0.10    0.50    0.90    0.95
#sumz-  18006   17904.00  17904   18006   18127   18127.0
#sumz+  18286   18071.50  18086   18206   18305   18305.0
#fsumz- 18006   17988.00  17988   18006   18024   18040.5
#fsumz+ 18267   18084.55  18127   18239   18267   18286.0


# Table of TITAN results for all taxa
# Zscores, IndVal scores, p-values, etc...
head(BS321W.titan$sppmax)
dim(BS321W.titan$sppmax) #1401x16 for ASVs, 1093x16 for OTUs
# SAVE THIS TABLE!
BS321W.titan.dt <- as.data.table(BS321W.titan$sppmax)
# Results table doesn't have informative rownames, so re-attach ASVs:
# transpose taxa table and then move rownames to the results table
BS321W.taxa <- as.data.frame(t(BS321W.taxa)) #BS.taxa was the input to TITAN that can also be seen as BS.titan$taxa
head(BS321W.taxa)
BS321W.taxa.rn <- rownames_to_column(BS321W.taxa, var="OTU")
BS321W.taxa.rn[1:5, 1:2]
# Paste OTUs onto the TITAN results table:
merge1 <- cbind(BS321W.taxa.rn$OTU, BS321W.titan.dt)
head(merge1)
names(merge1)[1] <- "GenusOTU"
# split the GenusASV column to make an ASV column that can be matched with the ASV column in tax.df
merge1[ , c("Genus", "OTU"):=tstrsplit(GenusOTU, "_", fixed=TRUE)]

tax.df[1:5, 1:5]
tax.df.rn <- rownames_to_column(tax.df, var="OTU")
tax.df.rn[1:5, 1:5]
# merge taxonomy table with results table:
BS321W.titan.AllResults <- merge(merge1, tax.df.rn, by="OTU")
head(BS321W.titan.AllResults)
# EXPORT RESULTS to your computer!
write.csv(BS321W.titan.AllResults, "~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/TITAN321W_ITS+16S.csv")
# filter the table for only the results that were ultimately plotted = where the value in the "filter" column is >0
BS321W.titan.TopResults <- BS321W.titan.AllResults[filter>0]

colnames(BS321W.titan.TopResults)
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
#
#
#PLOT RESULTS!
#
# Summary plot of total z scores...
plot_sumz_density(BS321W.titan)
# illustrate the filtered sumz's...
plot_sumz(BS321W.titan, xmin = 17800, xmax = 18400, col1 = "blue", col2 = "orange")
#

#plot z-scores of indicator taxa for the +/- change points! NEW METHOD (colored distributions)
plot_taxa_ridges(BS321W.titan, xlabel = "Time")
#plot plot z-scores of indicator taxa for the +/- change points! OLD METHOD (B&W circle+line "boxplots")
plot_taxa(BS321W.titan, xlabel = "Time",  xmin = 17800, xmax = 18400)
#plot Z-score with ggplot2!!
ggplot(BS321W.titan.TopResults, aes(x=zenv.cp, y=reorder(GenusOTU, -zenv.cp), color=Phylum))+
  geom_point(aes(size=zscore))+
  geom_errorbar(aes(xmax=`95%`, xmin=`5%`))+
  theme_light()+
  xlab("Time (days)")+
  ylab("Indicator Taxa")+
  ggtitle("Control2 TITAN results")


#plot abundance of each of the indicator taxa over time in a faceted manner:
plot_cps(BS321W.titan)
#plot an individual taxon from the above plot
#plot_cps(BS400.titan, taxaID = "Streptomyces_OTU208", cp.med=FALSE, xlabel = "Time")
#black circles = abundance at each time point (comes directly from input table)
#red circle = observed change point (if "cp.med=TRUE" = median change point across all bootstrap reps)
#blue = histogram of bootstrap replicates by probability density
# IF YOU GET AN ERROR THAT SAYS "FIGURE MARGINS TOO LARGE"... 
## ...clear your plotting window by clicking the broom!
#
#
########### comp321E TITAN! 
#
#
# (1) OTU Table -- Subset OTUtax.t.rn table for only the column names that match the SeqIDs in Metadata400
BS321E.taxa <- OTUtax.t.rn[OTUtax.t.rn$SeqID %in% Metadata321E$SeqID,]
BS321E.taxa[1:5,1:5]
rownames(BS321E.taxa) <- c() #explicitly assign nothing to the rownames, so that column_to_rownames has a blank slate to work with...
BS321E.taxa <- column_to_rownames(BS321E.taxa, var="SeqID")
dim(BS321E.taxa) #94, 110138
# remove OTUs that occur in only 3 or fewer samples (TITAN requirement)
BS321E.taxa <- BS321E.taxa[colSums(BS321E.taxa > 0) >= 3 ]
dim(BS321E.taxa) # 94, 8432
# "BS321E.taxa > 0" creates a TRUE/FALSE table for abundance values >0, 
# colSums() adds up all the values in each column (TRUE = 1, FALSE = 0)
# subset the table for where their are more than 3 TRUE (non-zero) values!

# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
# Subset the Metadata400 table to create the .env argument:
BS321E.env <- Metadata321E[ ,c("SeqID", "numdate")]
head(BS321E.env)
BS321E.env <- unique(BS321E.env)
dim(BS321E.env)
BS321E.env <- column_to_rownames(BS321E.env, var="SeqID")
###check that everything looks right before launching TITAN
# .env argument should be numeric:
class(BS321E.env$numdate)
# if it's not numeric, tell it to be numeric:
BS321E.env$numdate <- as.numeric(BS321E.env$numdate)
# .env and .taxa tables should have equivallent rownames:
head(BS321E.env)
BS321E.taxa[1:5,1:5]
# both tables should have the same number of rows (first number)
dim(BS321E.env)
dim(BS321E.taxa)

# Run TITAN2 with default parameters:
# WARNING - it can take hours to run!
#example:: glades.titan <- titan(glades.env, glades.taxa) 
BS321E.titan <- titan(BS321E.env, BS321E.taxa) 
# start @ 0918, next day @2129
# start @ 2:59pm, end @ 4:28! after 500 permutations
# MAKE NOTE THE OUTPUTS! ...(so you don't have to run it again!)
#
# "100% occurrence detected 1 times (0.1% of taxa)"
# 
# "Proportion of pure and reliable taxa = 0.0576"

##EXPLORE THE RESULTS OF TITAN!
#
# overview of what's contained within out TITAN object:
summary(BS321E.titan)
# Each item within the TITAN object can be accessed as if it were a column
# ...for example: BS321E.titan$sumz.cp
#
#             Length  Class  Mode   
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
BS321E.titan$sumz.cp
#output for ASVs:
#       cp      0.05      0.10      0.50    0.90      0.95
#sumz-  17904   17904.0   17904.0   18006   18127.0   18127 ##negatively responding taxa, "17904" = Jan2019 = the first post-fire sample!
#sumz+  18330   18267.0   18286.0   18305   18317.5   18330 ##positively responding taxa, "18330" = March2020
#fsumz- 18024   17988.0   17988.0   18024   18040.5   18057 ##filtered sumz-, "18024" = May2019
#fsumz+ 18305   18191.5   18191.5   18267   18286.0   18305 ##filtered sumz+, "18304" = Feb2020

#output for OTUs:
#       cp      0.05      0.10    0.50    0.90    0.95
#sumz-  18006   17904.00  17904   18006   18127   18127.0
#sumz+  18286   18071.50  18086   18206   18305   18305.0
#fsumz- 18006   17988.00  17988   18006   18024   18040.5
#fsumz+ 18267   18084.55  18127   18239   18267   18286.0


# Table of TITAN results for all taxa
# Zscores, IndVal scores, p-values, etc...
head(BS321E.titan$sppmax)
dim(BS321E.titan$sppmax) #1401x16 for ASVs, 1093x16 for OTUs, 8432 x 16 for combined
# SAVE THIS TABLE!
BS321E.titan.dt <- as.data.table(BS321E.titan$sppmax)
# Results table doesn't have informative rownames, so re-attach ASVs:
# transpose taxa table and then move rownames to the results table
BS321E.taxa <- as.data.frame(t(BS321E.taxa)) #BS.taxa was the input to TITAN that can also be seen as BS.titan$taxa
head(BS321E.taxa)
BS321E.taxa.rn <- rownames_to_column(BS321E.taxa, var="OTU")
BS321E.taxa.rn[1:5, 1:2]
# Paste OTUs onto the TITAN results table:
merge1 <- cbind(BS321E.taxa.rn$OTU, BS321E.titan.dt)
head(merge1)
names(merge1)[1] <- "GenusOTU"
# split the GenusASV column to make an ASV column that can be matched with the ASV column in tax.df
merge1[ , c("Genus", "OTU"):=tstrsplit(GenusOTU, "_", fixed=TRUE)]

#tax.df[1:5, 1:5]
#tax.df.rn <- rownames_to_column(tax.df, var="OTU")
tax.df.rn[1:5, 1:5]
# merge taxonomy table with results table:
BS321E.titan.AllResults <- merge(merge1, tax.df.rn, by="OTU")
head(BS321E.titan.AllResults)
# EXPORT RESULTS to your computer!
write.csv(BS321E.titan.AllResults, "~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/TITAN321E_ITS+16S.csv")
# filter the table for only the results that were ultimately plotted = where the value in the "filter" column is >0
BS321E.titan.AllResults <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/TITAN321E_ITS+16S.csv")
BS321E.titan.TopResults <- BS321E.titan.AllResults[filter>0]

colnames(BS321E.titan.TopResults)
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
#
#
#PLOT RESULTS!
#
# Summary plot of total z scores...
plot_sumz_density(BS321E.titan)
# illustrate the filtered sumz's...
plot_sumz(BS321E.titan, xmin = 17800, xmax = 18400, col1 = "blue", col2 = "orange")
#

#plot z-scores of indicator taxa for the +/- change points! NEW METHOD (colored distributions)
plot_taxa_ridges(BS321E.titan, xlabel = "Time")
#plot plot z-scores of indicator taxa for the +/- change points! OLD METHOD (B&W circle+line "boxplots")
plot_taxa(BS321E.titan, xlabel = "Time",  xmin = 17800, xmax = 18400)
#plot Z-score with ggplot2!!
ggplot(BS321E.titan.TopResults, aes(x=zenv.cp, y=reorder(GenusOTU, -zenv.cp), color=Phylum))+
  geom_point(aes(size=zscore))+
  geom_errorbar(aes(xmax=`95%`, xmin=`5%`))+
  theme_light()+
  xlab("Time (days)")+
  ylab("Indicator Taxa")+
  ggtitle("Low-Intensity TITAN results")


#plot abundance of each of the indicator taxa over time in a faceted manner:
plot_cps(BS321E.titan)
#plot an individual taxon from the above plot
#plot_cps(BS400.titan, taxaID = "Streptomyces_OTU208", cp.med=FALSE, xlabel = "Time")
#black circles = abundance at each time point (comes directly from input table)
#red circle = observed change point (if "cp.med=TRUE" = median change point across all bootstrap reps)
#blue = histogram of bootstrap replicates by probability density
# IF YOU GET AN ERROR THAT SAYS "FIGURE MARGINS TOO LARGE"... 
## ...clear your plotting window by clicking the broom!
#
#
########### comp240 TITAN! 
#
#
# (1) OTU Table -- Subset OTUtax.t.rn table for only the column names that match the SeqIDs in Metadata400
BS240.taxa <- OTUtax.t.rn[OTUtax.t.rn$SeqID %in% Metadata240$SeqID,]
BS240.taxa[1:5,1:5]
rownames(BS240.taxa) <- c() #explicitly assign nothing to the rownames, so that column_to_rownames has a blank slate to work with...
BS240.taxa <- column_to_rownames(BS240.taxa, var="SeqID")
dim(BS240.taxa) #48, 110138
# remove OTUs that occur in only 3 or fewer samples (TITAN requirement)
BS240.taxa <- BS240.taxa[colSums(BS240.taxa > 0) >= 3 ]
dim(BS240.taxa) # 48, 5418
# "BS240.taxa > 0" creates a TRUE/FALSE table for abundance values >0, 
# colSums() adds up all the values in each column (TRUE = 1, FALSE = 0)
# subset the table for where their are more than 3 TRUE (non-zero) values!

# (2) Evironmental gradient -- data.frame: rownames = samples, single column gradient variable
# Subset the Metadata400 table to create the .env argument:
BS240.env <- Metadata240[ ,c("SeqID", "numdate")]
head(BS240.env)
BS240.env <- unique(BS240.env)
dim(BS240.env)
BS240.env <- column_to_rownames(BS240.env, var="SeqID")
###check that everything looks right before launching TITAN
# .env argument should be numeric:
class(BS240.env$numdate)
# if it's not numeric, tell it to be numeric:
BS240.env$numdate <- as.numeric(BS240.env$numdate)
# .env and .taxa tables should have equivallent rownames:
head(BS240.env)
BS240.taxa[1:5,1:5]
# both tables should have the same number of rows (first number)
dim(BS240.env)
dim(BS240.taxa)

# Run TITAN2 with default parameters:
# WARNING - it can take hours to run!
#example:: glades.titan <- titan(glades.env, glades.taxa) 
BS240.titan <- titan(BS240.env, BS240.taxa) 
# start @ 0918
# start @ 2:59pm, end @ 4:28! after 500 permutations
# MAKE NOTE THE OUTPUTS! ...(so you don't have to run it again!)
#
# "100% occurrence detected 1 times (0.1% of taxa)"
# 
# "Proportion of pure and reliable taxa = 0.0576"

##EXPLORE THE RESULTS OF TITAN!
#
# overview of what's contained within out TITAN object:
summary(BS240.titan)
# Each item within the TITAN object can be accessed as if it were a column
# ...for example: BS321E.titan$sumz.cp
#
#             Length  Class  Mode   
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
BS240.titan$sumz.cp
#output for ASVs:
#       cp      0.05      0.10      0.50    0.90      0.95
#sumz-  17904   17904.0   17904.0   18006   18127.0   18127 ##negatively responding taxa, "17904" = Jan2019 = the first post-fire sample!
#sumz+  18330   18267.0   18286.0   18305   18317.5   18330 ##positively responding taxa, "18330" = March2020
#fsumz- 18024   17988.0   17988.0   18024   18040.5   18057 ##filtered sumz-, "18024" = May2019
#fsumz+ 18305   18191.5   18191.5   18267   18286.0   18305 ##filtered sumz+, "18304" = Feb2020

#output for OTUs:
#       cp      0.05      0.10    0.50    0.90    0.95
#sumz-  18006   17904.00  17904   18006   18127   18127.0
#sumz+  18286   18071.50  18086   18206   18305   18305.0
#fsumz- 18006   17988.00  17988   18006   18024   18040.5
#fsumz+ 18267   18084.55  18127   18239   18267   18286.0


# Table of TITAN results for all taxa
# Zscores, IndVal scores, p-values, etc...
head(BS240.titan$sppmax)
dim(BS240.titan$sppmax) #1401x16 for ASVs, 1093x16 for OTUs, 5418x16 for cobined data
# SAVE THIS TABLE!
BS240.titan.dt <- as.data.table(BS240.titan$sppmax)
# Results table doesn't have informative rownames, so re-attach ASVs:
# transpose taxa table and then move rownames to the results table
BS240.taxa <- as.data.frame(t(BS240.taxa)) #BS.taxa was the input to TITAN that can also be seen as BS.titan$taxa
head(BS240.taxa)
BS240.taxa.rn <- rownames_to_column(BS240.taxa, var="OTU")
BS240.taxa.rn[1:5, 1:2]
# Paste OTUs onto the TITAN results table:
merge1 <- cbind(BS240.taxa.rn$OTU, BS240.titan.dt)
head(merge1)
names(merge1)[1] <- "GenusOTU"
# split the GenusASV column to make an ASV column that can be matched with the ASV column in tax.df
merge1[ , c("Genus", "OTU"):=tstrsplit(GenusOTU, "_", fixed=TRUE)]

#tax.df[1:5, 1:5]
#tax.df.rn <- rownames_to_column(tax.df, var="OTU")
tax.df.rn[1:5, 1:5]
# merge taxonomy table with results table:
BS240.titan.AllResults <- merge(merge1, tax.df.rn, by="OTU")
head(BS240.titan.AllResults)
# EXPORT RESULTS to your computer!
write.csv(BS240.titan.AllResults, "~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/TITAN240_ITS+16S.csv")
# filter the table for only the results that were ultimately plotted = where the value in the "filter" column is >0
BS240.titan.TopResults <- BS240.titan.AllResults[filter>0]

colnames(BS240.titan.TopResults)
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
#
#
#PLOT RESULTS!
#
# Summary plot of total z scores...
plot_sumz_density(BS240.titan)
# illustrate the filtered sumz's...
plot_sumz(BS240.titan, xmin = 17800, xmax = 18400, col1 = "blue", col2 = "orange")
#

#plot z-scores of indicator taxa for the +/- change points! NEW METHOD (colored distributions)
plot_taxa_ridges(BS240.titan, xlabel = "Time")
#plot plot z-scores of indicator taxa for the +/- change points! OLD METHOD (B&W circle+line "boxplots")
plot_taxa(BS240.titan, xlabel = "Time",  xmin = 17800, xmax = 18400)
#plot Z-score with ggplot2!!
ggplot(BS240.titan.TopResults, aes(x=zenv.cp, y=reorder(GenusOTU, -zenv.cp), color=Phylum))+
  geom_point(aes(size=zscore))+
  geom_errorbar(aes(xmax=`95%`, xmin=`5%`))+
  theme_light()+
  xlab("Time (days)")+
  ylab("Indicator Taxa")+
  ggtitle("Control1 TITAN results")


#plot abundance of each of the indicator taxa over time in a faceted manner:
plot_cps(BS240.titan)
#plot an individual taxon from the above plot
#plot_cps(BS400.titan, taxaID = "Streptomyces_OTU208", cp.med=FALSE, xlabel = "Time")
#black circles = abundance at each time point (comes directly from input table)
#red circle = observed change point (if "cp.med=TRUE" = median change point across all bootstrap reps)
#blue = histogram of bootstrap replicates by probability density
# IF YOU GET AN ERROR THAT SAYS "FIGURE MARGINS TOO LARGE"... 
## ...clear your plotting window by clicking the broom!


############################################
#######TITAN Barplot & Dotplot Summary######
############################################

T400 <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/TITAN400_ITS+16S.csv", header=TRUE)
head(T400)

table (T400$maxgrp) # number of taxa in SumZ- = 3262 / SumZ+ = 3187
table (T400$filter) # 0 = 6328, 1 = 62, 2 = 59

T321W <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/TITAN321W_ITS+16S.csv", header = TRUE)
head(T321W)

table (T321W$maxgrp) # number of taxa in SumZ- = 2950 / SumZ+ = 2352
table (T321W$filter) # 0 = 5229, 1 = 47, 2 = 26

T321E <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/TITAN321E_ITS+16S.csv", header = TRUE)
head(T321E)

table(T321E$maxgrp) # number of taxa in SumZ- = 4185 / SumZ+ = 4247
table(T321E$filter) # 0 = 8276, 1 = 73, 2 = 83

T240 <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/TITAN240_ITS+16S.csv", header = TRUE)
head(T240)

table(T240$maxgrp) # number of taxa in SumZ- = 2871 / SumZ+ = 2547
table(T240$filter) # 0 = 5386, 1 = 28, 2 = 4

##TOTAL TAXA identified as Negative vs. Positive Responders##

#create a datafrome with the SumZ- and SumZ+ values
Plot <- c("C1", "C1", "C2", "C2", "Lo", "Lo", "Hi", "Hi")
SumZ <- c("Neg", "Pos", "Neg", "Pos", "Neg", "Pos", "Neg", "Pos")
TotalTaxa <- c(2871, 2547, 2950, 2352, 4185, 4247, 3262, 3187)

## V2 - taxa identified by TITAN
Plot <- c("C1", "C1", "C2", "C2", "Lo", "Lo", "Hi", "Hi")
SumZ <- c("Neg", "Pos", "Neg", "Pos", "Neg", "Pos", "Neg", "Pos")
TotalTaxa <- c(28, 4, 47, 26, 73, 83, 62, 59)


TAXA <- data.frame(Plot, SumZ, TotalTaxa)

print(TAXA)


ggplot(TAXA)+
  geom_bar(aes(y=TotalTaxa, x=Plot, fill = SumZ), stat = "identity", position = "dodge")+
  scale_fill_manual(values=c("#CC0066", "#3399FF"))+
  theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, color = "black"))+
  xlab("PLOT")+
  ylab("Total Taxa Count")+
  ggtitle("TITAN - Negative vs. Positive Responders")


##PROPORTION of Taxa identified as indicators vs. total taxa in Plot##

# Total TAXA input into TITAN
# 400 = 6449 vs. ID = 121
# Ratios. Z- = 0.5124, Z+ = 0.4876
# 321W = 5302 vs. ID =  73
# Ratios. Z- = 0.6438, Z+ = 0.3562
# 321E = 8432 vs. ID = 156
# Ratios. Z- = 0.4679, Z+ = 0.5321
# 240 = 5418 vs. ID = 32
# Ratios. Z- = 0.875, Z+ = 0.125

R <- data.frame(Plot = c("C1", "C1", "C2", "C2", "Lo", "Lo", "Hi", "Hi"), Response = c("neg", "pos", "neg", "pos", "neg", "pos", "neg", "pos"), TotalTaxa = c(28, 4, 47, 26, 73, 83, 62, 59),
                Date = c("17929", "18267", "17855", "18222.5", "17833", "18206", "17904", "18267"))
print(R)
#NumDate = c("17929", "18267", "17855", "18222.5", "17833", "18206", "17904", "18267")
#Date = c("2/2/19", "1/6/20", "11/20/18", "11/22/19", "10/29/18", "11/6/19", "1/8/19", "1/6/20")
#Date = c("2Feb19", "6Jan20", "20Nov18", "22Nov19", "29Oct18", "6Nov19", "8Jan19", "6Jan20")

ggplot(R, aes(x=Date, y=TotalTaxa, color=Response, fill=Plot))+ 
  geom_point(shape=21, stroke=2, size=6)+
  scale_fill_manual(values=c("#009acd", "#00cd9a", "#ff0102", "#ff7f01"))+
  scale_color_manual(values=c("black", "grey50"))+
  theme_classic()+
  theme(axis.text.x = element_text(hjust=0.5, vjust=0.5, angle=45, color = "black", size=14))+
  theme(axis.text.y = element_text(color = "black", size=14))+
  theme(axis.title = element_text(color = "black", size=14))+
  xlab("Date")+
  ylab("Total Number of Taxa")+
  ggtitle("Proportion of Indicator Response Idenitified by TITAN")


FR <- fread("~/Desktop/Box Sync/Blodgett_16S_Comm/TITAN_16S+ITS/FireResponsiveTaxa.csv", header = TRUE)
head(FR)

FR240 <- merge(T240, FR, by = "Genus", all = FALSE)
table(FR240$filter) # 0 = 254, 1 = 2, 2 = 1
# Total240 1 = 28, 2 = 4

FR321E <- merge(T321E, FR, by = "Genus", all = FALSE)
table(FR321E$filter) # 0 = 493, 1 = 16, 2 = 9
# Total321E 1 = 73, 2 = 83

FR321W <- merge(T321W, FR, by = "Genus", all = FALSE)
table(FR321W$filter) # 0 = 253, 1 = 5, 2 = 1
# Total321W 1 = 47, 2 = 26

FR400 <- merge(T400, FR, by = "Genus", all = FALSE)
table(FR400$filter) # 0 = 441, 1 = 23, 2 = 4
# Total400 1 = 62, 2 = 59


P240 <- c(26, 2, 3, 1)
pie(P240, labels = c("Negative", "Fire Reponsive (-)", "Positive", "Fire Reponsive (+)"), border = "white", col = c("#0B6BFC", "#0B6BFC", "#FC630B", "#FC630B"))

P321E <- c(57, 16, 74, 9)
pie(P321E, labels = c("Negative", "Fire Reponsive (-)", "Positive", "Fire Reponsive (+)"), border = "white", col = c("#0B6BFC", "#0B6BFC", "#FC630B", "#FC630B"))


P321W <- c(42, 5, 25, 1)
pie(P321W, labels = c("Negative", "Fire Reponsive (-)", "Positive", "Fire Reponsive (+)"), border = "white", col = c("#0B6BFC", "#0B6BFC", "#FC630B", "#FC630B"))


P400 <- c(39, 23, 55, 4)
pie(P400, labels = c("Negative", "Fire Reponsive (-)", "Positive", "Fire Reponsive (+)"), border = "white", col = c("#0B6BFC", "#0B6BFC", "#FC630B", "#FC630B"))






