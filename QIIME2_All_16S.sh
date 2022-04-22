#!/bin/bash
# Job name:
#SBATCH --job-name=qiime2
#
# Account:
#SBATCH --account=fc_traxspec
#
# Partition:
#SBATCH --partition=savio2
#
# Request one node:
#SBATCH --nodes=1
#
# Specify one task:
#SBATCH --ntasks-per-node=1
#
# Number of processors for single task needed for use case (example):
#SBATCH --cpus-per-task=8
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
#SBATCH --output=qiime2_All_16S_%j.out
#SBATCH --error=qiime2_All_16S_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=neempatel@berkeley.edu
#
## Command(s) to run:
mkdir /global/scratch/users/neempatel/tmp
export TMPDIR=/global/scratch/users/neempatel/tmp
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

##Start QIIME2 environment:
module load qiime2/2020.11

##CasavaOneEightSingleLanePerSampleDirFmt - the format of data that comes off the illumina platform. 

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format CasavaOneEightSingleLanePerSampleDirFmt --input-path /global/scratch/users/neempatel/BS_16S_Raw/All_16S/Data/ --output-path /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.qza

##qiime demux summarize --help

##to visualize the data
qiime demux summarize --i-data /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.qza --p-n 10000 --o-visualization BS_16S_Raw/All_16S/QIIMEobjects/demux.qzv


##dada2 integrates quality scores into its algorithm. 


qiime cutadapt trim-paired --i-demultiplexed-sequences /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.qza --p-cores 8 --p-front-f 'CCTACGGGNBGCASCAG' --p-front-r 'GACTACNVGGGTATCTAATCC' --p-match-adapter-wildcards --p-minimum-length 250 --p-discard-untrimmed --o-trimmed-sequences /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.trimmed.qza --verbose

##visualize trimmed data
qiime demux summarize --i-data /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.trimmed.qza --p-n 10000 --o-visualization /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.trimmed.qzv


##Make an OTU table


##To denoise and generate the OTU table:
qiime dada2 denoise-paired --i-demultiplexed-seqs /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.trimmed.qza --p-trunc-len-f 270 --p-trunc-len-r 200 --p-trim-left-f 0 --p-trim-left-r 0 --p-chimera-method consensus --p-pooling-method independent --p-n-reads-learn 1000000 --o-table /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.table.qza --o-representative-sequences /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.rep-seqs.qza --o-denoising-stats /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.stats.qza --verbose


##Turning the stats into qzv objects for visualization:
qiime feature-table summarize --i-table /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.table.qza --o-visualization /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.table.qzv

qiime feature-table tabulate-seqs --i-data /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.rep-seqs.qza --o-visualization /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.rep-seqs.qzv

qiime metadata tabulate --m-input-file /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.stats.qza --o-visualization /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.stats.qzv



##Taxonomy classification/assignment:
##qiime feature-classifier —-help
##classify-sklearn


qiime feature-classifier classify-sklearn --i-classifier /global/scratch/users/neempatel/Databases/silva-138-99-nb-classifier.qza --i-reads /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.rep-seqs.qza --o-classification /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.taxonomy.qza

##To visualize:
qiime metadata tabulate --m-input-file /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.taxonomy.qza --o-visualization /global/scratch/users/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.taxonomy.qzv 



##To make a taxa barplot:
##qiime taxa barplot --i-table /global/scratch/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.table.qza --i-taxonomy /global/scratch/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/demux.taxonomy.qza --m-metadata-file /global/scratch/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/stats-metadata.tsv --o-visualization /global/scratch/neempatel/BS_16S_Raw/All_16S/QIIMEobjects/taxa-p-bar-plots.qzv

