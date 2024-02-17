library(rstudioapi)
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set current dir as WD
# setwd("~/your_working_directory")

# -------------------------------------------------------------
# 1. Packages Imports
# -------------------------------------------------------------
library(dada2); packageVersion("dada2")

# -------------------------------------------------------------
# 2. Config
# -------------------------------------------------------------
gene="16S"  # ITS, 16S ou EF1
source("__config.R")

##forward reads
fnFs <- sort(list.files(rawdata, pattern="_R1.fastq", full.names = TRUE))

##extract the sample names from the forward reads
sample.names2 <- sapply(strsplit(basename(fnFs), "-"), `[`, 2)
sample.names3 <- sapply(strsplit(basename(fnFs), "-"), `[`, 3)
sample.names = paste(sample.names2, sample.names3, sep="")


##Place the filtered sequences in a 'Filtered' sub-folder
filtFs <- file.path(dossier_donnees, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

##Display of graphs representing the quality profiles 
plotQualityProfile(fnFs[1:2]) 


##Place the filtered sequences in a 'Filtered' sub-folder
filtFs <- file.path(dossier_donnees, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names


###########
if (file.exists(file.path(output_results,"out.rds"))) {
  out=readRDS(file.path(output_results,"out.rds"))
} else {
  #16s
  # out <- filterAndTrim(fnFs, filtFs, truncLen=10, trimLeft=270 , maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
  # its
  out <- filterAndTrim(fnFs, filtFs,
                       trimLeft=10,
                       truncLen =275,
                       maxN=0, maxEE=2, 
                       truncQ=2, rm.phix=TRUE, 
                       compress=TRUE, 
                       multithread=TRUE, verbose=FALSE) # On Windows set multithread=FALSE
  
  saveRDS(out, file.path(output_results,"out.rds"))
}


# -------------------------------------------------------------
##Calculation of the error rate ----

set.seed(100)
if (file.exists(file.path(output_results,"errF.rds"))) {
  errF=readRDS(file.path(output_results,"errF.rds"))
} else {
  errF <- learnErrors(filtFs, nbases = 1e8, multithread=TRUE, randomize = TRUE)
  saveRDS(errF,file.path(output_results,"errF.rds"))
}

# -------------------------------------------------------------
# dereplication/Infer sequence variants 

dadaFs <- vector("list", length(sample.names))
names(dadaFs) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepFs <- derepFastq(filtFs[[sam]])
  dadaFs[[sam]] <- dada(derepFs, err=errF, multithread=TRUE)
}

saveRDS(dadaFs, file.path(output_results, "dadaFs.rds"))

# -------------------------------------------------------------
# Construct sequence table and write to disk ----

if (file.exists(file.path(output_results,"seqtab.rds"))) {
  seqtab=readRDS(file.path(output_results,"seqtab.rds"))
} else {
  seqtab <- makeSequenceTable(dadaFs)
  saveRDS(seqtab, file.path(output_results,"seqtab.rds"))
}

# -------------------------------------------------------------
# Remove chimeras ----

if (file.exists(file.path(output_results,"seqtab_nochim.rds"))) {
  seqtab.nochim=readRDS(file.path(file.path(output_results,"seqtab_nochim.rds")))
} else {
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=FALSE)
  saveRDS(seqtab.nochim,file.path(output_results,"seqtab_nochim.rds"))
}

# -------------------------------------------------------------
# Assign taxonomy

if (file.exists(file.path(output_results,"taxa_species.rds"))) {
  taxa.species=readRDS(file.path(output_results,"taxa_species.rds"))
} else {
  taxa.species <- assignTaxonomy(seqtab.nochim, file.path(dbs[gene]), multithread=TRUE)
  saveRDS(taxa.species, file.path(output_results,"taxa_species.rds"))
}

# track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
track = as.data.frame(track)
write.csv2(track, file.path(output_results,"track.csv"))

# -------------------------------------------------------------