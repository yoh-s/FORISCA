library("phyloseq")
library("dada2")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library(data.table)
library("tibble")       # Needed for converting column to row names
library("readxl")
library("DECIPHER")
library("phangorn")
library("msa")
library("microbiome")
library("ggpubr")
library("magrittr")

# -------------------------------------------------------------
# 2. Config
# -------------------------------------------------------------
gene="16S"  # ITS, 16S ou EF1
source("__config.R")

seqtab.nochim <- seqtab_nochim
taxa.species <- taxa_species

#combine data into a phyloseq object
samdf <- read_excel("~/your_working_directory/project_forisca_metadata_test.xlsx", sheet = "forisca_metadata")
samdf <- samdf %>% tibble::column_to_rownames("sampleID")

##create a phyloseq object
psoriginal <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(samdf), tax_table(taxa.species))
psoriginal
saveRDS(psoriginal, file.path(output_results,"psoriginal.rds"))

#Subsetting the data (keep=1)
##dataset in keep =1 has complete metadata value
#keep = 0 shows several missing or no metadata value, thus will be removed from further analysis
pssubset = subset_samples(psoriginal, keep != "0")
pssubset
saveRDS(pssubset, file.path(output_results,"pssubset.rds"))

#remove 0 abundance
sorted_taxasums =  sort(taxa_sums(pssubset), TRUE)
non_zero_taxa = sorted_taxasums[sorted_taxasums!=0]
names(non_zero_taxa)
psabun = prune_taxa(names(non_zero_taxa), pssubset)
psabun
saveRDS(psabun, file.path(output_results, "psabun.rds"))

#prevalence calculation
prev = apply(X=otu_table(pssubset),
             MARGIN=ifelse(taxa_are_rows(pssubset), yes=1, no=2),
             FUN=function(x){sum(x>0)})
prev

prevdf=data.frame(Prevalence=prev,
                  TotalAbundance=taxa_sums(pssubset),
                  tax_table(pssubset))
prevdf

abund = apply(X=otu_table(pssubset),
              MARGIN=ifelse(taxa_are_rows(pssubset), yes=1, no=2),
              FUN=function(x){sum(x)})
abund

# execute prevalence filter using prune_taxa()
ps_prevabun =  prune_taxa((prev>0) & (abund>10), pssubset)
ps_prevabun
saveRDS(ps_prevabun, file.path(output_results, "ps_prevabun.rds"))
sequences_prevabun = taxa_names(ps_prevabun)

# alternatively skip if the above step is followed
sequences_prevabun =getSequences(tax_table(ps_prevabun))


####################tree_from_sequence function#######################
tree_from_sequences <- function(seqs, seqtab) {
  # seqs = sequences2
  
  names(seqs) <- seqs
  
  library("phangorn")
  library("msa")
  library("dada2")
  
  # multiple alignment
  mult <- msaClustalOmega(seqs, type="dna", order="input")
  
  # creation of phylogenetic trees
  phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
  dm<-dist.ml(phang.align)
  
  treeNJ <- NJ(dm)
  
  
  fit = pml(treeNJ, data=phang.align)
  
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", opInv=TRUE, optGamma=TRUE, rearrangement="stochastic", control=pml.control(trace=0))
  
  return(fitGTR$tree)
  
}

##################end of tree_from_sequence function_____________________________

#Run the above function first and follow the next step to construct phylogenetic tree

# creation and insert tree for this subset
tree = tree_from_sequences(sequences_prevabun, seqtab.nochim)

# save tree in nwk format
write.tree(tree, file = "treefile", append = FALSE,
           digits = 10, tree.names = FALSE)

# add tree to phyloseq object
prevabun_tree <- merge_phyloseq(ps_prevabun , tree)

#save phyloseq object
saveRDS(prevabun_tree, file.path(output_results, "prevabun_tree.rds"))
# plot trees

#create table, number of features for each phyla
table(tax_table(prevabun_tree) [, "Phylum"], exclude = NULL)

prevabun_tree0 <- subset_taxa(prevabun_tree, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
saveRDS(prevabun_tree0, file.path(output_results, "prevabun_tree0.rds"))

###Rooting a phylogenetic tree from the phyloseq object
###############function to pick_new_outgroup#########################
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}
###############end of the function################################

#rooting the phylogenetic tree tree but difficult to implement
new.outgroup = pick_new_outgroup(tree)
rootedTree = ape::root(tree, outgroup=new.outgroup, resolve.root=TRUE)
rootedTree ###check if it is rooted

#save rooted tree
write.tree(rootedTree, file = "rootedTree", append = FALSE,
           digits = 10, tree.names = FALSE)

##################################################################

#create a phyloseq object with rooted tree and do the rarefaction, null phylum removal
# add tree to phyloseq object
prevabun_rooted <- merge_phyloseq(ps_prevabun , rootedTree)

saveRDS(prevabun_rooted, file.path(output_results, "prevabun_rooted.rds"))

#create table, number of features for each phyla
table(tax_table(prevabun_rooted) [, "Phylum"], exclude = NULL)

prevabun_tree0_rooted <- subset_taxa(prevabun_rooted, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
saveRDS(prevabun_tree0_rooted, file.path(output_results, "prevabun_tree0_rooted.rds"))

################################################################################