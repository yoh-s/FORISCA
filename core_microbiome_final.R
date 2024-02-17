###prevabun_tree0_rooted is subseted to keep only taxa with positive sums

##prevabun_tree is original data 
library(microbiome)

ps <- prevabun_tree0_rooted

# keep only taxa with positive sums
ps <- prune_taxa(taxa_sums(prevabun_tree0_rooted) > 0, prevabun_tree0_rooted)
print(ps)

####Extract only the baseline data
ps_Base <- subset_samples(ps, Baseline_Endline=="Baseline")
ps_Base

ps_End <- subset_samples(ps, Baseline_Endline == "Endline")
ps_End

# Calculate compositional version of the data
# (relative abundances)
psBase.rel <- microbiome::transform(ps_Base, "compositional")
psEnd.rel <- microbiome::transform(ps_End, "compositional")
###Outputs of deblur/dada2 will most likely have seqs as rownames instead of OTU ids or taxa names

taxa_names(psBase.rel)[1:3]

##################################################################

psBase.relG <-tax_glom(psBase.rel, taxrank = "Genus")
psEnd.relG <- tax_glom(psEnd.rel, taxrank = "Genus")

##caluclate the core taxa

psge95 <- core(psBase.relG, detection = 0, prevalence = 95/100)
psend95 <- core(psEnd.relG, detection = 0, prevalence = 95/100)

####################change to data frame
corebase95 <- psmelt(psge95)
coreend95 <- psmelt(psend95)

###save data frame as csv

write.table(corebase95, "corebase95t.csv", row.names = FALSE)
write.table(coreend95, "coreend95.csv", row.names = FALSE)

##export to GraphPad prism for plotting a graph

#Baseline core microbiota based on intervention groups, melt phyloseq
#export to GraphPad prism for plotting a graph
# Create a separate phyloseq object for each intervention group at the baseline
psgeplacebo95 <- core(subset_samples(psBase.relG, Intervention_groups == "Placebo"), detection = 0, prevalence = 95/100)
corebaseplacebo95 <- psmelt(psgeplacebo95)
write.table(corebaseplacebo95, "core_march/corebaseplacebo95.csv", row.names = FALSE)

psgeoriginal95 <- core(subset_samples(psBase.relG, Intervention_groups == "UR Original"), detection = 0, prevalence = 95/100)
corebaseoriginal95 <- psmelt(psgeoriginal95)
write.table(corebaseoriginal95, "core_march/corebaseoriginal95.csv", row.names = FALSE)

psgeimproved95 <- core(subset_samples(psBase.relG, Intervention_groups == "UR Improved"), detection = 0, prevalence = 95/100)

coreimproved95 <- psmelt(psgeimproved95)
write.table(coreimproved95, "core_march/coreimproved95.csv", row.names = FALSE)

#Endline core microbiota based on intervention groups, melt phyloseq
#export to GraphPad prism for plotting a graph
psEnd.relG

psenplacebo95 <- core(subset_samples(psEnd.relG, Intervention_groups == "Placebo"), detection = 0, prevalence = 95/100)
coreendplacebo95 <- psmelt(psenplacebo95)
write.table(coreendplacebo95, "core_march/coreendplacebo95.csv", row.names = FALSE)

psenoriginal95 <- core(subset_samples(psEnd.relG, Intervention_groups == "UR Original"), detection = 0, prevalence = 95/100)
coreendoriginal95 <- psmelt(psenoriginal95)
write.table(coreendoriginal95, "core_march/coreendoriginal95.csv", row.names = FALSE)

psenimproved95 <- core(subset_samples(psEnd.relG, Intervention_groups == "UR Improved"), detection = 0, prevalence = 95/100)
coreendimproved95 <- psmelt(psenimproved95)
write.table(coreendimproved95, "core_march/coreendimproved95.csv", row.names = FALSE)
