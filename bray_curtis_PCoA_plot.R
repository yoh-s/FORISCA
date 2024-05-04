library(Biostrings)
library(magrittr)
library(vegan)
library(phyloseq.extended)
dna <- Biostrings::DNAStringSet(taxa_names(prevabun_tree0_rooted))
names(dna) <- taxa_names(prevabun_tree0_rooted)
prevabun_tree0_rooted <- merge_phyloseq(prevabun_tree0_rooted, dna)
taxa_names(prevabun_tree0_rooted) <- paste0("ASV", seq(ntaxa(prevabun_tree0_rooted)))
# now check again
prevabun_tree0_rooted
# save the modified phyloseq object to a new name
# look for the new name of the phyloseq object "prevabun_treeasv0.rds"

saveRDS(prevabun_tree0_rooted, file = "prevabun_treeasv0.rds")

###Check the metadata in the prevabun_tree0_rooted to see the number of variables and to try to create
##a new metadata for the analysis of Procrustes analysis
##to check the 

prevabuntreeorooted_metadata <- sample_data(prevabun_tree0_rooted)
is.data.frame(prevabuntreeorooted_metadata)
write.csv(prevabuntreeorooted_metadata, file = "metadata_feb_2023/prevabuntree0rootedme.csv")

###go to procrustes analysis where the change on the metadata was made

# modify the sample data to include new variables and remove unwanted ones.
# find the number of reads where the sequences were rarefied 
# rarefied at 3500 reads per sample

#############note from previous work#######################
View(tax_table(prevabun_tree0_rooted))
ps <- prevabun_tree0_rooted
microbiome::summarize_phyloseq(prevabun_tree0_rooted)
##Average number of reads = 39355.6
##Median number of reads = 35378.5
################################################################################
# Assuming you have a phyloseq object 'prevabun_tree0_rooted'
mean_sample_depth <- 39355.6
cutoff <- round(0.001 * mean_sample_depth)

# Calculate total counts for each ASV
total_counts <- rowSums(otu_table(prevabun_tree0_rooted))

# Identify ASVs with total counts greater than or equal to the cutoff
asvs_to_keep <- names(total_counts[total_counts >= cutoff])

# Filter the phyloseq object to keep only the selected ASVs
prevabun_tree0_rooted_filtered <- prune_samples(asvs_to_keep, prevabun_tree0_rooted)

# Subset the OTU table to keep only the selected ASVs
otu_table_filtered <- otu_table(prevabun_tree0_rooted)[asvs_to_keep, ]

# Create a new phyloseq object with the filtered OTU table
prevabun_tree0_rooted_filtered <- phyloseq(
  otu_table(otu_table_filtered, taxa_are_rows = TRUE),
  sample_data(prevabun_tree0_rooted),
  tax_table(prevabun_tree0_rooted)
)
###lot beta diversity - unweighted unifrac
bray_dist = phyloseq::distance(ps.rarefied, method="bray")
bray_ordination = ordinate(ps.rarefied, method="PCoA", distance=bray_dist)

##Extract metadata from the rarefied phyloseq object
ps.rarefied_df <- data.frame(sample_data(ps.rarefied))
write.csv(ps.rarefied_df, "ps.rarefied_df.csv")

# Define ggplot2 custom colors:
allGroupColors <- c("brown3", "blue3")

plot_ordination(ps.rarefied, bray_ordination, color="Baseline_Endline") + theme(aspect.ratio=1) +
  theme_classic() +
  geom_point(size = 2) +
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13), axis.line = element_line(size = 1, linetype = "solid"),
        axis.title.x = element_text(size=13.5, face="bold"),
        axis.title.y = element_text(size=13.5, face="bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12)) +
  labs(x = "PCoA 1 (14.2%)", y = "PCoA 2 (5.7%)") +
  stat_ellipse(aes(group = Baseline_Endline)) +
  scale_color_manual(values = allGroupColors) +
  guides(color = guide_legend(title = "   "))

ggsave('bray_curtis-final-sepetmber.tiff', dpi=300)

##Test significance 
adonis2(bray_dist ~ Baseline_Endline, data = ps.rarefied_df, strata = ps.rarefied_df$SubjectID)

