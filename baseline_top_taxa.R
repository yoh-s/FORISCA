## load library
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(magrittr)
tidyverse_update(recursive = FALSE, repos = getOption("repos"))

#Agglomerate taxa
psBase_resfinal_genus <- phyloseq::tax_glom(psBase_resfinal, "Genus")

phyloseq::taxa_names(psBase_resfinal_genus) <- phyloseq::tax_table(psBase_resfinal_genus)[, "Genus"]
psBasefilter_relabun <- transform_sample_counts(psBase_filter, function(OTU) OTU/sum(OTU) * 100)

##Selbal

df_psBasefilter_relabun <- psmelt(psBasefilter_relabun)
write.csv(df_psBasefilter_relabun, "df_psBasefilter_relabun.csv")

psBase_phylum <- tax_glom(psBasefilter_relabun, "Phylum", NArm = TRUE)

##Change the phyloseq object to a dataframe

taxa_abundance_table_phylum <- psmelt(psBase_phylum)

BoxPlot_phylumBASE <- taxa_abundance_table_phylum %>% 
  ggplot(aes(x =Abundance, y = reorder(Phylum, +Abundance), fill = Phylum)) +
  geom_boxplot(lwd=1, outlier.shape = NA, alpha=0.5, width=0.6,coef=0, show.legend = FALSE) +
  geom_jitter(show.legend = FALSE, shape=21, color="black", position=position_jitter(width=0.15, height=0.15)) +
  labs(x = "Relative Abundance (%)",
       y = "",
       title = "Baseline Phylum Relative Abundance") +
  theme(
    axis.text.x = element_text(size = 11, face = "bold"),
    plot.title=element_text(hjust=0.5),
    axis.text.y = element_text(size = 11),
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 11, face = "bold")
  ) +
  theme_classic()
BoxPlot_phylumBASE
ggsave('BoxPlot_phylumBASE.tiff', dpi=600)

##
taxa_abundance_table_phylum %>% 
  ggplot(aes(x =Abundance, y = reorder(Phylum, +Abundance), fill = Phylum)) +
  geom_boxplot(lwd=1, outlier.shape = NA, alpha=0.5, width=0.6,coef=0, show.legend = FALSE) +
  geom_jitter(show.legend = FALSE, shape=21, color="black", position=position_jitter(width=0.15, height=0.15)) +
  labs(x = "Relative Abundance (%)",
       y = "") +
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    plot.title=element_text(hjust=0.5),
    axis.text.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 10, face = "bold")
  )

BoxPlot_phylumBASE
ggsave('base10/BoxPlot_phylumBASE10.tiff', dpi=600)

#############CLASS RELATIVE ABUNDANCE

psBase_class <- tax_glom(psBaseFi, "Class", NArm = TRUE)

# Transform Taxa counts to relative abundance
psBase_class_relabun <- transform_sample_counts(psBase_class, function(OTU) OTU/sum(OTU) * 100)

# Convert into dataframe
taxa_abundance_table_class <- psmelt(psBase_class_relabun)

###And now we can plot a box plot: ordered top to bottom b the abundance

BoxPlot_classBASE <- taxa_abundance_table_class %>% 
  ggplot(aes(x =Abundance, y = reorder(Class, +Abundance), fill = Class)) +
  geom_boxplot(lwd=1, outlier.shape = NA, alpha=0.5, width=0.6,coef=0) +
  geom_jitter(show.legend = FALSE, shape=21, color="black", position=position_jitter(width=0.15, height=0.15)) +
  labs(x = "Relative Abundance (%)",
       y = "",
       title = "Baseline Class Relative Abundance") +
  theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    plot.title=element_text(hjust=0.5),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

BoxPlot_classBASE + theme(legend.position = "none")
ggsave('BoxPlot_classBASE.tiff', dpi=600)

########
taxa_abundance_table_phylum %>% 
  ggplot(aes(x =Abundance, y = reorder(Phylum, +Abundance), fill = Phylum)) +
  geom_boxplot(lwd=1, outlier.shape = NA, alpha=0.5, width=0.6,coef=0, show.legend = FALSE) +
  geom_jitter(show.legend = FALSE, shape=21, color="black", position=position_jitter(width=0.15, height=0.15)) +
  labs(x = "Relative Abundance (%)",
       y = "") +
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    plot.title=element_text(hjust=0.5),
    axis.text.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 10, face = "bold")
  )

####ORDER RELATIVE ABUNDANCE AT THE THE BASELINE

# Merge everything to the phylum level
psBase_order <- tax_glom(psBaseFi, "Order", NArm = TRUE)

# Transform Taxa counts to relative abundance
ps1_order_relabun <- transform_sample_counts(psBase_order, function(OTU) OTU/sum(OTU) * 100)

###Then extract the data from the phyloseq object:

taxa_abundance_table_order <- psmelt(ps1_order_relabun)

###And now we can plot a box plot: ordered top to bottom b the abundance

BoxPlot_orderBASE <- taxa_abundance_table_order %>% 
  ggplot(aes(x =Abundance, y = reorder(Order, +Abundance), fill = Order)) +
  geom_boxplot(lwd=1, outlier.shape = NA, alpha=0.5, width=0.6,coef=0) +
  geom_jitter(show.legend = FALSE, shape=21, color="black", position=position_jitter(width=0.15, height=0.15)) +
  labs(x = "Relative Abundance (%)",
       y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1, face = "bold"),
    plot.title=element_text(hjust=0.5),
    axis.text.y = element_text(size = 10,face = "bold"),
    strip.text = element_text(size = 12)
  )

BoxPlot_orderBASE + theme(legend.position = "none")
ggsave('BoxPlot_orderALLBASE.tiff', dpi=600)



taxa_abundance_table_phylum %>% 
  ggplot(aes(x =Abundance, y = reorder(Phylum, +Abundance), fill = Phylum)) +
  geom_boxplot(lwd=1, outlier.shape = NA, alpha=0.5, width=0.6,coef=0, show.legend = FALSE) +
  geom_jitter(show.legend = FALSE, shape=21, color="black", position=position_jitter(width=0.15, height=0.15)) +
  labs(x = "Relative Abundance (%)",
       y = "") +
  theme_classic() + 
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),
    plot.title=element_text(hjust=0.5),
    axis.text.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 10, face = "bold")
  )


##############FAMILY RELATIVE ABUNDANCE (Top 20)

ps1_family <- tax_glom(psBaseFi, "Family", NArm = TRUE)

# Get top 20 genera
top20_family <- names(sort(taxa_sums(ps1_family), decreasing=TRUE))[1:20]

# Extract the top 10 taxa and Regular Diet Samples
ps1_family_top20 <- prune_taxa(top20_family, ps1_family_relabun)

# Convert into dataframe
taxa_abundance_table_family20 <- psmelt(ps1_family_top20)

BoxPlot_familyBASE20 <- taxa_abundance_table_family20 %>% 
  ggplot(aes(x =Abundance, y = reorder(Family, +Abundance), fill = Family)) +
  geom_boxplot(lwd=1, outlier.shape = NA, alpha=0.5, width=0.6,coef=0) +
  geom_jitter(show.legend = FALSE, shape=21, color="black", position=position_jitter(width=0.15, height=0.15)) +
  labs(x = "Relative Abundance (%)",
       y = "",
       title = "Baseline top 20 family Relative Abundance") +
  theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    plot.title=element_text(hjust=0.5),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

BoxPlot_familyBASE20 + theme(legend.position = "none")
ggsave('BoxPlot_familyBASE20.tiff', dpi=600)

##########################GENUS RELATIVE ABUNDANCE (Top 20)
ps1_genus <- tax_glom(psBaseFi, "Genus", NArm = TRUE)

#Get top 20 genera
top20_genera <- names(sort(taxa_sums(ps1_genus), decreasing=TRUE))[1:20]

# Extract the top 10 taxa and Regular Diet Samples
ps1_genus_top20 <- prune_taxa(top20_genera, ps1_genus_relabun)

# Convert into dataframe
taxa_abundance_table_genus20 <- psmelt(ps1_genus_top20)

BoxPlot_genusBASE20 <- taxa_abundance_table_genus20 %>% 
  ggplot(aes(x =Abundance, y = reorder(Genus, +Abundance), fill = Genus)) +
  geom_boxplot(lwd=1, outlier.shape = NA, alpha=0.5, width=0.6,coef=0) +
  geom_jitter(show.legend = FALSE, shape=21, color="black", position=position_jitter(width=0.15, height=0.15)) +
  labs(x = "Relative Abundance (%)",
       y = "",
       title = "Baseline top 20 genus Relative Abundance") +
  theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    plot.title=element_text(hjust=0.5),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12)
  )

BoxPlot_genusBASE20 + theme(legend.position = "none")
ggsave('BoxPlot_genusBASE20.tiff', dpi=600)
