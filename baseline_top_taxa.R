## load library
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(microViz)
library(tidyverse)
library(magrittr)
tidyverse_update(recursive = FALSE, repos = getOption("repos"))

ps

ps_metadata <- data.frame(sample_data(ps))
write.csv(ps_metadata, "ps_metadata.csv")
######Filter taxa that are found only in 5% of the population
ps_filter <- tax_filter(ps, min_prevalence = 0.05)   ##38/760 samples
sum(taxa_sums(ps_filter) < 100)

saveRDS(ps_filter, "ps_filter.RDS")

####ps_Base
####Extract only the baseline data
psBase_filter <- subset_samples(ps_filter, Baseline_Endline == "Baseline")
psBase_filter
sum(taxa_sums(ps_filter) == 0)  ##0

microbiome::summarize_phyloseq(ps_filter)

####PHYLUM RELATIVE ABUNDANCE AT THE THE BASELINE
###Before we can plot phylum relative abundance, we need to merge all ASV's together that are within the same Phylum:

##Baseline filtered data
df_psBase_filter <- data.frame(sample_data(psBase_filter))

write.csv(df_psBase_filter, "df_psBase_filter.csv")

##Calculate relative abundance for the baseline 
psBasefilter_relabun <- 

##read the responder and non-responder group
library(readr)
responder_non_responder_zinc <- read_csv("responder_non_responder_zinc.csv")
View(responder_non_responder_zinc)

##Read the responder and non-responder
responder_non_responder_df <- sample_data(responder_non_responder_zinc)
sample_names(responder_non_responder_df) <- responder_non_responder_df$sampleid

##create a new phylose object for response
psBase_response <- phyloseq(otu_table(psBase_filter), tax_table(psBase_filter),
                            sample_data(responder_non_responder_df))

##complete object 

library(ANCOMBC)
library(mia)
##Select responders and non-responder groups
psBase_resfinal <- subset_samples(psBase_response, select_response == "1")
psBase_resfinal <- phyloseq::prune_taxa(taxa_sums(psBase_resfinal) > 0, psBase_resfinal)  ##all of them are above zero
sum(taxa_sums(psBase_resfinal) == 0)

##assuming you have a phyloseq object called physeq
tax_table(psBase_resfinal) <- tax_table(psBase_resfinal)[,1:6]

psBase_resfinalGE <- tax_glom(psBase_resfinal, taxrank = "Genus")

##change to tree summarized experiment
tse_psBase_resfinal <- mia::makeTreeSummarizedExperimentFromPhyloseq(psBase_resfinalGE)

# Run ANCOM-BC at the genus level and only including the prevalent genera
ancombc_response1SD <- ancombc2(data = tse_psBase_resfinal,
                               assay_name = "counts",
                               fix_formula = "Response_1SD",
                               tax_level = "Genus",
                               p_adj_method = "holm",
                               group = "Response_1SD",
                               struc_zero = TRUE,
                               neg_lb = TRUE,
                               # multi group comparison is deactivated automatically
                               global = TRUE)

ancombc_1SD <- ancombc_response1SD$res
ancombc_1SD_1 <- ancombc_1SD %>%
  dplyr::select(taxon, `lfc_Response_1SDResponder`,`q_Response_1SDResponder`, se_Response_1SDResponder) %>%
  filter(`q_Response_1SDResponder` < 0.05)

psBase_resfinalGE

###add phylum information 
tax_resfinalGE <- tax_table(psBase_resfinalGE)
tax_resfinalGE <- data.frame(tax_resfinalGE)
tax_resfinalGE1 <- tax_resfinalGE %>% select(-Species)

response_baseline_all <- left_join(ancombc_1SD_1, tax_resfinalGE1, by = c('taxon' = 'Genus')) 

response_baseline_all_1 <- response_baseline_all %>% dplyr::mutate(Direction = ifelse(lfc_Response_1SDResponder > 0, "Positive", "Negative"))

#Add signifcance star
response_baseline_all_1 <- response_baseline_all_1 %>%
  mutate(Significance = ifelse(q_Response_1SDResponder <= 0.001, '***',
                               ifelse(q_Response_1SDResponder <= 0.01, '**',
                                      ifelse(q_Response_1SDResponder <= 0.05, '*', ''))))
View(response_baseline_all_1)

##plot fold change at the baseline between responder and non-responder groups
response_baseline_all_1 %>%
  arrange(lfc_Response_1SDResponder) %>% mutate(label_y = ifelse(lfc_Response_1SDResponder < 0, 0.08, -0.08),
                                                    label_hjust = ifelse(lfc_Response_1SDResponder < 0, 0, 1)) %>% 
  ggplot(aes(x = lfc_Response_1SDResponder, y = reorder(taxon, lfc_Response_1SDResponder), fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_errorbar(aes(xmin = lfc_Response_1SDResponder - se_Response_1SDResponder, 
                    xmax = lfc_Response_1SDResponder + se_Response_1SDResponder), width = 0.2,
                position = position_dodge(0.05)) +
  theme_void() +
  theme(legend.text = element_text(size = 11, face = "italic"), legend.title = element_text(size = 10, face = "bold"), legend.title.align = 0.5, legend.margin = margin(t = 0, b = 0, r = 0, l = -30),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 12, face = "bold"),
        axis.ticks.y = element_blank(), text = element_text(size = 2),
        panel.grid.major.x = element_line(colour = "gray", linewidth = 0.4, linetype = "dotted"),
        plot.margin = margin(l = 10, r = 5, b = 10, t = 10, unit = "pt")) +
  geom_text(aes(x = label_y, label = taxon, hjust = label_hjust), size = 3.5) + #label text based on value
  labs(x = "Log fold change", y = "") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4))
ggsave("Responder_non-responder_final/baseline_response.tiff", bg = "white", dpi = 600)

psBase_resfinalGE

# assuming you have a phyloseq object called physeq
tax_table(psBase_resfinalGE) <- tax_table(psBase_resfinalGE)[,1:6]

##remove the species rank from the phyloseq object
tax_table(ps_responseall) <- tax_table(ps_responseall)[,1:6]

#glom taxa at the genus level
ps_responseallGE <- tax_glom(ps_responseall, taxrank = "Genus")

sum(taxa_sums(ps_responseall) == 0)

library(Maaslin2)
longitudinal_response <- Maaslin2(
  input_data = data.frame(otu_table(ps_responseallGE)),
  input_metadata = data.frame(sample_data(ps_responseallGE)),
  output = "~/final_project_forisca/Responder_non-responder_final",
  min_abundance = 0.0,
  min_prevalence = 0.0,
  normalization = "CSS",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  fixed_effects = c("Response_1SD", "Baseline_Endline", "Intervention_groups"),
  random_effects = c("SubjectID"),
  correction = "BH",
  standardize = FALSE,
  cores = 1)

longitudinal_responseNEGBIN <- Maaslin2(
  input_data = data.frame(otu_table(ps_responseallGE)),
  input_metadata = data.frame(sample_data(ps_responseallGE)),
  output = "~/final_project_forisca/Responder_non-responder_final",
  min_abundance = 0.0,
  min_prevalence = 0.0,
  normalization = "TSS",
  transform = "NONE",
  analysis_method = "NEGBIN",
  max_significance = 0.05,
  fixed_effects = c("Response_1SD", "Baseline_Endline", "Intervention_groups"),
  random_effects = c("SubjectID"),
  correction = "BH",
  standardize = FALSE,
  cores = 1)

##Baseline genus relative abundance
psBaseresponse_rel <- transform_sample_counts(psBase_response, function(OTU) OTU/sum(OTU))

##remove the species rank from the phyloseq object
tax_table(psBaseresponse_rel) <- tax_table(psBaseresponse_rel)[,1:6]

psBasefilter_relabunGE <- tax_glom(psBaseresponse_rel, taxrank = "Genus")

# Create a vector of taxon names to filter on
taxa_base_response <- response_baseline_all_1$taxon
taxa_base_response <- as.character(taxa_base_response)

# Use subset_samples to create a new phyloseq object with matching taxa
filtered_base_response <- subset_taxa(psBasefilter_relabunGE, Genus %in% taxa_base_response)

filtered_base_response_table <- psmelt(filtered_base_response)
write.csv(filtered_base_response_table, "filtered_base_response_table.csv")

##Plot the top genus
filtered_base_response_table %>% 
  ggplot(aes(x = reorder(Genus, -Abundance),
             y = Abundance,
             fill = Response_1SD)) +
  geom_boxplot() + 
  labs(x = "", y = "Relative Abundance (%)") +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2", name = "Response") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        strip.text = element_text(size = 11, face = "bold"),
        axis.title.y = element_text(size = 11, face = "bold"), 
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(), axis.line = element_line(linewidth = 0.5, colour = "black"), legend.title.align = 0.5,
        legend.margin = margin(t = 0, b = 0, r = 0, l = -10))
ggsave("area_plot/week24_community_typing.tiff", dpi = 600)

###############################################################################
BiocManager::install("lefser")
library(lefser)

res <- lefser(zeller14, groupCol = "study_condition", blockCol = "age_category")

##Lefse microbiome
sum(taxa_sums(psBase_resfinal) == 0)  ##0

taxa_are_rows(psBase_resfinal)  ##FALSE

psBase_resfinalT <- phyloseq(t(otu_table(psBase_resfinal)), sample_data(psBase_resfinal), tax_table(psBase_resfinal))
taxa_are_rows(psBase_resfinalT) ##TRUE 

tree_resfinalT <- mia::makeTreeSummarizedExperimentFromPhyloseq(psBase_resfinalT)

##Factorize
tree_resfinalT$Response_1SD <- factor(tree_resfinalT$Response_1SD, levels = c("Responder", "Non-responder"))

##Relative abundance
tree_resfinalTRB <- relativeAb(tree_resfinalT)
tree_resfinalTRB <- relativeAb(tree_resfinalT)

##Run LEfSe on the responder
res <- lefser(tree_resfinalTRB, groupCol = "Response_1SD")

##remove species rank from the phyloseq object
tax_table(psBase_resfinal) <- tax_table(psBase_resfinal)[, 1:6]
#LEfSe microbiome 
mm_lefse <- run_lefse(
  psBase_resfinalT,
  wilcoxon_cutoff = 0.01,
  norm = "CPM",
  group = "Response_1SD",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 2
)

##After intervention what happened
library(readr)
responder_non_responder_zinc <- read_csv("responder_non_responder_zinc.csv")
View(responder_non_responder_zinc)

responder_non_responder_zinc1 <- sample_data(responder_non_responder_zinc)

sample_names(responder_non_responder_zinc1) <- responder_non_responder_zinc1$sampleid

##create a new phylose object for response
psfilter_response <- phyloseq(otu_table(ps_filter), tax_table(ps_filter),
                            sample_data(responder_non_responder_zinc1))
response_base_end

ps_responseall <- subset_samples(psfilter_response, response_base_end == "1")

ps_responsend <- subset_samples(psfilter_response, Baseline_Endline == "Endline")

response_base_end

ps_responsend1 <- subset_samples(ps_responsend, response_base_end == "1")

sum(taxa_sums(ps_responsend1) == 0)  ##12

ps_responsend1 <- phyloseq::prune_taxa(taxa_sums(ps_responsend1) > 0, ps_responsend1)  ##all of them are above zero

ps_responsend1GE <- tax_glom(ps_responsend1, taxrank = "Genus")

##change to tree summarized experiment
tse_ps_responsend1 <- mia::makeTreeSummarizedExperimentFromPhyloseq(ps_responsend1GE)

# Run ANCOM-BC at the genus level and only including the prevalent genera
ancombc_respend1 <- ancombc2(data = tse_ps_responsend1,
                               assay_name = "counts",
                               fix_formula = "Response_1SD",
                               tax_level = "Genus",
                               p_adj_method = "holm",
                               group = "Response_1SD",
                               struc_zero = TRUE,
                               neg_lb = TRUE,
                               # multi group comparison is deactivated automatically
                               global = TRUE)

ancombc_1endSD <- ancombc_respend1$res
ancombc_1endSD_1 <- ancombc_1endSD %>%
  dplyr::select(taxon, `lfc_Response_1SDResponder`,`q_Response_1SDResponder`, se_Response_1SDResponder) %>%
  filter(`q_Response_1SDResponder` < 0.05)

###add phylum information 
tax_end_resfinalGE <- tax_table(ps_responsend1GE)
tax_end_resfinalGE <- data.frame(tax_end_resfinalGE)
tax_end_resfinalGE1 <- tax_end_resfinalGE %>% select(-Species)

response_endline_all <- left_join(ancombc_1endSD_1, tax_end_resfinalGE1, by = c('taxon' = 'Genus')) 

response_endline_all_1 <- response_endline_all %>% dplyr::mutate(Direction = ifelse(lfc_Response_1SDResponder > 0, "Positive", "Negative"))

#Add signifcance star
response_endline_all_1 <- response_endline_all_1 %>%
  mutate(Significance = ifelse(q_Response_1SDResponder <= 0.001, '***',
                               ifelse(q_Response_1SDResponder <= 0.01, '**',
                                      ifelse(q_Response_1SDResponder <= 0.05, '*', ''))))
View(response_endline_all_1)

##plot fold change at the baseline between responder and non-responder groups
response_endline_all_1 %>%
  arrange(lfc_Response_1SDResponder) %>% mutate(label_y = ifelse(lfc_Response_1SDResponder < 0, 0.08, -0.08),
                                                label_hjust = ifelse(lfc_Response_1SDResponder < 0, 0, 1)) %>% 
  ggplot(aes(x = lfc_Response_1SDResponder, y = reorder(taxon, lfc_Response_1SDResponder), fill = Phylum)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_errorbar(aes(xmin = lfc_Response_1SDResponder - se_Response_1SDResponder, 
                    xmax = lfc_Response_1SDResponder + se_Response_1SDResponder), width = 0.2,
                position = position_dodge(0.05)) +
  theme_void() +
  theme(legend.text = element_text(size = 11, face = "italic"), legend.title = element_text(size = 10, face = "bold"), legend.title.align = 0.5, legend.margin = margin(t = 0, b = 0, r = 0, l = -10),
        axis.text.y = element_blank(), axis.title.y = element_blank(), axis.text.x = element_text(size = 11), axis.title.x = element_text(size = 12, face = "bold"),
        axis.ticks.y = element_blank(), text = element_text(size = 2),
        panel.grid.major.x = element_line(colour = "gray", linewidth = 0.4, linetype = "dotted"),
        plot.margin = margin(l = 10, r = 5, b = 10, t = 10, unit = "pt")) +
  geom_text(aes(x = label_y, label = taxon, hjust = label_hjust), size = 3.5) + #label text based on value
  labs(x = "Log fold change", y = "") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2, 3, 4))
ggsave("Responder_non-responder_final/endline_response.tiff", bg = "white", dpi = 600)


#Agglomerate taxa
psBase_resfinal_genus <- phyloseq::tax_glom(psBase_resfinal, "Genus")

phyloseq::taxa_names(psBase_resfinal_genus) <- phyloseq::tax_table(psBase_resfinal_genus)[, "Genus"]
#Run selbal
cv_sebal <- selbal::selbal.cv(x = data.frame(t(data.frame(phyloseq::otu_table(psBase_resfinal_genus)))), 
                              y = phyloseq::sample_data(psBase_resfinal_genus)$Status, 
                              n.fold = 5, n.iter = 1) 

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