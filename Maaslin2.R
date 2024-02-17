###Maaslin2
BiocManager::install("Maaslin2")
library(Maaslin2)
library(magrittr)
library(tidyverse)

ps_february ##see data_modification.R

ps <- ps_february
ps
ps <- phyloseq(tax_table(ps), otu_table(ps), sample_data(ps))

# keep only taxa with positive sums
ps1 <- prune_taxa(taxa_sums(ps) > 0, ps)

######Filter taxa that are found only in 5% of the particpants
# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps1),
               MARGIN = ifelse(taxa_are_rows(ps1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps1),
                    tax_table(ps1))

# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps1)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
psFilt = prune_taxa(keepTaxa, ps1)

saveRDS(psFilt, "psFilt.rds")

#calculate relative abundance
ps_relative <- transform_sample_counts(psFilt, function(OTU) OTU/sum(OTU) * 100)

psFilt  ##will be used for downstream analyses

#Filter the baseline and endline samples
psFilt_base <- subset_samples(ps_relative, Baseline_Endline == "Baseline")
psFilt_end <- subset_samples(ps_relative, Baseline_Endline == "Endline")

#Agglomerate taxa of the same type

psFilt_baseG <- tax_glom(psFilt_base, taxrank = "Genus")

#save relative frequency per sample as .tsv file
write.table(psFilt_baseG %>% psmelt() %>% 
              select(OTU, Sample, Kingdom, Phylum, Class, Order, Family, Genus, Abundance) %>% pivot_wider(names_from = Sample, values_from = Abundance) %>% t(),
            file = "ps.relative_abundance.genus.tsv", sep = "\t", quote = F, row.names = T, col.names = F)

#save sample data from the psFilt_baseG
baseline_metadata <- sample_data(psFilt_baseG)
write.csv(baseline_metadata, "baseline_metadata.csv")

#input file 
baseline_relative <- ps.relative_abundance.genus #transposed data frame
baseline_stunting <- baseline_metadata %>% drop_na(Stunting)

###Maaslin2
maas_stunting_3 <- Maaslin2(
  input_data = baseline_relative,
  input_metadata = baseline_stunting,
  min_prevalence = 0.0,
  min_abundance = 0.0,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  standardize = "FALSE",
  random_effects = c('Sex', 'Age_groups'),
  fixed_effects = c('Stunting'),
  output = "~/final_project_forisca/maaslin2_stunting_2")


