### work on foresca
###We can change the name of the sequence to ASVIDs

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
################################################################################
sum(taxa_sums(prevabun_tree0_rooted) < 5) ##0
sum(taxa_sums(prevabun_tree0_rooted) < 11) ##0  ##it was zero in the original phyloseq object
sum(taxa_sums(prevabun_tree0_rooted) < 21) ###4496
sum(taxa_sums(prevabun_tree0_rooted) < 31) ###6229
sum(taxa_sums(prevabun_tree0_rooted) < 41) ###7149
sum(taxa_sums(prevabun_tree0_rooted) < 101) ##9567
ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=3500, replace=F)

######################Plot beta diversity - unweighted unifrac
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

##summarize phyloseq object
prevabun_tree0_rooted  ##14032
microbiome::summarize_phyloseq(prevabun_tree0_rooted)
#Sparsity = 0.988821228170197
#Min 2948 and Max 314227

sum(taxa_sums(prevabun_tree0_rooted) < 50)
View(tax_table(prevabun_tree0_rooted))

ps.rarefied
microbiome::summarize_phyloseq(ps.rarefied)
##Number of singletons 2058
library(remotes)
remotes::install_github("mahendra-mariadassou/phyloseq-extended", ref = "dev")

##Linear mixed model
##Extract baseline
psrarefied_base

##Alpha diversity
library(phyloseq.extended)

base_rareedge <- edge_pca(psrarefied_base)

##alpha diversity all
baseline_alpha <- microbiome::alpha(psrarefied_base, index = "all")

##Phylogenetic diversity
baseline_faithpd <- phylodiv(psrarefied_base)

##change rownames to sampleid2
baseline_alpha <- baseline_alpha %>% rownames_to_column(var = "sampleid2")

baseline_alphdiversity_all <- right_join(baseline_faithpd, baseline_alpha, by = "sampleid2")

##Factorize
baseline_alphdiversity_all$Age_groups <- factor(baseline_alphdiversity_all$Age_groups)
baseline_alphdiversity_all$Anemia_class <- factor(baseline_alphdiversity_all$Anemia_class)
baseline_alphdiversity_all$Vitamin_A_deficiency <- factor(baseline_alphdiversity_all$Vitamin_A_deficiency)
baseline_alphdiversity_all$Iron_deficiency_anemia <- factor(baseline_alphdiversity_all$Iron_deficiency_anemia)
baseline_alphdiversity_all$Iron_deficiency <- factor(baseline_alphdiversity_all$Iron_deficiency)
baseline_alphdiversity_all$Systemic_inflammation <- factor(baseline_alphdiversity_all$Systemic_inflammation)
baseline_alphdiversity_all$Zinc_deficiency <- factor(baseline_alphdiversity_all$Zinc_deficiency)
baseline_alphdiversity_all$Hemoglobinopathy <- factor(baseline_alphdiversity_all$Hemoglobinopathy)
baseline_alphdiversity_all$Parasite_infection <- factor(baseline_alphdiversity_all$Parasite_infection)
baseline_alphdiversity_all$Stunting <- factor(baseline_alphdiversity_all$Stunting)
baseline_alphdiversity_all$Underweight <- factor(baseline_alphdiversity_all$Underweight)
baseline_alphdiversity_all$Gastrointestinal_inflammation <- factor(baseline_alphdiversity_all$Gastrointestinal_inflammation)

####Age
##evenness_pielou  Pielou's evenness
##diversity_shannon  Shannon diversity
##Baseline alpha diversity (pairwise difference)
lm(data = baseline_alphdiversity_all ~ evenness_pielou+Sex)
age_evenness <- lm(evenness_pielou ~ Age_groups + Sex, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_age = cbind(coefficients(summary(age_evenness)), confint(age_evenness))
#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    0.661|      0.013|  51.612|              0.000|  0.635|  0.686|
#|Age_groups6-9_years |   -0.039|      0.015|  -2.657|              0.008| -0.068| -0.010|
#|SexGirl             |   -0.008|      0.015|  -0.570|              0.569| -0.037|  0.020|

##The regression suggests that age group "6 - 9 years" has a statistically significant negative effect on the dependent variable,
##while being a girl does not have a statistically significant effect. The intercept represents the estimated value of the
##dependent variable when all predictors are zero or at their reference levels.

##Anemia
anemia_evenness <- lm(evenness_pielou ~ Anemia_class + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_anemia = cbind(coefficients(summary(anemia_evenness)), confint(anemia_evenness))
knitr::kable(res_evenness_anemia, digits = 3)

#|                      | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:---------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)           |    0.622|      0.019|  32.238|              0.000|  0.585|  0.660|
#|Anemia_classNo anemia |    0.047|      0.018|   2.616|              0.009|  0.012|  0.083|
#|SexGirl               |   -0.011|      0.014|  -0.750|              0.454| -0.039|  0.018|
#|Age_groups6-9_years   |   -0.036|      0.015|  -2.501|              0.013| -0.065| -0.008|

##Vitamin A deficiency
vitA_evenness <- lm(evenness_pielou ~ Vitamin_A_deficiency + Age_groups + Sex, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_vitA = cbind(coefficients(summary(vitA_evenness)), confint(vitA_evenness))
knitr::kable(res_evenness_vitA, digits = 3)

#|                        | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-----------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)             |    0.661|      0.013|  51.733|              0.000|  0.636|  0.686|
#|Vitamin_A_deficiencyYes |   -0.043|      0.027|  -1.566|              0.118| -0.097|  0.011|
#|Age_groups6-9_years     |   -0.034|      0.015|  -2.289|              0.023| -0.063| -0.005|
#|SexGirl                 |   -0.008|      0.015|  -0.522|              0.602| -0.036|  0.021|

##Iron deficiency anemia
ironda_evenness <- lm(evenness_pielou ~ Iron_deficiency_anemia + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_ironda = cbind(coefficients(summary(ironda_evenness)), confint(ironda_evenness))
knitr::kable(res_evenness_ironda, digits = 3)

#|                          | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)               |    0.668|      0.013|  51.190|              0.000|  0.643|  0.694|
#|Iron_deficiency_anemiaYes |   -0.057|      0.022|  -2.615|              0.009| -0.101| -0.014|
#|SexGirl                   |   -0.012|      0.014|  -0.795|              0.427| -0.040|  0.017|
#|Age_groups6-9_years       |   -0.038|      0.015|  -2.591|              0.010| -0.066| -0.009|

##Iron_deficiency
irond_evenness <- lm(evenness_pielou ~ Iron_deficiency + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_id = cbind(coefficients(summary(irond_evenness)), confint(irond_evenness))
knitr::kable(res_evenness_id, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    0.672|      0.015|  44.427|              0.000|  0.643|  0.702|
#|Iron_deficiencyYes  |   -0.021|      0.015|  -1.470|              0.142| -0.050|  0.007|
#|SexGirl             |   -0.010|      0.015|  -0.668|              0.504| -0.038|  0.019|
#|Age_groups6-9_years |   -0.040|      0.015|  -2.707|              0.007| -0.068| -0.011|
##Systemic inflammation
sysinfla_evenness <- lm(evenness_pielou ~ Systemic_inflammation + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_sysinfla = cbind(coefficients(summary(sysinfla_evenness)), confint(sysinfla_evenness))
knitr::kable(res_evenness_sysinfla, digits = 3)

#|                         | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:------------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)              |    0.661|      0.014|  48.617|              0.000|  0.634|  0.688|
#|Systemic_inflammationYes |   -0.002|      0.015|  -0.141|              0.888| -0.032|  0.028|
#|SexGirl                  |   -0.008|      0.015|  -0.556|              0.579| -0.037|  0.021|
#|Age_groups6-9_years      |   -0.039|      0.015|  -2.641|              0.009| -0.068| -0.010|

##Zinc_deficiency
zincde_evenness <- lm(evenness_pielou ~ Zinc_deficiency + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_zincda = cbind(coefficients(summary(zincde_evenness)), confint(zincde_evenness))
knitr::kable(res_evenness_zincda, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    0.646|      0.027|  24.318|              0.000|  0.594|  0.699|
#|Zinc_deficiencyYes  |    0.007|      0.025|   0.291|              0.771| -0.041|  0.055|
#|SexGirl             |    0.009|      0.015|   0.588|              0.557| -0.021|  0.039|
#|Age_groups6-9_years |   -0.030|      0.015|  -1.946|              0.053| -0.060|  0.000|

##Hemoglobinopathy
hemoglobin_evenness <- lm(evenness_pielou ~ Hemoglobinopathy + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_hemoglo = cbind(coefficients(summary(hemoglobin_evenness)), confint(hemoglobin_evenness))
knitr::kable(res_evenness_hemoglo, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    0.674|      0.015|  43.529|              0.000|  0.643|  0.704|
#|HemoglobinopathyYes |   -0.022|      0.015|  -1.512|              0.131| -0.051|  0.007|
#|SexGirl             |   -0.007|      0.015|  -0.502|              0.616| -0.036|  0.021|
#|Age_groups6-9_years |   -0.041|      0.015|  -2.811|              0.005| -0.070| -0.012|

##Parasite_infection
parainfec_evenness <- lm(evenness_pielou ~ Parasite_infection + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_parainfec = cbind(coefficients(summary(parainfec_evenness)), confint(parainfec_evenness))
knitr::kable(res_evenness_parainfec, digits = 3)

#|                      | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:---------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)           |    0.656|      0.014|  48.171|              0.000|  0.629|  0.683|
#|Parasite_infectionYes |    0.016|      0.016|   1.007|              0.315| -0.016|  0.049|
#|SexGirl               |   -0.009|      0.015|  -0.627|              0.531| -0.038|  0.020|
#|Age_groups6-9_years   |   -0.038|      0.015|  -2.580|              0.010| -0.067| -0.009|

##Stunting
stunt_evenness <- lm(evenness_pielou ~ Stunting + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_stunt = cbind(coefficients(summary(stunt_evenness)), confint(stunt_evenness))
knitr::kable(res_evenness_stunt, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    0.672|      0.016|  42.674|              0.000|  0.642|  0.703|
#|StuntingYes         |   -0.020|      0.015|  -1.300|              0.194| -0.049|  0.010|
#|SexGirl             |   -0.008|      0.015|  -0.584|              0.560| -0.037|  0.020|
#|Age_groups6-9_years |   -0.044|      0.015|  -2.912|              0.004| -0.074| -0.014|

#Underweight
under_evenness <- lm(evenness_pielou ~ Underweight + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_under = cbind(coefficients(summary(under_evenness)), confint(under_evenness))
knitr::kable(res_evenness_under, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    0.775|      0.076|  10.189|              0.000|  0.625|  0.925|
#|UnderweightYes      |   -0.038|      0.020|  -1.854|              0.065| -0.078|  0.002|
#|SexGirl             |   -0.020|      0.020|  -0.970|              0.333| -0.060|  0.020|
#|Age_groups6-9_years |   -0.130|      0.076|  -1.727|              0.086| -0.279|  0.018|

##Hemoglobinopathy
hemoglo_evenness <- lm(evenness_pielou ~ Hemoglobinopathy + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_hemoglo = cbind(coefficients(summary(hemoglo_evenness)), confint(hemoglo_evenness))
knitr::kable(res_evenness_hemoglo, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    0.674|      0.015|  43.529|              0.000|  0.643|  0.704|
#|HemoglobinopathyYes |   -0.022|      0.015|  -1.512|              0.131| -0.051|  0.007|
#|SexGirl             |   -0.007|      0.015|  -0.502|              0.616| -0.036|  0.021|
#|Age_groups6-9_years |   -0.041|      0.015|  -2.811|              0.005| -0.070| -0.012|

##Gastrointestinal_inflammation
gastroinf_evenness <- lm(evenness_pielou ~ Gastrointestinal_inflammation + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_evenness_gastroinf = cbind(coefficients(summary(gastroinf_evenness)), confint(gastroinf_evenness))
knitr::kable(res_evenness_gastroinf, digits = 3)

#|                                 | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:--------------------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)                      |    0.661|      0.014|  47.691|              0.000|  0.634|  0.689|
#|Gastrointestinal_inflammationYes |   -0.005|      0.027|  -0.172|              0.864| -0.058|  0.049|
#|SexGirl                          |   -0.010|      0.016|  -0.647|              0.518| -0.041|  0.021|
#|Age_groups6-9_years              |   -0.042|      0.016|  -2.717|              0.007| -0.073| -0.012|

##Shannon diversity
#diversity_shannon
age_shannon <- lm(diversity_shannon ~ Age_groups + Sex, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_shannon_age = cbind(coefficients(summary(age_shannon)), confint(age_shannon))
knitr::kable(res_shannon_age, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|(Intercept)         |    3.280|      0.075|  43.823|              0.000|  3.133|  3.427|
#|Age_groups6-9_years |   -0.216|      0.086|  -2.521|              0.012| -0.384| -0.047|
#|SexGirl             |   -0.079|      0.085|  -0.926|              0.355| -0.246|  0.089|

vitA_shannon <- lm(diversity_shannon ~ Vitamin_A_deficiency + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_shannon_vita = cbind(coefficients(summary(vitA_shannon)), confint(vitA_shannon))
knitr::kable(res_shannon_vita, digits = 3)

#|                        | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-----------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)             |    3.282|      0.075|  43.916|              0.000|  3.135|  3.429|
#|Vitamin_A_deficiencyYes |   -0.239|      0.161|  -1.488|              0.137| -0.555|  0.077|
#|SexGirl                 |   -0.075|      0.085|  -0.881|              0.379| -0.242|  0.092|
#|Age_groups6-9_years     |   -0.189|      0.087|  -2.171|              0.031| -0.361| -0.018|

##Iron deficiency anemia
irond_shannon <- lm(diversity_shannon ~ Iron_deficiency_anemia + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_shannon_irond = cbind(coefficients(summary(irond_shannon)), confint(irond_shannon))
knitr::kable(res_shannon_irond, digits = 3)

#|                          | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)               |    3.322|      0.076|  43.437|              0.000|  3.172|  3.473|
#|Iron_deficiency_anemiaYes |   -0.307|      0.129|  -2.382|              0.018| -0.560| -0.054|
#|SexGirl                   |   -0.096|      0.085|  -1.131|              0.259| -0.263|  0.071|
#|Age_groups6-9_years       |   -0.209|      0.085|  -2.457|              0.014| -0.376| -0.042|

##Anemia
##Iron deficiency anemia
anemia_shannon <- lm(diversity_shannon ~ Anemia_class + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_shannon_anemia = cbind(coefficients(summary(anemia_shannon)), confint(anemia_shannon))
knitr::kable(res_shannon_anemia, digits = 3)

#|                      | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:---------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)           |    3.082|      0.113|  27.240|              0.000|  2.859|  3.304|
#|Anemia_classNo anemia |    0.247|      0.106|   2.326|              0.021|  0.038|  0.455|
#|SexGirl               |   -0.092|      0.085|  -1.086|              0.278| -0.259|  0.075|
#|Age_groups6-9_years   |   -0.203|      0.085|  -2.378|              0.018| -0.370| -0.035|

##Iron_deficiency
irond_shannon <- lm(diversity_shannon ~ Iron_deficiency + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_shannon_irond = cbind(coefficients(summary(irond_shannon)), confint(irond_shannon))
knitr::kable(res_shannon_irond, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    3.348|      0.089|  37.817|              0.000|  3.174|  3.522|
#|Iron_deficiencyYes  |   -0.122|      0.085|  -1.434|              0.153| -0.289|  0.045|
#|SexGirl             |   -0.087|      0.085|  -1.022|              0.308| -0.255|  0.080|
#|Age_groups6-9_years |   -0.219|      0.085|  -2.569|              0.011| -0.387| -0.051|

##Systemic_inflammation
sysinfla_shannon <- lm(diversity_shannon ~ Systemic_inflammation + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_sysinfla_irond = cbind(coefficients(summary(sysinfla_shannon)), confint(sysinfla_shannon))
knitr::kable(res_sysinfla_irond, digits = 3)

#|                         | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:------------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)              |    3.288|      0.080|  41.352|              0.000|  3.132|  3.445|
#|Systemic_inflammationYes |   -0.028|      0.089|  -0.321|              0.749| -0.203|  0.146|
#|SexGirl                  |   -0.077|      0.086|  -0.895|              0.371| -0.245|  0.092|
#|Age_groups6-9_years      |   -0.214|      0.086|  -2.496|              0.013| -0.383| -0.045|

##Zinc_deficiency
zincd_shannon <- lm(diversity_shannon ~ Zinc_deficiency + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_zincd_irond = cbind(coefficients(summary(zincd_shannon)), confint(zincd_shannon))
knitr::kable(res_zincd_irond, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    3.253|      0.157|  20.756|              0.000|  2.945|  3.561|
#|Zinc_deficiencyYes  |   -0.029|      0.145|  -0.198|              0.843| -0.313|  0.256|
#|SexGirl             |    0.013|      0.089|   0.147|              0.883| -0.163|  0.189|
#|Age_groups6-9_years |   -0.152|      0.090|  -1.694|              0.091| -0.329|  0.025|

##Hemoglobinopathy
hemoglo_shannon <- lm(diversity_shannon ~ Hemoglobinopathy + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_hemoglo_irond = cbind(coefficients(summary(hemoglo_shannon)), confint(hemoglo_shannon))
knitr::kable(res_hemoglo_irond, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    3.357|      0.091|  37.083|              0.000|  3.179|  3.535|
#|HemoglobinopathyYes |   -0.129|      0.086|  -1.505|              0.133| -0.298|  0.040|
#|SexGirl             |   -0.073|      0.085|  -0.858|              0.391| -0.240|  0.094|
#|Age_groups6-9_years |   -0.230|      0.086|  -2.674|              0.008| -0.399| -0.061|

##Parasite_infection
parasinfe_shannon <- lm(diversity_shannon ~ Parasite_infection + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_parasinfe_irond = cbind(coefficients(summary(parasinfe_shannon)), confint(parasinfe_shannon))
knitr::kable(res_parasinfe_irond, digits = 3)

#|                      | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:---------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)           |    3.249|      0.080|  40.824|              0.000|  3.093|  3.406|
#|Parasite_infectionYes |    0.107|      0.096|   1.122|              0.262| -0.081|  0.295|
#|SexGirl               |   -0.084|      0.085|  -0.990|              0.323| -0.252|  0.083|
#|Age_groups6-9_years   |   -0.209|      0.086|  -2.437|              0.015| -0.377| -0.040|

##Stunting
stunt_shannon <- lm(diversity_shannon ~ Stunting + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_stunt_irond = cbind(coefficients(summary(stunt_shannon)), confint(stunt_shannon))
knitr::kable(res_stunt_irond, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    3.314|      0.092|  35.900|              0.000|  3.133|  3.496|
#|StuntingYes         |   -0.056|      0.089|  -0.636|              0.525| -0.231|  0.118|
#|SexGirl             |   -0.079|      0.085|  -0.932|              0.352| -0.247|  0.088|
#|Age_groups6-9_years |   -0.231|      0.089|  -2.597|              0.010| -0.406| -0.056|

##Underweight
underwe_shannon <- lm(diversity_shannon ~ Underweight + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_underw_irond = cbind(coefficients(summary(underwe_shannon)), confint(underwe_shannon))
knitr::kable(res_underw_irond, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:-------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)         |    3.998|      0.438|   9.135|              0.000|  3.135|  4.860|
#|UnderweightYes      |   -0.178|      0.118|  -1.516|              0.131| -0.410|  0.054|
#|SexGirl             |   -0.146|      0.118|  -1.242|              0.216| -0.378|  0.086|
#|Age_groups6-9_years |   -0.817|      0.434|  -1.882|              0.061| -1.673|  0.039|

##Gastrointestinal_inflammation
gastroin_shannon <- lm(diversity_shannon ~ Gastrointestinal_inflammation + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_gastroin_irond = cbind(coefficients(summary(gastroin_shannon)), confint(gastroin_shannon))
knitr::kable(res_gastroin_irond, digits = 3)

#|                                 | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:--------------------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)                      |    3.286|      0.081|  40.726|              0.000|  3.128|  3.445|
#|Gastrointestinal_inflammationYes |   -0.079|      0.158|  -0.499|              0.618| -0.391|  0.233|
#|SexGirl                          |   -0.091|      0.091|  -1.004|              0.316| -0.269|  0.087|
#|Age_groups6-9_years              |   -0.237|      0.091|  -2.599|              0.010| -0.416| -0.058|

####Faith phylogenetic diversity

#Iron_deficiency
#Systemic_inflammation
#Zinc_deficiency
#Hemoglobinopathy
#Parasite_infection
#Stunting
#Underweight
#Gastrointestinal_inflammation

##Age_groups
age_faithpd <- lm(pd ~ Age_groups + Sex, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_age_faithpd = cbind(coefficients(summary(age_faithpd)), confint(age_faithpd))
knitr::kable(res_age_faithpd, digits = 3)

#|                    | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|(Intercept)         |   23.627|      0.606|  39.008|              0.000| 22.436| 24.818|
#|Age_groups6-9_years |   -0.863|      0.692|  -1.247|              0.213| -2.224|  0.498|
#|SexGirl             |   -0.617|      0.689|  -0.896|              0.371| -1.972|  0.737|

#Anemia_class
anemia_faithpd <- lm(pd ~ Anemia_class + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_anemia_faithpd = cbind(coefficients(summary(anemia_faithpd)), confint(anemia_faithpd))
knitr::kable(res_anemia_faithpd, digits = 3)

#|                      | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|:---------------------|--------:|----------:|-------:|------------------:|------:|------:|
#|(Intercept)           |   23.654|      0.922|  25.651|              0.000| 21.840| 25.467|
#|Anemia_classNo anemia |   -0.032|      0.864|  -0.038|              0.970| -1.732|  1.667|
#|SexGirl               |   -0.616|      0.691|  -0.891|              0.374| -1.975|  0.744|
#|Age_groups6-9_years   |   -0.865|      0.694|  -1.245|              0.214| -2.230|  0.501|

#Vitamin_A_deficiency
vitaA_faithpd <- lm(pd ~ Vitamin_A_deficiency + Sex + Age_groups, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_vitA_faithpd = cbind(coefficients(summary(vitaA_faithpd)), confint(vitaA_faithpd))
knitr::kable(res_vitA_faithpd, digits = 3)
#|                        | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|(Intercept)             |   23.224|      0.483|  48.092|              0.000| 22.275| 24.174|
#|Vitamin_A_deficiencyYes |   -0.963|      1.278|  -0.754|              0.451| -3.475|  1.549|
#|SexGirl                 |   -0.618|      0.690|  -0.896|              0.371| -1.975|  0.739|

#Iron_deficiency_anemia
irond_faithpd <- lm(pd ~ Iron_deficiency_anemia + Sex, baseline_alphdiversity_all)

#Combine the summary of the model with the 95% confidence interval of the estimates
res_irond_faithpd = cbind(coefficients(summary(irond_faithpd)), confint(irond_faithpd))
knitr::kable(res_irond_faithpd, digits = 3)

#|                          | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|  2.5 %| 97.5 %|
#|(Intercept)               |   23.263|      0.500|  46.519|              0.000| 22.279| 24.246|
#|Iron_deficiency_anemiaYes |   -0.705|      1.050|  -0.672|              0.502| -2.769|  1.359|
#|SexGirl                   |   -0.675|      0.692|  -0.976|              0.330| -2.036|  0.686|

##Linear mixed model on longitudinal data
library(phyloseq.extended)
ps.rarefied

##alpha diversity all
base_end_alpha <- microbiome::alpha(ps.rarefied, index = "all")
rownames(base_end_alpha)

##Phylogenetic diversity
base_end_faithpd <- phylodiv(ps.rarefied)
rownames(base_end_faithpd)

##Merge base_end_alpha and base_end_faithpd
base_end_alpadiversity_all <- merge(base_end_alpha, base_end_faithpd, by = 'row.names', all= TRUE)
View(base_end_alpadiversity_all)

base_end_alpadiversity_all$Intervention_groups <- factor(
  base_end_alpadiversity_all$Intervention_groups,
  levels = c("Placebo", "UR Improved", "UR Original")
)
base_end_alpadiversity_all$Baseline_Endline <- factor(base_end_alpadiversity_all$Baseline_Endline)


##Phylogenetic diversity
faithpd_intervention <- lm(pd ~ Baseline_Endline + Intervention_groups + Sex + Age_groups+ (1|SubjectID), base_end_alpadiversity_all)

res_faithintervention_faithpd = coefficients(summary(faithpd_intervention))


knitr::kable(res_intervention_faithpd, digits = 3)
#|                                                               | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|
#|(Intercept)                                                    |   22.747|      0.538|  42.318|              0.000|
#|factor(Baseline_Endline)Endline                                |   -5.969|      0.760|  -7.852|              0.000|
#|Intervention_groupsUR Improved                                 |   -0.448|      0.769|  -0.583|              0.560|
#|Intervention_groupsUR Original                                 |    0.784|      0.771|   1.017|              0.310|
#|factor(Baseline_Endline)Endline:Intervention_groupsUR Improved |    0.925|      1.089|   0.850|              0.396|
#|factor(Baseline_Endline)Endline:Intervention_groupsUR Original |    2.863|      1.092|   2.621|              0.009|

##Shannon diversity
shannon_intervention <- lm(diversity_shannon ~ factor(Baseline_Endline)*Intervention_groups + (1|SubjectID), base_end_alpadiversity_all)

#Summary of the model
res_intervention_shannon = coefficients(summary(shannon_intervention))
knitr::kable(res_intervention_shannon, digits = 3)

#|                                                               | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|
#|(Intercept)                                                    |    3.082|      0.058|  52.819|              0.000|
#|factor(Baseline_Endline)Endline                                |    0.392|      0.083|   4.755|              0.000|
#|Intervention_groupsUR Improved                                 |   -0.168|      0.084|  -2.006|              0.045|
#|Intervention_groupsUR Original                                 |    0.293|      0.084|   3.506|              0.000|
#|factor(Baseline_Endline)Endline:Intervention_groupsUR Improved |    0.226|      0.118|   1.914|              0.056|
#|factor(Baseline_Endline)Endline:Intervention_groupsUR Original |   -0.124|      0.119|  -1.047|              0.296|

##Pielou evenness
evenness_intervention <- lm(evenness_pielou ~ factor(Baseline_Endline)*Intervention_groups + (1|SubjectID), base_end_alpadiversity_all)

#Summary of the model
res_intervention_evenness = coefficients(summary(evenness_intervention))
knitr::kable(res_intervention_evenness, digits = 3)

#|                                                               | Estimate| Std. Error| t value| Pr(>&#124;t&#124;)|
#|(Intercept)                                                    |    0.628|      0.010|  62.864|              0.000|
#|factor(Baseline_Endline)Endline                                |    0.118|      0.014|   8.366|              0.000|
#|Intervention_groupsUR Improved                                 |   -0.027|      0.014|  -1.880|              0.061|
#|Intervention_groupsUR Original                                 |    0.049|      0.014|   3.406|              0.001|
#|factor(Baseline_Endline)Endline:Intervention_groupsUR Improved |    0.032|      0.020|   1.570|              0.117|
#|factor(Baseline_Endline)Endline:Intervention_groupsUR Original |   -0.041|      0.020|  -2.034|              0.042|

##Beta diversity measures
##Baseline
psrarefied_base

##Factorize
baseline_alphdiversity_all$Age_groups <- factor(baseline_alphdiversity_all$Age_groups)
baseline_alphdiversity_all$Anemia_class <- factor(baseline_alphdiversity_all$Anemia_class)
baseline_alphdiversity_all$Vitamin_A_deficiency <- factor(baseline_alphdiversity_all$Vitamin_A_deficiency)
baseline_alphdiversity_all$Iron_deficiency_anemia <- factor(baseline_alphdiversity_all$Iron_deficiency_anemia)
baseline_alphdiversity_all$Iron_deficiency <- factor(baseline_alphdiversity_all$Iron_deficiency)
baseline_alphdiversity_all$Systemic_inflammation <- factor(baseline_alphdiversity_all$Systemic_inflammation)
baseline_alphdiversity_all$Zinc_deficiency <- factor(baseline_alphdiversity_all$Zinc_deficiency)
baseline_alphdiversity_all$Hemoglobinopathy <- factor(baseline_alphdiversity_all$Hemoglobinopathy)
baseline_alphdiversity_all$Parasite_infection <- factor(baseline_alphdiversity_all$Parasite_infection)
baseline_alphdiversity_all$Stunting <- factor(baseline_alphdiversity_all$Stunting)
baseline_alphdiversity_all$Underweight <- factor(baseline_alphdiversity_all$Underweight)
baseline_alphdiversity_all$Gastrointestinal_inflammation <- factor(baseline_alphdiversity_all$Gastrointestinal_inflammation)

psrarefied_base_jsd <- phyloseq::distance(psrarefied_base, method = "jsd")
rarefied_base_df <- data.frame(sample_data(psrarefied_base))

md1_p0_bigfinal_jsd_genus_subject <- adonis2(
  psrarefied_base_jsd ~ Baseline_Endline * Intervention_groups, strata = rarefied_base_df$SubjectID,
  rarefied_base_df)


##Effect of intervention
ps.rarefied

##read the csv file for the complete metadata
library(readr)
responder_nonresponder_rarefied <- read_csv("responder_nonresponder_rarefied.csv")
View(responder_nonresponder_rarefied)

responder_nonresponder_rarefied1 <- sample_data(responder_nonresponder_rarefied)
sample_names(responder_nonresponder_rarefied1) <- responder_nonresponder_rarefied1$sample

##New phyloseq object with responder and non-responder

psrare_response <- phyloseq(otu_table(ps.rarefied), tax_table(ps.rarefied), phy_tree(ps.rarefied),
                            refseq(ps.rarefied), sample_data(responder_nonresponder_rarefied1))

psrare_response1 <- subset_samples(psrare_response, select_endline == "1")

sum(taxa_sums(psrare_response1) == 0)  ##3788

psrare_response1 <- phyloseq::prune_taxa(taxa_sums(psrare_response1) > 0, psrare_response1)  ##all of them are above zero

##Baseline response 
psrare_respon_base <- subset_samples(psrare_response, select_baseline == "1")

sum(taxa_sums(psrare_respon_base) == 0)  ##6105

psrare_respon_base <- phyloseq::prune_taxa(taxa_sums(psrare_respon_base) > 0, psrare_respon_base)  ##all of them are above zero

##Alpha diversity basline 
##calculate alpha diversity
psrare_resbase_alpha <- phyloseq::estimate_richness(psrare_respon_base, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson"))
psrare_resbase_alpha1 <- rownames_to_column(psrare_resbase_alpha, "sample")

##Phylogenetic diversity
psrare_resbase_phylo <- phylodiv(psrare_respon_base)

##merge faith pd and other diversity metrics
response_alphadiversity <- merge(psrare_resbase_alpha1, psrare_resbase_phylo, by = "sample")

##Comparsion for alpha diversity
response_base <- list(c("Responder", "Non-responder"))
p_value_cut = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

###Boxplot Shannon
ggplot(data = response_alphadiversity, aes(x = Response_1SD, y = Shannon, fill = Response_1SD)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8) +
  labs(x = "", y = "Shannon Diversity\n") +
  geom_jitter(width = 0.3) +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 11),
        strip.placement = "outside", strip.background = element_blank(), strip.text.x = element_text(size = 11)) +
stat_compare_means(comparisons = response_base, label = "p.signif", symnum.args = p_value_cut)

###Boxplot faith phylogenetic tree
ggplot(data = response_alphadiversity, aes(x = Response_1SD, y = pd, fill = Response_1SD)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8) +
  labs(x = "", y = "Shannon Diversity\n") +
  geom_jitter(width = 0.3) +
  theme_classic() +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 11),
        strip.placement = "outside", strip.background = element_blank(), strip.text.x = element_text(size = 11)) +
  stat_compare_means(comparisons = response_base, label = "p.signif", symnum.args = p_value_cut)

######################Plot beta diversity - unweighted unifrac
base_respo_dist = phyloseq::distance(psrare_respon_base, method="bray")
base_respo_ordination = ordinate(psrare_respon_base, method="PCoA", distance=base_respo_dist)

##Extract metadata from the rarefied phyloseq object
ps.rarefied_df <- data.frame(sample_data(ps.rarefied))
write.csv(ps.rarefied_df, "ps.rarefied_df.csv")

# Define ggplot2 custom colors:
allGroupColors <- c("brown3", "blue3")

plot_ordination(psrare_respon_base, base_respo_ordination, color="Response_1SD") + theme(aspect.ratio=1) +
  theme_classic() +
  geom_point(size = 2) +
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13), axis.line = element_line(size = 1, linetype = "solid"),
        axis.title.x = element_text(size=13.5, face="bold"),
        axis.title.y = element_text(size=13.5, face="bold"),
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12)) +
  labs(x = "PCoA 1 (14.2%)", y = "PCoA 2 (5.7%)") +
  stat_ellipse(aes(group = Response_1SD)) +
  scale_color_manual(values = allGroupColors) +
  guides(color = guide_legend(title = "   "))

ggsave('bray_curtis-final-sepetmber.tiff', dpi=300)

df_rare_respon_base <- data.frame(sample_data(psrare_respon_base))
base_respo_dist

#Run PERMANOVA on distances.
response_bray <- adonis2(base_respo_dist ~ Response_1SD,
                              data = df_rare_respon_base, permutations = 1000)




##Jaccard
ps.rarefied_jsd <- phyloseq::distance(ps.rarefied, method = "jsd")
rarefied_all_df <- data.frame(sample_data(ps.rarefied))

##Effect of intervention on jaccard distance
intervention_jaccard <- adonis2(
  ps.rarefied_jsd ~ Baseline_Endline * Intervention_groups, strata = rarefied_all_df$SubjectID,
  rarefied_all_df)

##Weighted unifrac
ps.rarefied_wunifrac <- phyloseq::distance(ps.rarefied, method = "wunifrac")

intervention_weiunifrac <- adonis2(
  ps.rarefied_wunifrac ~ Baseline_Endline * Intervention_groups, strata = rarefied_all_df$SubjectID,
  rarefied_all_df)

##Unweighted UniFrac
ps.rarefied_unifrac <- phyloseq::distance(ps.rarefied, method = "unifrac")

intervention_unifrac <- adonis2(
  ps.rarefied_unifrac ~ Baseline_Endline * Intervention_groups, strata = rarefied_all_df$SubjectID,
  rarefied_all_df)

##Bray curtis
ps.rarefied_bray <- phyloseq::distance(ps.rarefied, method = "bray")

intervention_bray <- adonis2(
  ps.rarefied_bray ~ Baseline_Endline * Intervention_groups, strata = rarefied_all_df$SubjectID,
  rarefied_all_df)


##Jaccard



##Answer to reviewers summary
prevabun_tree0_rooted

microbiome::summarize_phyloseq(prevabun_tree0_rooted)
##Min number of reads = 2948
##Maximum number of reads = 314227
##Average number of reads = 39355.5657894737
##Median number of reads = 35378.5


