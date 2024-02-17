set.seed(2023)
#Load necessary packages
library(stats)
library(data.table)
library(phyloseq)
library(dplyr)
library(stringr)
library(Tjazi)
library(curatedMetagenomicData)
library(curatedMetagenomicAnalyses)

psBase_filter
taxa_are_rows(psBase_filter) ##FALSE

ps_baseline_nozero <- phyloseq::prune_taxa(taxa_sums(ps_baseline) > 0, ps_baseline)  ##3479
taxa_are_rows(ps_baseline_nozero) ##FALSE

#Biomarkers of anemia identified by lefse
biomarker_anemia_genus <- c("Prevotella_7", "Prevotellaceae NK3B31 group", "Limosilactobacillus", "Anaerostipes", 
  "Lachnospira", "[Eubacterium] hallii group", "[Ruminococcus] torques group","Faecalibacterium", "Klebsiella")

biomarker_anemia_family <- c("Erysipelatoclostridiaceae", "Lachnospiraceae")

##Genus level glom
psbase_zerorgenus <- tax_glom(ps_baseline_nozero, taxrank = "Genus")
psbase_zerorgenusclr <- microbiome::transform(psbase_zerorgenus, 'clr')

##Family level glom
psbase_zerorfamily <- tax_glom(ps_baseline_nozero, taxrank = "Family")
psbase_zerorfamilyclr <- microbiome::transform(psbase_zerorfamily, 'clr')

##filter the significant features from the family
anemia_family <- subset_taxa(psbase_zerorfamilyclr, Family %in% biomarker_anemia_family)
#filter the significant features from the genus
anemia_genus <- subset_taxa(psbase_zerorgenusclr, Genus %in% biomarker_anemia_genus)

##extract otu table
otu_genus_anemia_t <- t(data.frame(otu_table(anemia_genus)))
otu_family_anemia_t <- t(data.frame(otu_table(anemia_family)))
#merge genus and family
anemia_sig_all <- rbind(otu_genus_anemia_t, otu_family_anemia_t)

##extract sample_data
sample_data_glm <- data.frame(sample_data(ps_baseline))
sample_data_glm$Anemia_class <- factor(sample_data_glm$Anemia_class)
nrow(otu_genus_anemia_t)
#9
nrow(sample_data_glm)
#366

sample_data_glm$Anemia_class <- factor(sample_data_glm$Anemia_class)

# assuming that 'otu_table' and 'sample_data' are the extracted data
anemia_glm <- fw_glm(x = anemia_sig_all, f = ~ Anemia_class + Sex + Age_groups, metadata = sample_data_glm, adjust.method = "BH", format = "brief")
glimpse(anemia_glm)

#Plot the features that show a group effect at q < 0.1
anemia_BH <- anemia_sig_all[anemia_glm[anemia_glm$`Anemia_classNo anemia Pr(>|t|).BH` < 0.1,"feature"],]
View(glimpse(anemia_glm))

dir.create("glm_model")
df_anemia_glm <- data.frame(glimpse(anemia_glm))
write.csv(df_anemia_glm, "glm_model/df_anemia_glm.csv")
rownames(anemia_BH)

##"ASV4"   "ASV13"  "ASV82"  "ASV88"  "ASV109" "ASV23"  "ASV24" 

###################################################################################################################################
##Iron deficiency anemia
#Biomarkers of anemia identified by lefse
biomarker_ida_genus <- c("Prevotella_7", "Prevotellaceae NK3B31 group", "HT002", "Lachnospira", "Faecalibacterium", "Ruminococcus",
                            "Lachnospiraceae UCG-001", "[Eubacterium] hallii group","Lachnospiraceae UCG-010","Monoglobus")

biomarker_ida_family <- c("Lachnospiraceae", "Monoglobaceae", "Ruminococcaceae")

biomarker_ida_order <- c("Lachnospirales", "Monoglobales", "Oscillospirales")

##Genus level glom
psbase_zerorgenus 
psbase_zerorgenusclr  ##clr transformed

##Family level glom
psbase_zerorfamily
psbase_zerorfamilyclr

##Order level glom
psbase_zerorrder <- tax_glom(ps_baseline_nozero, taxrank = "Order")
psbase_zerorderclr <- microbiome::transform(psbase_zerorrder, 'clr')

##filter the significant features from the family
ida_family <- subset_taxa(psbase_zerorfamilyclr, Family %in% biomarker_ida_family)

#filter the significant features from the genus
ida_genus <- subset_taxa(psbase_zerorgenusclr, Genus %in% biomarker_ida_genus)

#filter the significant features from the genus
ida_order <- subset_taxa(psbase_zerorderclr, Order %in% biomarker_ida_order)

##extract otu table
otu_genus_ida_t <- t(data.frame(otu_table(ida_genus)))
otu_family_ida_t <- t(data.frame(otu_table(ida_family)))
otu_order_ida_t <- t(data.frame(otu_table(ida_order)))

#merge genus and family
ida_sig_all <- rbind(otu_genus_ida_t, otu_family_ida_t, otu_order_ida_t)

##extract sample_data
sample_data_glm <- data.frame(sample_data(ps_baseline))
sample_data_glm$Iron_deficiency_anemia <- factor(sample_data_glm$Iron_deficiency_anemia)


# assuming that 'otu_table' and 'sample_data' are the extracted data
ida_glm <- fw_glm(x = ida_sig_all, f = ~ Iron_deficiency_anemia + Sex + Age_groups, metadata = sample_data_glm, adjust.method = "BH", format = "brief")

##save as csv
df_ida_glm <- data.frame(glimpse(ida_glm))
write.csv(df_ida_glm, "glm_model/df_ida_glm.csv")


#Plot the features that show a group effect at q < 0.1
# Check the condition
condition_ida <- ida_glm$`Iron_deficiency_anemiaYes Pr(>|t|).BH` < 0.1
print(condition_ida)

# Subset ida_glm based on the condition
subset_glm_ida <- ida_glm[condition_ida, ]

######################################################################################################################
##Vitamin A deficiency
#Biomarkers of anemia identified by lefse
biomarker_vitaA_genus <- c("Gordonibacter", "Roseburia", "[Eubacterium] ruminantium group", "UBA1819", "Romboutsia",
                         "[Eubacterium] siraeum group")

biomarker_vitaA_family <- c("Desulfovibrionaceae", "Lactobacillaceae", "Lachnospiraceae")

##Genus level glom
psbase_zerorgenus 
psbase_zerorgenusclr  ##clr transformed

##Family level glom
psbase_zerorfamily
psbase_zerorfamilyclr


##filter the significant features from the family
vitA_family <- subset_taxa(psbase_zerorfamilyclr, Family %in% biomarker_vitaA_family)

#filter the significant features from the genus
vitA_genus <- subset_taxa(psbase_zerorgenusclr, Genus %in% biomarker_vitaA_genus)


##extract otu table
otu_genus_vitA_t <- t(data.frame(otu_table(vitA_genus)))
otu_family_vitA_t <- t(data.frame(otu_table(vitA_family)))


#merge genus and family
vitA_sig_all <- rbind(otu_genus_vitA_t, otu_family_vitA_t)

##extract sample_data
sample_data_glm <- data.frame(sample_data(ps_baseline))
sample_data_glm$Vitamin_A_deficiency <- factor(sample_data_glm$Vitamin_A_deficiency)


# assuming that 'otu_table' and 'sample_data' are the extracted data
vitA_glm <- fw_glm(x = vitA_sig_all, f = ~ Vitamin_A_deficiency + Sex + Age_groups, metadata = sample_data_glm, adjust.method = "BH", format = "brief")
View(glimpse(vitA_glm))

df_vitA_glm <- data.frame(glimpse(vitA_glm))
write.csv(df_vitA_glm, "glm_model/df_vitA_glm.csv")

#Plot the features that show a group effect at q < 0.1
# Check the condition
condition_vitA <- vitA_glm$`Vitamin_A_deficiencyYes Pr(>|t|).BH` < 0.1
print(condition_vitA)

# Subset vitA_glm based on the condition
subset_glm_vitA <- vitA_glm[condition_vitA, ]

#############################################################################################
##Stunting
#Biomarkers of anemia identified by lefse
biomarker_stunting_genus <- c("Raoultibacter", "Odoribacter", "Prevotella_7", "Alistipes", "Holdemanella", "Family XIII AD3011 group",
                              "Succinivibrio")

biomarker_stunting_family <- c("Succinivibrionaceae", "Erysipelotrichaceae")

##Genus level glom
psbase_zerorgenus 
psbase_zerorgenusclr  ##clr transformed

##Family level glom
psbase_zerorfamily
psbase_zerorfamilyclr


##filter the significant features from the family
stunting_family <- subset_taxa(psbase_zerorfamilyclr, Family %in% biomarker_stunting_family)

#filter the significant features from the genus
stunting_genus <- subset_taxa(psbase_zerorgenusclr, Genus %in% biomarker_stunting_genus)

##extract otu table
otu_genus_stunt_t <- t(data.frame(otu_table(stunting_genus)))
otu_family_stunt_t <- t(data.frame(otu_table(stunting_family)))
write.csv(otu_family_stunt_t, "otu_family_stunt_t.csv")

# Assuming your matrix is named 'your_matrix'
row_to_remove <- which(rownames(otu_family_stunt_t) == "ASV28")

# Remove the row with the specified index
otu_family_stunt_tf <- otu_family_stunt_t[-row_to_remove, , drop = FALSE]

#merge genus and family
stunting_sig_all <- rbind(otu_genus_stunt_t, otu_family_stunt_tf)

##extract sample_data
sample_data_glm <- data.frame(sample_data(ps_baseline))
sample_data_glm$Stunting <- factor(sample_data_glm$Stunting)

# assuming that 'otu_table' and 'sample_data' are the extracted data
stunting_glm <- fw_glm(x = vitA_sig_all, f = ~ Stunting + Sex + Age_groups, metadata = sample_data_glm, adjust.method = "BH", format = "brief")
View(glimpse(stunting_glm))

#Plot the features that show a group effect at q < 0.1
# Check the condition
condition_stunting <- stunting_glm$`StuntingYes Pr(>|t|).BH` < 0.1
print(condition_stunting)

# Subset stunting_glm based on the condition
subset_glm_stunting <- stunting_glm[condition_stunting, ]  ##No significant value

#############################################################################################
##Systemic inflammation 
#Biomarkers of anemia identified by lefse
biomarker_sysinfla_genus <- c("Butyricimonas", "Alloprevotella", "Holdemania", "Enterococcus", "Lapidilactobaillus", "Lentilactobacillus",
                              "CAG-352", "Romboutsia", "Enterobacter", "Klebsiella", "Succinivibrio")

biomarker_sysinfla_family <- c("Enterococcaceae", "Peptostreptococcaceae", "Succinivibrionaceae")

biomarker_sysinfla_order <- c("Gastranaerophilales", "RF39")

##Genus level glom
psbase_zerorgenus 
psbase_zerorgenusclr  ##clr transformed

##Family level glom
psbase_zerorfamily
psbase_zerorfamilyclr

##Order level glom
##Order level glom
psbase_zerorrder
psbase_zerorderclr 

##filter the significant features from the family
sysinfla_family <- subset_taxa(psbase_zerorfamilyclr, Family %in% biomarker_sysinfla_family)

#filter the significant features from the genus
sysinfla_genus <- subset_taxa(psbase_zerorgenusclr, Genus %in% biomarker_sysinfla_genus)

#filter the significant features from the order
sysinfla_order <- subset_taxa(psbase_zerorderclr, Order %in% biomarker_sysinfla_order)

##extract otu table
otu_genus_sysinfla_t <- t(data.frame(otu_table(sysinfla_genus)))
otu_family_sysinfla_t <- t(data.frame(otu_table(sysinfla_family)))
otu_order_sysinfla_t <- t(data.frame(otu_table(sysinfla_order)))

#merge genus and family
sysinfla_sig_all <- rbind(otu_genus_sysinfla_t, otu_family_sysinfla_t, otu_order_sysinfla_t)

##extract sample_data
sample_data_glm <- data.frame(sample_data(ps_baseline))
sample_data_glm$Systemic_inflammation <- factor(sample_data_glm$Systemic_inflammation)

sample_data_glm$Systemic_inflammation <- factor(sample_data_glm$Systemic_inflammation)
sample_data_glm$Sex <- factor(sample_data_glm$Sex)
sample_data_glm$Age_groups <- factor(sample_data_glm$Age_groups)


# assuming that 'otu_table' and 'sample_data' are the extracted data
sysinfla_glm <- fw_glm(x = sysinfla_sig_all, f = ~ Systemic_inflammation + Sex + Age_groups, metadata = sample_data_glm, adjust.method = "BH", format = "brief")
View(glimpse(sysinfla_glm))

df_sysinfla_glm <- data.frame(glimpse(sysinfla_glm))
write.csv(df_sysinfla_glm, "glm_model/df_sysinfla_glm.csv")
#Plot the features that show a group effect at q < 0.1
# Check the condition
condition_sysinfla <- sysinfla_glm$`Systemic_inflammationYes Pr(>|t|).BH` < 0.1
print(condition_sysinfla)

# Subset stunting_glm based on the condition
subset_glm_sysinfla <- sysinfla_glm[condition_sysinfla, ]  ##No significant value

#####################################################################################
##Parasitic infection
#Biomarkers of anemia identified by lefse
biomarker_parainfec_genus <- c("Gordonibacter", "Kurthia", "UCG-004", "Holdemanella", "Weissella", "Sarcina",
                              "Fastidiosipila")

biomarker_parainfec_family <- c("Desulfovibrionaceae", "Erysipelotrichaceae", "Hungateiclostridiaceae")

##Genus level glom
psbase_zerorgenus 
psbase_zerorgenusclr  ##clr transformed

##Family level glom
psbase_zerorfamily
psbase_zerorfamilyclr

##filter the significant features from the family
parainfec_family <- subset_taxa(psbase_zerorfamilyclr, Family %in% biomarker_parainfec_family)

#filter the significant features from the genus
parainfec_genus <- subset_taxa(psbase_zerorgenusclr, Genus %in% biomarker_parainfec_genus)

##extract otu table
otu_genus_parainfec_t <- t(data.frame(otu_table(parainfec_genus)))
otu_family_parainfec_t <- t(data.frame(otu_table(parainfec_family)))

#merge genus and family
parainfec_sig_all <- rbind(otu_genus_parainfec_t, otu_family_parainfec_t)

##extract sample_data
sample_data_glm <- data.frame(sample_data(ps_baseline))
sample_data_glm$Parasite_infection <- factor(sample_data_glm$Parasite_infection)

# assuming that 'otu_table' and 'sample_data' are the extracted data
parainfec_glm <- fw_glm(x = parainfec_sig_all, f = ~ Parasite_infection + Sex + Age_groups, metadata = sample_data_glm, adjust.method = "BH", format = "brief")
View(glimpse(parainfec_glm))

df_parainfec_glm <- data.frame(glimpse(parainfec_glm))
write.csv(df_parainfec_glm, "glm_model/df_parainfec_glm.csv")

#Plot the features that show a group effect at q < 0.1
# Check the condition
condition_parainfec <- parainfec_glm$`Parasite_infectionYes Pr(>|t|).BH` < 0.1
print(condition_parainfec)

# Subset stunting_glm based on the condition
subset_glm_parainfec <- parainfec_glm[condition_parainfec, ]  ##No significant value

######################################################################################
##Age

#Biomarkers of anemia identified by lefse
biomarker_age_genus <- c("Prevotella", "Prevotella_7", "Agrilactobacillus", "Lacticaseibacillus", "Lactobacillus", "Agathobacter",
                               "Lachnospira", "Lachnospiraceae NK4A136 group", "Marvinbryantia", "Roseburia", "[Eubacterium] eligens group",
                               "Butyricicoccus", "Colidextribacter", "Anaerococcus", "Dialister", "Cetobacterium")

biomarker_age_family <- c("Fusobacteriaceae", "Veillonellaceae", "Family XI", "Butyricicoccaceae")

##Genus level glom
psbase_zerorgenus 
psbase_zerorgenusclr  ##clr transformed

##Family level glom
psbase_zerorfamily
psbase_zerorfamilyclr

##filter the significant features from the family
age_family <- subset_taxa(psbase_zerorfamilyclr, Family %in% biomarker_age_family)

#filter the significant features from the genus
age_genus <- subset_taxa(psbase_zerorgenusclr, Genus %in% biomarker_age_genus)

##extract otu table
otu_genus_age_t <- t(data.frame(otu_table(age_genus)))
otu_family_age_t <- t(data.frame(otu_table(age_family)))

#merge genus and family
age_sig_all <- rbind(otu_genus_age_t, otu_family_age_t)

##extract sample_data
sample_data_glm <- data.frame(sample_data(ps_baseline))
sample_data_glm$Age_groups <- factor(sample_data_glm$Age_groups)

# assuming that 'otu_table' and 'sample_data' are the extracted data
age_glm <- fw_glm(x = age_sig_all, f = ~ Age_groups + Sex + Age_groups, metadata = sample_data_glm, adjust.method = "BH", format = "brief")
glimpse(age_glm)

#Plot the features that show a group effect at q < 0.1
# Check the condition
condition_age <- age_glm$`Age_groups6-9_years Pr(>|t|).BH` < 0.1
print(condition_age)

# Subset stunting_glm based on the condition
subset_glm_age <- age_glm[condition_age, ]  ##No significant value
