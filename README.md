# FORISCA
- Faecal microbiota of schoolchildren is associated with nutritional status, micronutrient 2 status and markers of inflammation
## **Scripts**
- _**baseline_top_taxa.R**_ is the R-script used to plot the top taxa at the baseline.  
- _**bray_curtis_PCoA_plot.R**_ is the R-script used to generate PCoA plot of the bray_curtis at the Baseline and Endline.  
- _**core_microbiome_final.R**_ is the R-script used to calculate the core microbiota at the baseline  
- _**forward_dada2.R**_ is the R-script used to process the raw sequence.  
- _**__config.R**_ is a configuration file that must be located in the same working directory as the _forward_dada2.R_ script.
- _**glm_lefse.R**_ implements a linear mixed model on microbial features identified as significant by Linear discriminant analysis Effect Size (LEfSe).
- _**linear_model_alphadiversity.R**_ is a linear mixed model on alpha diversity metrics at the baseline.
- _**Maaslin2.R**_ is the R-script used to generate the linear mixed models.  
- _**phyloseq_final.R**_ is the R-script used to generate the phyloseq object and prilimary processing.  
- _**QIIME2_alpha_beta_diversity_baseline**_ folder contains command used to generate alpha and beta diversity differences at the baseline, metadata, feature-table, tree and README file.  
- _**QIIME2_longitudinal_alpha_diversity**_ folder contains command used to generate longitudinal differences in alpha diversity, metadata, feature-table, tree and README file.  
- _**QIIME2_longitudinal_beta_diversity**_ folder contains command used to generate longitudinal differences in beta diversity, metadata, feature-table, tree and README file.  
- _**test_data**_ folder contains test data (13-sample (.fastq file)and metadata (.xlsx)).
## **Phyloseq object**
- _**prevabun_tree0_rooted.rds**_ is the final phyloseq object generated using phyloseq_final.R script  
- please modify the _**"dbs"**_ value inside the _"__config.R"_ file as follows to be able to run the _forward_dad2.R_ file to run properly
    - we assume you downloaded the [silva_nr99_v138.1_wSpecies_train_set.fa.gz] (https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1) file and saved it in a folder called "database" in your R working directory (your_R_working_directory) in a folder named "database"
    - dbs = list("16S" = "~/your_R_working_directory/database/silva_nr99_v138.1_wSpecies_train_set.fa.gz")
