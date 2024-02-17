#! /bin/bash
set -e
set -u

####Calculation of core diversity metrics
####subjects without endline or baseline information were excluded from longitudinal analysis
###Qiime longitudinal pairwise does not work without matching baseline and endline subjects
##metadata is modified
###feature table (i-table) is modified

qiime diversity core-metrics-phylogenetic \
  --i-table prevabun-paired-filtered-table.qza \
  --i-phylogeny prevabun-tree-rooted.qza \
  --p-sampling-depth 3500 \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --output-dir BaseEndline-metrics-results\

#####################################################
##Alpha diversity significance test

##Shannon vector
#Sex
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Sex \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Sex.qzv

#Intervention_groups
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Intervention_groups \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Intervention_groups.qzv

#Age_group
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Age_group \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Age_groups.qzv

#Stunting
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Stunting \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Stunting.qzv

#Stunting_level
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Stunting_level \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Stunting_level.qzv  

#Underweight
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Underweight \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Underweight.qzv    
 
#Anemia_class
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Anemia_class \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Anemia_class.qzv 

#Systemic_inflammation
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Systemic_inflammation \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Inflammation.qzv  

#Iron_deficiency
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Iron_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Iron_deficiency.qzv

#Iron_deficiency_anemia
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Iron_deficiency_anemia \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Iron_deficiency_anemia.qzv

#Hemoglobinopathy  
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Hemoglobinopathy \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Hemoglobinopathy.qzv

#Vitamin_A_deficiency
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Vitamin_A_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Vitamin_A_deficiency.qzv

#Zinc_deficiency
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Zinc_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Zinc_deficiency.qzv

#Parasite_infection
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/shannon_vector.qza \
  --p-metric shannon_entropy \
  --p-group-column Parasite_infection \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-shannon-Parasite_infection.qzv

#####################################################
##Faith_pd_vector

#Sex
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Sex \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Sex.qzv

#Intervention_groups
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Intervention_groups \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Intervention_groups.qzv
  
#Age_group
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Age_group \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Age_group.qzv

#Stunting
  
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Stunting \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Stunting.qzv

#Stunting_level

qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Stunting_level \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Stunting_level.qzv  

#Underweight
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Underweight \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Underweight.qzv 

#Anemia_class
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Anemia_class \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Anemia_class.qzv 

#Systemic_inflammation
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Systemic_inflammation \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Systemic_inflammation.qzv  

#Iron_deficiency
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Iron_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Iron_deficiency.qzv

#Iron_deficiency_anemia
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Iron_deficiency_anemia \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Iron_deficiency_anemia.qzv

#Hemoglobinopathy  
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Hemoglobinopathy \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Hemoglobinopathy.qzv

#Vitamin_A_deficiency
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Vitamin_A_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Vitamin_A_deficiency.qzv

#Zinc_deficiency
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Zinc_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Zinc_deficiency.qzv

#Parasite_infection
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/faith_pd_vector.qza \
  --p-metric faith_pd \
  --p-group-column Parasite_infection \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-faith_pd_vector-Parasite_infection.qzv

#####################################################
##pielou_evenness

#Sex
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Sex \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Sex.qzv

#Intervention_groups
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Intervention_groups \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Intervention_groups.qzv

#Age_group
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Age_group \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Age_group.qzv

#Stunting
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Stunting \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Stunting.qzv

#Stunting_level
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Stunting_level \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Stunting_level.qzv  

#Anemia_class
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Anemia_class \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Anemia_class.qzv 

#Systemic_inflammation
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Systemic_inflammation \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Systemic_inflammation.qzv  

#Iron_deficiency
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Iron_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Iron_deficiency.qzv

#Iron_deficiency_anemia
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Iron_deficiency_anemia \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Iron_deficiency_anemia.qzv

#Hemoglobinopathy  
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Hemoglobinopathy \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Hemoglobinopathy.qzv

#Vitamin_A_deficiency
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Vitamin_A_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Vitamin_A_deficiency.qzv

#Zinc_deficiency
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Zinc_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Zinc_deficiency.qzv

#Parasite_infection
qiime longitudinal pairwise-differences \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --m-metadata-file BaseEndline-metrics-results/evenness_vector.qza \
  --p-metric pielou_evenness \
  --p-group-column Parasite_infection \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-differences-pielou_evenness-Parasite_infection.qzv
