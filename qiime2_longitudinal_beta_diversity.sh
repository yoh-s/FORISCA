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
##Beta diversity significance test  

##Weighted_unifrac_distance_matrix

#Sex
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Sex \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Sex.qzv

#Intervention groups
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Intervention_groups \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Intervention_groups.qzv

#Age_group
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Age_group \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Age_group.qzv

#Stunting
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Stunting \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Stunting.qzv

#Stunting_level
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Stunting_level \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Stunting_level.qzv

#Underweight
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Underweight \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Underweight.qzv  

#Anemia_class
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Anemia_class \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-differences-Anemia_class.qzv 

#Systemic_inflammation
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Systemic_inflammation \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-differences-Systemic_inflammation.qzv  
  
#Iron_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Iron_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Iron_deficiency.qzv

#Iron_deficiency_anemia
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Iron_deficiency_anemia \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distance-weighted_unifrac-Iron_deficiency_anemia.qzv

#Hemoglobinopathy
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Hemoglobinopathy \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Hemoglobinopathy.qzv

#Vitamin_A_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Vitamin_A_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Vitamin_A_deficiency.qzv

#Zinc_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Zinc_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Zinc_deficiency.qzv

#Parasite_infection
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/weighted_unifrac_distance_matrix.qza \
  --p-group-column Parasite_infection \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-weighted_unifrac-Parasite_infection.qzv

##UNweighted_unifrac_distance_matrix

#Sex
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Sex \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Sex.qzv

#Intervention groups
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Intervention_groups \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Intervention_groups.qzv

#Age_group
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Age_group \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Age_group.qzv

#Stunting
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Stunting \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Stunting.qzv

#Stunting_level
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Stunting_level \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Stunting_level.qzv

#Underweight
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Underweight \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Underweight.qzv  

#Anemia_class
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Anemia_class \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-differences-Anemia_class.qzv 

#Systemic_inflammation
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Systemic_inflammation \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-differences-Systemic_inflammation.qzv  
  
#Iron_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Iron_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Iron_deficiency.qzv

#Iron_deficiency_anemia
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Iron_deficiency_anemia \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distance-UNweighted_unifrac-Iron_deficiency_anemia.qzv

#Hemoglobinopathy
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Hemoglobinopathy \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Hemoglobinopathy.qzv

#Vitamin_A_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Vitamin_A_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Vitamin_A_deficiency.qzv

#Zinc_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Zinc_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Zinc_deficiency.qzv

#Parasite_infection
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --p-group-column Parasite_infection \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-UNweighted_unifrac-Parasite_infection.qzv

##bray_curtis_distance_matrix 

#Sex
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Sex \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Sex.qzv

#Intervention groups
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Intervention_groups \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Intervention_groups.qzv

#Age_group
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Age_group \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Age_group.qzv

#Stunting
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Stunting \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Stunting.qzv

#Stunting_level
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Stunting_level \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Stunting_level.qzv

#Underweight
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Underweight \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Underweight.qzv  

#Anemia_class
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Anemia_class \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Anemia_class.qzv 

#Systemic_inflammation
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Systemic_inflammation \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Systemic_inflammation.qzv  
  
#Iron_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Iron_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Iron_deficiency.qzv

#Iron_deficiency_anemia
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Iron_deficiency_anemia \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distance-bray_curtis-Iron_deficiency_anemia.qzv

#Hemoglobinopathy
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Hemoglobinopathy \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Hemoglobinopathy.qzv

#Vitamin_A_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Vitamin_A_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Vitamin_A_deficiency.qzv

#Zinc_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Zinc_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Zinc_deficiency.qzv

#Parasite_infection
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/bray_curtis_distance_matrix.qza \
  --p-group-column Parasite_infection \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-bray_curtis-Parasite_infection.qzv

##jaccard_distance_matrix.qza

#Sex
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Sex \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Sex.qzv

#Intervention groups
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Intervention_groups \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Intervention_groups.qzv

#Age_group
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Age_group \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Age_group.qzv

#Stunting
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Stunting \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Stunting.qzv

#Stunting_level
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Stunting_level \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Stunting_level.qzv

#Underweight
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Underweight \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Underweight.qzv  

#Anemia_class
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Anemia_class \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Anemia_class.qzv 

#Systemic_inflammation
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Systemic_inflammation \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Systemic_inflammation.qzv  
  
#Iron_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Iron_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Iron_deficiency.qzv

#Iron_deficiency_anemia
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Iron_deficiency_anemia \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distance-jaccard-Iron_deficiency_anemia.qzv

#Hemoglobinopathy
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Hemoglobinopathy \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Hemoglobinopathy.qzv

#Vitamin_A_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Vitamin_A_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Vitamin_A_deficiency.qzv

#Zinc_deficiency
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Zinc_deficiency \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Zinc_deficiency.qzv

#Parasite_infection
qiime longitudinal pairwise-distances \
  --m-metadata-file metadata-modified-qiime2-paired.txt \
  --i-distance-matrix BaseEndline-metrics-results/jaccard_distance_matrix.qza \
  --p-group-column Parasite_infection \
  --p-state-column Baseline_Endline \
  --p-state-1 Baseline \
  --p-state-2 Endline \
  --p-individual-id-column SubjectID \
  --p-replicate-handling random \
  --o-visualization pairwise-distances-jaccard-Parasite_infection.qzv