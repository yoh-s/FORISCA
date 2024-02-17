#! /bin/bash
set -e
set -u

#NA_removed_sample-metadata-added-separater-ID.txt ##NA assigned to missing data is removed from R metadata

##calculation of core diversity metrics
qiime diversity core-metrics-phylogenetic \
    --i-table baseline-prevabun_35k.qza \
    --i-phylogeny prevabun-tree-rooted.qza \
    --p-sampling-depth 3500 \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --output-dir baseline-new-0321-metrics-results\

##alpha diversity metrics calculation	
#shannon_vector
qiime diversity alpha-group-significance \
    --i-alpha-diversity baseline-new-0321-metrics-results/shannon_vector.qza \
    --m-metadata-file  NA_removed_sample-metadata-added-separater-ID.txt \
    --o-visualization baseline-new-0321-metrics-results/shannon-group-significance.qzv

#faith_pd_vector	
qiime diversity alpha-group-significance \
  --i-alpha-diversity baseline-new-0321-metrics-results/faith_pd_vector.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --o-visualization baseline-new-0321-metrics-results/faith-pd-group-significance.qzv

#evenness_vector
qiime diversity alpha-group-significance \
  --i-alpha-diversity baseline-new-0321-metrics-results/evenness_vector.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --o-visualization baseline-new-0321-metrics-results/evenness-group-significance.qzv	

##Baseline unweighted unifrac distance matrix

#Sex
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Sex \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Sex-significance.qzv \
  --p-pairwise

#Intervention_groups
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Intervention_groups \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Intervention_groups-significance.qzv \
  --p-pairwise

#Age groups
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Age_groups \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Age_groups-significance.qzv \
  --p-pairwise

#Stunting
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Stunting \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Stunting-significance.qzv \
  --p-pairwise
 
#Stunting_level 
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Stunting_level \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Stunting_level-significance.qzv \
  --p-pairwise

#Underweight  
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Underweight \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Underweight-significance.qzv \
  --p-pairwise  
 
#Anemia_class
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Anemia_class \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Anemia_class-significance.qzv \
  --p-pairwise

#Systemic_inflammation 
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Systemic_inflammation \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Systemic_inflammation-significance.qzv \
  --p-pairwise
  
#Iron_deficiency
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Iron_deficiency \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Iron_deficiency-significance.qzv \
  --p-pairwise
  
#Iron_deficiency_anemia
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Iron_deficiency_anemia \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Iron_deficiency_anemia-significance.qzv \
  --p-pairwise
  
#Hemoglobinopathy  
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Hemoglobinopathy \
  --o-visualization baseline-new-0321-metrics-results/unweighted-unifrac-Hemoglobinopathy-significance.qzv \
  --p-pairwise
  
#Vitamin_A_deficiency
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Vitamin_A_deficiency \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/uw_unifrac-Vitamin_A_deficiency_Total-significance.qzv

#Zinc_deficiency
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Zinc_deficiency \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/uw_unifrac-Zinc_deficiency_Total-significance.qzv

#Parasite_infection
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Parasite_infection \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/uw_unifrac-Parasite_infection_Total-significance.qzv

#Gastrointestinal_inflammation
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Gastrointestinal_inflammation \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/uw_unifrac-Gastrointestinal_inflammation_Total-significance.qzv 

#Baseline Weighted unifrac distance

#Sex
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Sex \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Sex-significance.qzv \
  --p-pairwise

#Intervention_groups
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Intervention_groups \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Intervention_groups-significance.qzv \
  --p-pairwise

#Age_groups
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Age_groups \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Age_groups-significance.qzv \
  --p-pairwise

#Stunting
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Stunting \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Stunting-group-significance.qzv \
  --p-pairwise
 
#Stunting_level 
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Stunting_level \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Stunting_level-significance.qzv \
  --p-pairwise

#Underweight
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Underweight \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Underweight-significance.qzv \
  --p-pairwise

#Anemia_class
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Anemia_class \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Anemia_class-significance.qzv \
  --p-pairwise

#Systemic_inflammation
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Systemic_inflammation \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Systemic_inflammation-significance.qzv \
  --p-pairwise
  
#Iron_deficiency
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Iron_deficiency \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Iron_deficiency-significance.qzv \
  --p-pairwise
  
#Iron_deficiency_anemia
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Iron_deficiency_anemia \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Iron_deficiency_anemia-significance.qzv \
  --p-pairwise
  
#Hemoglobinopathy  
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Hemoglobinopathy \
  --o-visualization baseline-new-0321-metrics-results/weighted-unifrac-Hemoglobinopathy-significance.qzv \
  --p-pairwise
  
#Vitamin_A_deficiency
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Vitamin_A_deficiency \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/weighted_unifrac-Vitamin_A_deficiency_Total-significance.qzv

#Zinc_deficiency
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Zinc_deficiency \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/weighted_unifrac-Zinc_deficiency_Total-significance.qzv

#Parasite_infection
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Parasite_infection \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/weighted_unifrac-Parasite_infection_Total-significance.qzv

#Gastrointestinal_inflammation
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Gastrointestinal_inflammation \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/weighted_unifrac-Gastrointestinal_inflammation_Total-significance.qzv 

##Baseline bray-curtis distance

#Sex
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Sex \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix.Sex-significance.qzv \
  --p-pairwise

#Intervention_groups
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Intervention_groups \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix.Intervention_groups-significance.qzv \
  --p-pairwise

#Age_groups
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Age_groups \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix.Age_groups-significance.qzv \
  --p-pairwise

#Stunting
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Stunting \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-Stunting-group-significance.qzv \
  --p-pairwise
 
#Stunting_level 
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Stunting_level \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-unifrac-Stunting_level-significance.qzv \
  --p-pairwise

#Underweight
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Underweight \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-unifrac-Underweight-significance.qzv \
  --p-pairwise
 
#Anemia_class
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Anemia_class \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-Anemia_class-significance.qzv \
  --p-pairwise

#Systemic_inflammation 
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Systemic_inflammation \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-Systemic_inflammation-significance.qzv \
  --p-pairwise
  
#Iron_deficiency
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Iron_deficiency \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-Iron_deficiency-significance.qzv \
  --p-pairwise
  
#Iron_deficiency_anemia
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Iron_deficiency_anemia \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-Iron_deficiency_anemia-significance.qzv \
  --p-pairwise
  
#Hemoglobinopathy  
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Hemoglobinopathy \
  --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-Hemoglobinopathy-significance.qzv \
  --p-pairwise

#Vitamin_A_deficiency
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Vitamin_A_deficiency \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-Vitamin_A_deficiency_Total-significance.qzv

#Zinc_deficiency
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Zinc_deficiency \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-Zinc_deficiency-significance.qzv

#Parasite_infection
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Parasite_infection \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-Parasite_infection_Total-significance.qzv

#Gastrointestinal_inflammation
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/bray_curtis_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Gastrointestinal_inflammation \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/bray_curtis_distance_matrix-Gastrointestinal_inflammation_Total-significance.qzv 

##Baseline jaccard distance

#Sex
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Sex \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Sex-significance.qzv \
  --p-pairwise

#Intervention_groups
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Intervention_groups \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Intervention_groups-significance.qzv \
  --p-pairwise

#Age_groups
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Age_groups \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Age_groups-significance.qzv \
  --p-pairwise  

#Stunting
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Stunting \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Stunting-group-significance.qzv \
  --p-pairwise

#Stunting_level 
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Stunting_level \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-unifrac-Stunting_level-significance.qzv \
  --p-pairwise

#Underweight
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Underweight \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-unifrac-Underweight-significance.qzv \
  --p-pairwise

#Anemia_class
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Anemia_class \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Anemia_class-significance.qzv \
  --p-pairwise

#Systemic_inflammation
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Systemic_inflammation \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Systemic_inflammation-significance.qzv \
  --p-pairwise
  
#Iron_deficiency
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Iron_deficiency \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Iron_deficiency-significance.qzv \
  --p-pairwise
  
#Iron_deficiency_anemia
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Iron_deficiency_anemia \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Iron_deficiency_anemia-significance.qzv \
  --p-pairwise
  
#Hemoglobinopathy  
qiime diversity beta-group-significance \
  --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
  --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
  --m-metadata-column Hemoglobinopathy \
  --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Hemoglobinopathy-significance.qzv \
  --p-pairwise
  
#Vitamin_A_deficiency
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Vitamin_A_deficiency \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Vitamin_A_deficiency_Total-significance.qzv

#Zinc_deficiency
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Zinc_deficiency \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Zinc_deficiency-significance.qzv

#Parasite_infection
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Parasite_infection \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Parasite_infection_Total-significance.qzv

#Gastrointestinal_inflammation
qiime diversity beta-group-significance \
    --i-distance-matrix baseline-new-0321-metrics-results/jaccard_distance_matrix.qza \
    --m-metadata-file NA_removed_sample-metadata-added-separater-ID.txt \
    --m-metadata-column Gastrointestinal_inflammation \
    --p-pairwise \
    --o-visualization baseline-new-0321-metrics-results/jaccard_distance_matrix-Gastrointestinal_inflammation_Total-significance.qzv 