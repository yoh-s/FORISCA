#Load libraries
library(phyloseq)
library(phyloseq.extended)
library(microbiome)
library(stats)
library(dplyr)

# Assuming psrarefied_base is a preloaded Phyloseq object
# Calculate alpha diversity measures
baseline_alpha <- microbiome::alpha(psrarefied_base, index = "all")

# Calculate Phylogenetic diversity (Faith's PD)
baseline_faithpd <- phylodiv(psrarefied_base)

# Convert row names to a column for merging
baseline_alpha <- baseline_alpha %>% rownames_to_column(var = "sampleid2")

# Merge Faith PD with other alpha diversity metrics
baseline_alphdiversity_all <- right_join(baseline_faithpd, baseline_alpha, by = "sampleid2")

# Convert specified columns to factors for appropriate modeling
columns_to_factorize <- c("Age_groups", "Anemia_class", "Vitamin_A_deficiency",
                          "Iron_deficiency_anemia", "Iron_deficiency", "Systemic_inflammation",
                          "Zinc_deficiency", "Hemoglobinopathy", "Parasite_infection",
                          "Stunting", "Underweight", "Gastrointestinal_inflammation")

baseline_alphdiversity_all[columns_to_factorize] <- lapply(baseline_alphdiversity_all[columns_to_factorize], factor)

# List of health and demographic variables to analyze
variables_to_analyze <- c("Anemia_class", "Vitamin_A_deficiency", "Iron_deficiency_anemia", 
                          "Iron_deficiency", "Systemic_inflammation", "Zinc_deficiency", 
                          "Hemoglobinopathy", "Parasite_infection", "Stunting", 
                          "Underweight", "Gastrointestinal_inflammation")

#Perform linear modeling of Pielou evenness and summarize results
analyze_evenness <- function(dependent_variable, data_frame) {
  model <- lm(formula = paste("evenness_pielou ~", dependent_variable, "+ Sex + Age_groups"), data = data_frame)
  results <- cbind(coefficients(summary(model)), confint(model))
  return(results)
}

# Apply the analysis function to each variable and print the results
results_pielou <- lapply(variables_to_analyze, function(var) {
  results <- analyze_evenness(var, baseline_alphdiversity_all)
  cat("Results for:", var, "\n")
  print(knitr::kable(results, digits = 3))
})

#Perform linear modeling of Shannon diversity and summarize results
analyze_shannon <- function(dependent_variable, data_frame) {
  model <- lm(formula = paste("diversity_shannon ~", dependent_variable, "+ Sex + Age_groups"), data = data_frame)
  results <- cbind(coefficients(summary(model)), confint(model))
  return(results)
}

# Apply the analysis function to each variable and print shannon diversity results
results_shannon <- lapply(variables_to_analyze, function(var) {
  results <- analyze_evenness(var, baseline_alphdiversity_all)
  cat("Results for:", var, "\n")
  print(knitr::kable(results, digits = 3))
})

#######################
#Perform linear modeling of Faith's PD and summarize results
analyze_faithpd <- function(dependent_variable, data_frame) {
  model <- lm(formula = paste("pd ~", dependent_variable, "+ Sex + Age_groups"), data = data_frame)
  results <- cbind(coefficients(summary(model)), confint(model))
  return(results)
}

# Apply the analysis function to each variable and print Faith's PD diversity
results_faithpd <- lapply(variables_to_analyze, function(var) {
  results <- analyze_evenness(var, baseline_alphdiversity_all)
  cat("Results for:", var, "\n")
  print(knitr::kable(results, digits = 3))
})