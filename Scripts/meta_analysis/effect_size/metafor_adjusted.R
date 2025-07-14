# setwd("~/switchdrive/humsub/Scripts/meta_analysis/effect_size/") # Please set the path correctly

# Load libraries
library(metafor)
library(tidyverse)
library(progress)

# Read data
relab <- readr::read_tsv('relab_species.tsv')
metadata <- readr::read_tsv('metadata_age_BMI.tsv')

# Long format
relab_long <- relab %>% 
  tidyr::pivot_longer(-sample_id, names_to = "species", values_to = "abundance")

# Merge with metadata
df <- merge(relab_long, metadata, by = "sample_id")

# Transform abundance
df$abundance <- asin(sqrt(df$abundance))

# Initialize result containers
df_results <- data.frame()
df_study_results_all <- data.frame()

# Progress bar
pb <- progress::progress_bar$new(total = length(unique(df$species)))

# Loop through species
for (species in unique(df$species)) {
  df_species <- df[df$species == species, ]
  
  df_study_results <- data.frame()
  
  for (study in unique(df_species$study)) {
    df_study <- df_species[df_species$study == study, ]
    
    # Skip if both groups not present
    if (!all(c("CRC", "control") %in% df_study$study_condition)) next
    
    # Fit covariate-adjusted linear model
    lm_model <- lm(abundance ~ study_condition + age + BMI, data = df_study)
    
    # Get adjusted means for CRC and control
    adjusted_means <- predict(lm_model, newdata = data.frame(
      study_condition = c("CRC", "control"),
      age = mean(df_study$age, na.rm = TRUE),
      BMI = mean(df_study$BMI, na.rm = TRUE)
    ))
    
    # Calculate effect size from adjusted means
    effect_size <- escalc(measure = "SMD", 
                          m1i = adjusted_means[1], 
                          m2i = adjusted_means[2], 
                          sd1i = sd(df_study[df_study$study_condition == 'CRC', 'abundance']), 
                          sd2i = sd(df_study[df_study$study_condition == 'control', 'abundance']), 
                          n1i = sum(df_study$study_condition == 'CRC'), 
                          n2i = sum(df_study$study_condition == 'control'))
    
    df_study_results <- rbind(df_study_results, 
                              data.frame(species = species, study = study, 
                                         effect_size = effect_size$yi, 
                                         variance = effect_size$vi))
  }
  
  df_study_results_all <- rbind(df_study_results_all, df_study_results)
  
  # Fit random-effects meta-analysis if there is sufficient data
  rma_success <- FALSE
  result <- tryCatch({
    re_model <- rma(yi = df_study_results$effect_size, 
                    vi = df_study_results$variance, method = "DL")
    rma_success <<- TRUE
    data.frame(species = species,
               mean_effect_size = re_model$beta[1],
               ci_lower = re_model$ci.lb[1],
               ci_upper = re_model$ci.ub[1],
               p_value = re_model$pval[1])
  }, error = function(e) {
    message(paste("Error for species:", species))
    return(data.frame(species = species, mean_effect_size = NA,
                      ci_lower = NA, ci_upper = NA, p_value = NA))
  })
  
  df_results <- rbind(df_results, result)
  pb$tick()
}

# Adjust p-values
df_results$p_value_adj <- p.adjust(df_results$p_value, method = "BH")

# Save results
write.table(df_results, "new_results/random_effects_sp.tsv", sep = "\t", row.names = FALSE)
write.table(df_study_results_all, "per_study_sp.tsv", sep = "\t", row.names = FALSE)
