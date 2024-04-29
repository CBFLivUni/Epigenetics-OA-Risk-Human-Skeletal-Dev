
##### PREPARATION #####

# Libraries
library(foreach)
library(doParallel)        

# Read in arguments
args <- commandArgs(trailingOnly=TRUE)

all_ncuts <- args[[1]]
n_clusters <- args[[2]]

# Define number of cores to use
cl <- parallel::makeCluster(n_clusters)

##### PREPARATION #####











##### READ IN DATA #####

# Read in data
for_step.data <- read.csv(paste0(res_dir,"stepfunction_regr.data.tsv"), sep="\t", header=TRUE, row.names=FALSE)

##### READ IN DATA #####












##### MAIN CODE #####

# Define CPGs of interest
cpgs_of_interest <- colnames(for_step.data)[grepl("^cg",colnames(colnames(for_step.data)))]

# Define right hand formula (of predictors)
predictors <- c("0","cut(Z_MAT_BMI, n_cut)", "Batch", "Z_STAGE_WEEKS", "sex", "Z_MAT_AGE")
predictor_formula <- paste(predictors, collapse=" + ")

# Initialise empty data frames for model parameters and metrics
all.modelparams_df <- data.frame()
all.modelmetric_df <- data.frame()

# Loop over number of cuts
for (n_cut  in all_ncuts) {
  
  # Report cut
  print(paste0("### Fitting LM with NCuts = ",n_cut))
  
  # Initialise empty data frames for model parameters and metrics
  cut.modelparams_df <- data.frame()
  cut.modelmetric_df <- data.frame()
  
  # Initialise number of cores
  doParallel::registerDoParallel(cl)
  
  # Loop over CpGs
  foreach(i=1:length(cpgs_of_interest)) %dopar% {
    
    # Get cpg name
    cpg <- cpgs_of_interest[[i]]
    
    # Fit step function
    fit_step <- lm(paste0(cpg, " ~ ", predictor_formula), data=for_step.data)
    
    # Get model output
    model_df <- tidy(fit_step)
    param_df <- glance(fit_step)
    
    # Assemble to larger dataframes
    cut.modelmetric_df <- rbind(CPG=cpg, NCUT=n_cut, all.modelmetric_df, model_df)
    cut.modelparams_df <- rbind(CPG=cpg, NCUT=n_cut, all.modelparams_df, model_df)
    
  }
  
  # Stop
  stopCluster(n_clusters)
  
  # Save
  saveRDS(all.modelmetric_df, paste0(res_dir,"step_function_regession.model_metric_df.rds"))
  saveRDS(all.modelparams_df, paste0(res_dir,"step_function_regession.model_params_df.rds"))
  
  # Add to larger dataframes
  all.modelmetric_df <- rbind(all.modelmetric_df, cut.modelmetric_df)
  all.modelparams_df <- rbind(all.modelparams_df, cut.modelparams_df)
  
}

##### MAIN CODE #####






