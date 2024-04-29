

GetExpData <- function(dmr, meta) {
  
  # Dependencies
  require(BBmisc)
  
  #get the stage of each sample
  expData <- dmr
  colnames(expData) <- as.character(meta$stage_weeks)
  cn <- colnames(expData)
  
  #take the mean M values for samples from the same stage
  expData<-data.frame(sapply(unique(cn), function(g) rowMeans(expData[,cn==g,drop=FALSE])))
  colnames(expData) <- gsub("^X", "",  colnames(expData))
  expData<-expData[order(as.numeric(names(expData)))]
  
  #standardise the data to focus on the shape of the CpG profiles across stage
  prot_norm_t<-as.data.frame(t(expData))
  sample_norm<-BBmisc::normalize(prot_norm_t, method="standardize")
  expData <- ExpressionSet(as.matrix(t(sample_norm)))
  
  # Return
  return(expData)
}


PostProcessMFuzz <- function(mfz) {
  
  # Dependencies
  require(tidyr)
  
  # Convert centroid data into long dataframe
  centroids <- mfz$centers
  centroids_df <- data.frame(centroids)
  centroids_df$cluster <- row.names(centroids_df)
  
  # Convert dataframe into long format
  centroids_long <- tidyr::pivot_longer(centroids_df, names_to = "sample", values_to = "value", 1:length(colnames(centroids)))
  
  # Remove "X" prefix from sample names
  centroids_long$sample <- gsub("X", "", centroids_long$sample)
  
  # Return
  return(centroids_long)
  
}
