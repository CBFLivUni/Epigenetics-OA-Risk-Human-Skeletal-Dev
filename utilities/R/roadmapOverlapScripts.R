
ReadInOverlapDMRs <- function(featursOfInterest, dmrSet, in_dir="./") {
  
  # Loop over ids and read in overlapped feature sets
  ovlps <- lapply(featursOfInterest, function(feat) { 
    
    # Read in
    print(feat)
    out_df <- try(read.table(paste0(in_dir, dmrSet, "/overlap",feat, "/unique_",dmrSet,"_",feat,"_overlaps.bed"), sep="\t", header=FALSE))
    
    # If any overlaps
    if (is.data.frame(out_df) > 0) {
      
      # Add feature
      out_df <- cbind(out_df, FEATURE=feat)
      
      # Return
      return(out_df)
      
      # Else return NULL
    } else { 
      
      # Return
      return(NULL)
      
    }
    
  })
  
  # Combine
  ovlp_df <- do.call("rbind",ovlps)
  
  # Return
  return(ovlp_df)
  
}


CalcFreqFeatures <- function(df, set) {
  
  # Get background frequencies
  out_df <- data.frame(table(df$FEATURE))
  colnames(out_df) <- c("Feature","Count")
  
  # Add proportions
  out_df["Prop"] <- out_df$Count / sum(out_df$Count)
  
  # Add that it is background
  out_df <- cbind(Set=set, out_df)
  
  # Return
  return(out_df)
}

CalcFeatureWidths <- function(df, set, startCol="V2", endCol="V3") {
  
  # Get background frequencies
  df["WIDTH"] <- df[[endCol]] - df[[startCol]] 
  
  # Look over features
  df["Count"] <- 0
  for (feat in unique(df$FEATURE)) {
    
    # Add proportions
    df[df$FEATURE == feat,]["Count"] <- sum(df[df$FEATURE == feat,]$WIDTH)
    
    
  }
  
  # Add to output
  out_df <- cbind(Set=set, df[!duplicated(df[,c("FEATURE")]),c("FEATURE","Count")])
  
  # Calculate proportion
  out_df["Prop"] <- out_df$Count / sum(out_df$Count)
  
  # Return
  colnames(out_df)[[2]] <- "Feature"
  rownames(out_df) <- NULL
  return(out_df)
}


TestChiSqEnrichment <- function(inDF, testFeature, backFeature, ...) {
  
  require(broom)
  
  # Subset features
  testDF <- inDF[inDF$Set == testFeature,]; rownames(testDF) <- testDF$Feature
  backDF <- inDF[inDF$Set == backFeature,]; rownames(backDF) <- backDF$Feature
  
  # Get common
  commonFeatures <- intersect(testDF$Feature, backDF$Feature)
  
  # Combine
  testDF <- cbind(testDF[commonFeatures,"Count",drop=FALSE], 
                  backDF[commonFeatures,"Count",drop=FALSE])
  
  # Test and get output
  outDF <- broom::tidy(chisq.test(testDF, ...))
  outDF <- cbind(Test=testFeature, Background=backFeature, outDF)
  
  # Return
  return(outDF)
  
}


FeatureWiseFishers <- function(inDF, testFeature, backFeature, maxItems, ...) {
  
  # Dependencies
  require(broom)
  
  # Split to test/background
  ovlpTestData <- inDF[inDF$Set == testFeature,]
  ovlpBackData <- inDF[inDF$Set == backFeature,]
  
  # Get common features
  commonFeatures <- intersect(ovlpTestData$Feature, ovlpBackData$Feature)
  
  # Loop over features and run tests
  allFishDF <- do.call("rbind",lapply(commonFeatures, function(feature) {
    
    # Track feature
    print(feature)
    
    # Subset
    subOvlpBackData <- ovlpBackData[ovlpBackData$Feature == feature,]
    subOvlpTestData <- ovlpTestData[ovlpTestData$Feature == feature,]
    
    # If we want to look for number of features
    testM <- rbind(cbind(subOvlpBackData$Count, sum(ovlpBackData$Count) - subOvlpBackData$Count),
                   cbind(subOvlpTestData$Count, sum(ovlpTestData$Count) - subOvlpTestData$Count))
    
    # Get as proportion
    testP <- testM / rowSums(testM)
    
    # Run test and result
    fishDF <- broom::tidy(fisher.test(testM))
    fishDF <- cbind(Test=testFeature, Background=backFeature, Feature=feature, fishDF)
    fishDF <- cbind(fishDF, 
                    t(data.frame(c(testM[1,]))), t(data.frame(c(testM[2,]))),
                    t(data.frame(c(testP[1,]))), t(data.frame(c(testP[2,]))))
    colnames(fishDF)[10:17] <- c("FreqBackYesOverlap","FreqBackNoOverlap","FreqTestYesOverlap","FreqTestNoOverlap",
                                 "PropBackYesOverlap","PropBackNoOverlap","PropTestYesOverlap","PropTestNoOverlap")
    
    # Add fold-enrichment
    fishDF["FoldEnrichment"] <- fishDF$PropTestYesOverlap / ( fishDF$PropTestYesOverlap + fishDF$PropTestNoOverlap ) /
      fishDF$PropBackYesOverlap / ( fishDF$PropBackYesOverlap + fishDF$PropBackNoOverlap )
    # Return
    return(fishDF)
    
  }))
  
  # Rename p-value and tidy rownames
  colnames(allFishDF) <- gsub("p[.]value","rawPValue", colnames(allFishDF))
  rownames(allFishDF) <- NULL
  
  # Return
  return(allFishDF)
  
}