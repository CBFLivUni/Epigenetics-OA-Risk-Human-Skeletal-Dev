
ReadInOverlapDMRs <- function(featuresOfInterest, dmrSet, in_dir="./") {
  
  # Loop over ids and read in overlapped feature sets
  ovlps <- lapply(featuresOfInterest, function(feat) { 
    
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

CalcFreqFeatures <- function(df, set, cols="FEATURE") {
  
  # Get background frequencies
  out_df <- data.frame(table(df[,cols]))
  colnames(out_df) <- c(cols,"Count")
  
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
  
  browser()
  
  # Split to test/background
  ovlpTestData <- inDF[inDF$Set == testFeature,]
  ovlpBackData <- inDF[inDF$Set == backFeature,]
  
  # Get common features
  commonFeatures <- intersect(ovlpTestData$FEATURE, ovlpBackData$FEATURE)
  
  # Loop over features and run tests
  allFishDF <- do.call("rbind",lapply(commonFeatures, function(feature) {
    
    # Track feature
    print(feature)
    
    # Subset
    subOvlpBackData <- ovlpBackData[ovlpBackData$FEATURE == feature,]
    subOvlpTestData <- ovlpTestData[ovlpTestData$FEATURE == feature,]
    
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


PlotDMROverlapProportions <- function(df, fileNamePrefix, colPal, outDir="./") {
  
  # Properly order
  df$Set <- factor(df$Set, levels=rev(c("All Sig DMRs", "Up Sig DMRs", "Down Sig DMRs", "Genome")))
  
  # Plot
  DMR_plt <- ggplot(df, aes(y=Set, x=Prop, fill=FEATURE)) + 
    geom_bar(stat="identity", color="black") +
    scale_fill_manual(values=colPal) +
    ylab("") + xlab("Proportion of States") +
    ggtitle("Developmental Stage DMRs") +
    theme_minimal(base_size=20) + theme(plot.title=element_text(size=36),
                                        axis.title=element_text(size=32),
                                        legend.text=element_text(size=16),
                                        legend.title=element_blank(),
                                        panel.grid.major.y = element_blank() ,
                                        panel.grid.major.x = element_line( size=.1, color="darkgrey" ))
  
  # Save
  ggsave(here(outDir, paste0(fileNamePrefix, "_freqBarPlot.png")), DMR_plt, width=4000, height=2000, units="px", device="png",bg="white")
  ggsave(here(outDir, paste0(fileNamePrefix, "_freqBarPlot.pdf")), DMR_plt, width=4000, height=2000, units="px", device="pdf",bg="white")
  
  # Return
  return(DMR_plt)
  
}
