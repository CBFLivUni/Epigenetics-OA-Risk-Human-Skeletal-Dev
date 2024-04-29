

# Function to extract limma results, annotate with probe information and print # sig genes
GetLimmaResults <- function(fit, cont_m, annot_df, padj_thr=0.05, trt_toggle=FALSE) {
  
  # Dependencies
  require(limma)
  
  # Initialise empty output
  probes.df <- data.frame()
  
  # View top table
  ## Loop over BMI group contrasts
  for (coef_name in colnames(cont_m)) {
    
    # Toggle if using treat or not (Basically set this to anything but FALSE)
    if (trt_toggle) { 
      
      # Get results
      df <- topTreat(fit, coef=coef_name, adjust="BH", sort.by="logFC", number=Inf) 
      
    } else {
      
      # Get results
      df <- topTable(fit, coef=coef_name, adjust="BH", sort.by="B", number=Inf) 
      
    }
    
    # Get number sig
    print(paste0("### Contrasting: ",coef_name," // Sig probes = ", nrow(df[df$adj.P.Val < padj_thr,])))
    
    # Merge with annotations
    df <- merge(annot_df[,c("chr","pos","strand","Name","ProbeSeqA","ProbeSeqB","Type","NextBase","Probe_rs","Probe_maf","CpG_rs","CpG_maf", 
                            "Islands_Name", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")], 
                df, 
                by="Name", by.y="row.names")
    df["CONTRAST"] <- coef_name
    
    # Add to larger dataframe
    probes.df <- data.frame(rbind(df,probes.df))
    
  }
  
  # Return
  return(probes.df)
  
}