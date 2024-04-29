


GetRsIDs <- function(in_df, dir2save2="./", fileprefix="snps", filesuffix=".rsids_from_snpdb.tsv", exclude=NA, reGenerateFile=FALSE) {
  
  # Dependencies
  require(SNPlocs.Hsapiens.dbSNP155.GRCh37)
  require(colochelpR)
  require(rsnps)
  
  # Test if file exists, else re-run
  print("##### NOTE THAT THIS ONLY RETAINS SNPS WITH RSIDS")
  annotation_snpfile <- paste0(dir2save2, "/", fileprefix, filesuffix)
  if (!file.exists(annotation_snpfile) | reGenerateFile) {
    
    # Get correct format for chrm-base IDs
    chrbp_df <- in_df[,c("chr","pos")]
    colnames(chrbp_df) <- c("CHR","BP")
    
    # Get SNP IDs
    print("Getting SNP IDs from database......")
    snp_annot.df <- convert_loc_to_rs(chrbp_df, SNPlocs.Hsapiens.dbSNP155.GRCh37)
    colnames(snp_annot.df)[[3]] <- "SNP"
    
    # Keep only alleles with rsIDs
    print(paste0("Excluding this manySNPs due to no rsID = ", unique(nrow(snp_annot.df[!grepl("^rs",snp_annot.df$SNP),]))))
    snp_annot.df <- snp_annot.df[grepl("^rs",snp_annot.df$SNP),]
    
    # If excluding any SNPs by rsID
    if (!is.na(exclude)) { print(paste0("#### NOTE THAT YOU HAVE CHOSEN TO MANUALLY EXCLUDE :- ",paste(exclude, collapse=", "))) }
    snp_annot.df <- snp_annot.df[!snp_annot.df$SNP %in% exclude,]
    
    # Get SNP alleles
    print("Getting SNP alleles......")
    
    # Generate possible alleles
    allNucleotides <- expand.grid(c("A","G","C","T"),c("A","G","C","T"))
    allNucleotides <- allNucleotides[allNucleotides$Var1 != allNucleotides$Var2,]
    allNucleotides <- paste(allNucleotides$Var1, allNucleotides$Var2, sep=":")
    
    # Generate all combinations of ref/allele CHR:POS IDs
    out_df <- data.frame()
    for (nt in allNucleotides) { out_df <- rbind(out_df, cbind(snp_annot.df, CHRPOSID=paste(snp_annot.df$CHR, snp_annot.df$BP, nt, sep=":"))) }
    
    # Retain only those present
    out_df <- out_df[which(out_df$CHRPOSID %in% in_df$CHRPOSID),]
    print(paste0("We have lost ", length(unique(in_df$CHRPOSID)) - length(unique(out_df$CHRPOSID)), " SNPs"))
    
    # Save
    print("Saving.....")
    write.table(out_df, annotation_snpfile, sep="\t", quote=FALSE)
    
    # Else get SNPs from biomart
  } else {
    
    # Read in
    print("Reading in SNP annotation data from file......")
    out_df <- read.table(annotation_snpfile, sep="\t", quote="\"")
    
  }
  
  # Return
  return(out_df)
  
}


SNPChrmPos2rsID <- function(df, snpdb, merge_on="CHRPOS") {
  
  # Dependencies
  require(colochelpR)
  
  # Get granges
  print("## Getting RSIDs......")
  snpid_df <- colochelpR::convert_loc_to_rs(data.frame(CHR=df[!is.na(df$chr),]$chr, 
                                                       BP=df[!is.na(df$chr),]$pos), 
                                            snpdb)
  
  # Merge
  print("## Merging RSIDs with data......")
  snpid_df <- merge(data.frame(snpid_df[,!grepl("CHR|BP",colnames(snpid_df))], CHRPOSID=paste0(snpid_df$CHR,":",snpid_df$BP)), 
                    df[,!grepl("marker",colnames(df))], by=merge_on)
  colnames(snpid_df)[[2]] <- "marker"
  
  # Return SNPs 
  return(snpid_df)
  
}



ReadPLINKtoDF <- function(plinkname, chrs=paste0("chr",1:22), ...) {
  
  # Dependencies
  require(genio)
  
  # NOTE: Don't use argyle here as not only does it take an age akin to the 1st 
  # age of Arda, but it also seems to ignore chrm 20, 21 and 22.
  
  # Read in data
  plink <- genio::read_plink(plinkname)
  # plink <- argyle::read.plink(plinkname, chr=chrs, ...)
  
  # Get as dataframe and rename columns to argyle-style (this is what I 
  # originally built the script for).
  plink.df <- data.frame(merge(plink$bim, plink$X, by.x="id", by.y="row.names"))
  colnames(plink.df)[1:6] <- c("marker","chr","cM","pos","alt","ref")
  # plink.df <- data.frame(plink)
  
  # Return
  return(plink.df)
  
}



CorrTechReps <- function(in_df, meta) {
  
  # Loop over patients and correlate technical replicate genotypes
  cors_df <- do.call("rbind",lapply(as.character(unique(meta$Patient_ID)), function(id) { 
    
    # Subset data
    sub_df <- in_df[,grepl("^X",colnames(in_df)) & gsub("^X|_.*","",colnames(in_df)) == id] 
    
    # Drop NAs
    colnames(sub_df) <- c("Rep1","Rep2")
    sub_df$Rep1 <- as.numeric(sub_df$Rep1)
    sub_df$Rep2 <- as.numeric(sub_df$Rep2)
    sub_df <- sub_df[!is.na(sub_df$Rep1) & !is.na(sub_df$Rep2),]
    
    # Get correlation value
    corr_val <- cor(sub_df$Rep1, sub_df$Rep2, use="pairwise.complete.obs")
    
    # Output
    out_df <- data.frame(ID=id, CORR=corr_val)
    return(out_df) 
    
  }))
  
  # Return
  return(cors_df)
}



GetMismatchedGenoTechRep <- function(in_df, meta) {
  
  # Loop over patients and correlate technical replicate genotypes
  match_df <- do.call("rbind",lapply(as.character(unique(meta$Patient_ID)), function(id) { 
    
    # Subset data
    sub_df <- in_df[,grepl("^X",colnames(in_df)) & gsub("^X|_.*","",colnames(in_df)) == id] 
    sub_df <- cbind(MARKER=in_df$marker, sub_df)
    
    # Drop NAs
    colnames(sub_df) <- c("MARKER","Rep1","Rep2")
    sub_df <- sub_df[!is.na(sub_df$Rep1) & !is.na(sub_df$Rep2),]
    
    # Filter only inconsistent
    sub_df <- sub_df[sub_df$Rep1 != sub_df$Rep2,]
    
    # If any inconsistent
    if (nrow(sub_df) > 0) { 
      
      # Return only inconsistent
      inc_df <- cbind(ID=id, sub_df)
      return(inc_df)
      
      # Else just return NULL
    } else { return(NULL) }
    
  }))
  
  # Return
  return(match_df)
  
}