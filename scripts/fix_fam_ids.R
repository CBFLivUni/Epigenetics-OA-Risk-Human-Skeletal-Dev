
# Get command line parameters
args <- commandArgs(TRUE)
fam_filename <- args[[1]]
meta_filename <- args[[2]]

# Sort out discrepancies between home 
# setwd("~")
# fam_filename <- ifelse(grepl("Document",getwd()), gsub("[/]home[/][0-9A-Za-z_-]*[/]", "~/../", fam_filename), fam_filename)
# meta_filename <- ifelse(grepl("Document",getwd()), gsub("[/]home[/][0-9A-Za-z_-]*[/]", "~/../", meta_filename), meta_filename)

# Read in metadata
meta_df <- read.table(meta_filename, sep="\t", header=TRUE)[,1:3]
meta_df["ORIGINAL_ID"] <- paste(meta_df$SentrixPosition_A, meta_df$SentrixBarcode_A, sep="_")

# Read in .fam data
fam_df <- read.table(fam_filename, sep=" ", header=FALSE)
fam_df["ORIGINAL_ID"] <- paste(fam_df$V1, fam_df$V2, sep="_")

# Merge
out_df <- merge(meta_df, fam_df, by="ORIGINAL_ID")
out_df <- out_df[,c("Sample_ID","Sample_ID","V3","V4","V5","V6")]
colnames(out_df)[1:2] <- c("V1","V2")

# Save
write.table(out_df, fam_filename, col.names=FALSE, sep="\t", quote=FALSE, row.names=FALSE)
