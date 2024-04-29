
# Get command line parameters
args <- commandArgs(TRUE)
filename <- args[[1]]

# Get output filename
outname <- paste0(filename, ".rsids.tsv")

# Load library
require(SNPlocs.Hsapiens.dbSNP144.GRCh37)
require(mgsub)

# Read in SNP dataframe
snp_df <- read.table(filename, header=FALSE, sep="\t")

# Get corret chrm and position
snp_df["chr"] <- gsub(":.*","",snp_df$V2)
snp_df["pos"] <- gsub(".*:","",gsub(":[A-Z][-].*","",snp_df$V2))

# Switch out awkward chromosome names
snp_df[snp_df$chr == "23",]["chr"] <- "X"
snp_df[snp_df$chr == "24",]["chr"] <- "Y"
snp_df[snp_df$chr == "26",]["chr"] <- "MT"

# Convert chromosome positions to rsIDs 
snpid_df <- colochelpR::convert_loc_to_rs(data.frame(CHR=snp_df[!is.na(snp_df$chr),]$chr, 
                                                     BP=snp_df[!is.na(snp_df$chr),]$pos), 
                                          SNPlocs.Hsapiens.dbSNP144.GRCh37)
snpid_df$CHR <- as.character(snpid_df$CHR)

# Remove duplicate SNPs
snpid_df <- snpid_df[!duplicated(snpid_df$SNP) & !is.na(snpid_df$SNP),]

# Switch out awkward chromosome names
snpid_df[snpid_df$CHR == "X",]["CHR"] <- "23"
snpid_df[snpid_df$CHR == "Y",]["CHR"] <- "24"
snpid_df[snpid_df$CHR == "MT",]["CHR"] <- "26"

# Save
write.table(snpid_df, outname, sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
