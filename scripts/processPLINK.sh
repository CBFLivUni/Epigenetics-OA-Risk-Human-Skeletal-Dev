
# Input variables
MAC_FILT=1
GENO_FILT=0.025
MIND_FILT=0.025
HWE_FILT=0.000001
MAIN_DIR=/mnt/f/SarahRice_mQTL/
SCRIPT_DIR=$MAIN_DIR/scripts/
PLINK_DIR=$MAIN_DIR/qced_data/plink/
MERGED_DIR=$MAIN_DIR/qced_data/merged/
SAMPLE_DATA=$MAIN_DIR/rawdata/sample_metadata.for_plink.tsv

# Make output directories
mkdir -p $MERGED_DIR/qc_stats/
mkdir -p $MERGED_DIR/IBD/
mkdir -p $MERGED_DIR/pca/

# Update sample IDs
bash $SCRIPT_DIR/edit_ids.sh \
     *xymt*.bed \
     $PLINK_DIR \
     $SAMPLE_DATA

# Merge PLINK files - this takes a while!
bash merge_plinks.sh \
     *xymt*.bed \
     all_individuals.merged \
     $PLINK_DIR \
     $MERGED_DIR

# Add sample metadata
plink \
      --keep-allele-order \
      --bfile $MERGED_DIR/all_individuals.merged \
      --update-sex $SAMPLE_DATA \
      --covar $SAMPLE_DATA \
      --make-bed --out $MERGED_DIR/all_individuals.merged.with_covar

# Calculate pre-filtration statistics
plink \
      --bfile $MERGED_DIR/all_individuals.merged.with_covar \
      --freq \
      --out $MERGED_DIR/qc_stats/preFilt
plink \
      --bfile $MERGED_DIR/all_individuals.merged.with_covar \
      --missing \
      --out $MERGED_DIR/qc_stats/preFilt
plink \
      --bfile $MERGED_DIR/all_individuals.merged.with_covar \
      --hardy \
      --out $MERGED_DIR/qc_stats/preFilt

# Filter data
plink \
     --keep-allele-order --nonfounders \
     --split-x b37 \
     --mac $MAC_FILT --geno $GENO_FILT --mind $MIND_FILT --hwe $HWE_FILT \
     --bfile $MERGED_DIR/all_individuals.merged.with_covar \
     --make-bed --out $MERGED_DIR/all_individuals.merged.q_filt
     
# Calculate pre-filtration statistics
plink \
      --bfile $MERGED_DIR/all_individuals.merged.q_filt \
      --freq \
      --out $MERGED_DIR/qc_stats/postFilt
plink \
      --bfile $MERGED_DIR/all_individuals.merged.q_filt \
      --missing \
      --out $MERGED_DIR/qc_stats/postFilt
plink \
      --bfile $MERGED_DIR/all_individuals.merged.q_filt \
      --hardy \
      --out $MERGED_DIR/qc_stats/postFilt

# Inter-relatedness check
plink --keep-allele-order --genome \
      --bfile $MERGED_DIR/all_individuals.merged.q_filt \
      --out $MERGED_DIR/IBD/ibd_report

# Run PCA and MDS
plink --keep-allele-order --pca \
      --bfile $MERGED_DIR/all_individuals.merged.q_filt \
      --out $MERGED_DIR/pca/all_individuals.pca
plink --keep-allele-order --cluster --mds-plot 20 \
      --bfile $MERGED_DIR/all_individuals.merged.q_filt \
      --out $MERGED_DIR/pca/all_individuals.mds