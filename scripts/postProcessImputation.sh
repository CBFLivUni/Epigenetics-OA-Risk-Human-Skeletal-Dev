
# Input variables
MAC_FILT=2
GENO_FILT=0.025
MIND_FILT=0.025
HWE_FILT=0.000001

# Input variables
MAIN_DIR=/mnt/f/SarahRice_mQTL/
SCRIPT_DIR=$MAIN_DIR/scripts/
MERGED_DIR=$MAIN_DIR/qced_data/merged/
IMPUT_DIR=$MAIN_DIR/imputation/
IC_SCRIPT=$MAIN_DIR/utilities/genotyping/imputation/ic.v1.0.9/IC/ic.pl
FASTA_FILE=$MAIN_DIR/genome/GRCh37/human_g1k_v37.fasta
HRC_SITES=$MAIN_DIR/utilities/genotyping/imputation/hrc/HRC.r1-1.GRCh37.wgs.mac5.sites.tab
SAMPLE_DATA=$MAIN_DIR/rawdata/sample_metadata.for_plink.tsv

# Make outputs
mkdir -p $IMPUT_DIR/imputed/qc/

# Run IC imputation QC
perl $IMPUT_DIR \
    -r $HRC_SITES \
    -d $IMPUT_DIR/imputed/ \
    -o $IMPUT_DIR/imputed/qc/

# Preprocess imputed VCFs
bash preprocess_imputed_vcfs.sh \
    chr*.vcf.gz \
    $FASTA_FILE \
    $IMPUT_DIR/imputed/ \
    hrc \

# Remove SNPs that fail imputation
zcat $IMPUT_DIR/imputed/imputed.hrc.vcf.gz | awk 'IFS="\t"; $7 == PASS {print $7}' | grep -v "^#" | awk '$0 != "" {print $0}' > $IMPUT_DIR/imputed/imputed.hrc.vcf.noheader
cat $IMPUT_DIR/imputed/imputed.hrc.vcf.noheader | sed 's/.*R2=//g; s/\s.*//g' | sort  | head # Check min R2
zcat $IMPUT_DIR/imputed/imputed.hrc.vcf.gz | grep "^#" > $IMPUT_DIR/imputed/imputed.hrc.vcf.header
cat $IMPUT_DIR/imputed/imputed.hrc.vcf.header $IMPUT_DIR/imputed/imputed.hrc.vcf.noheader | awk '$0 != "" {print $0}' > $IMPUT_DIR/imputed/imputed.hrc.imppass_filt.vcf

# Convert VCF to PLINK
bash vcf2plink.sh \
    imputed*.imppass_filt.vcf \
    $IMPUT_DIR/imputed/ \
    $IMPUT_DIR/plink/

# Add covariates
plink \
    --keep-allele-order \
    --bfile $IMPUT_DIR/plink/imputed.hrc.imppass_filt.plink \
    --update-sex ../rawdata/sample_metadata.for_plink.tsv \
    --covar $SAMPLE_DATA \
    --make-bed --out $IMPUT_DIR/plink/imputed.hrc.plink.with_covar

# Filter
plink \
    --keep-allele-order \
    --nonfounders \
    --mac $MAC_FILT --geno $GENO_FILT --mind $MIND_FILT --hwe $HWE_FILT \
    --bfile $IMPUT_DIR/plink/imputed.hrc.plink.with_covar \
    --make-bed --out $IMPUT_DIR/plink/imputed.hrc.plink.q_filt

# IBD checks
plink \
    --keep-allele-order \
    --genome \
    --bfile $IMPUT_DIR/plink/imputed.hrc.plink.q_filt \
    --out $IMPUT_DIR/plink/imputed.ibd_report

# Flip strand
bash flip_strand.sh \
    imputed.hrc.plink.q_filt.bim \
    $FASTA_FILE \
    $IMPUT_DIR/plink/ \
    $IMPUT_DIR/plink/

# PCA and MDS
plink \
    --keep-allele-order \
    --pca \
    --bfile $IMPUT_DIR/plink/imputed.hrc.plink.q_filt.strand_flip \
    --out $IMPUT_DIR/plink/pca/imputed.pca
plink \
    --keep-allele-order \
    --cluster --mds-plot 20 \
    --bfile $IMPUT_DIR/plink/imputed.hrc.plink.q_filt.strand_flip \
    --out $IMPUT_DIR/plink/pca/imputed


