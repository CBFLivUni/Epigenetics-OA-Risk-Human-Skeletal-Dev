
# Input variables
MAIN_DIR=/mnt/f/SarahRice_mQTL/
SCRIPT_DIR=$MAIN_DIR/scripts/
MERGED_DIR=$MAIN_DIR/qced_data/merged/
IMPUT_DIR=$MAIN_DIR/imputation/
FASTA_FILE=$MAIN_DIR/genome/GRCh37/human_g1k_v37.fasta
CHECKVCF_SCRIPT=$MAIN_DIR/utilities/genotyping/imputation/checkVCF/checkVCF.py

# Sort PARS region
plink --keep-allele-order \
      --merge-x no-fail \
      --bfile $MERGED_DIR/all_individuals.merged.q_filt \
      --make-bed --out $IMPUT_DIR/input/all_individuals.merged.merge_pars_x

# Recode to VCF
plink --keep-allele-order \
      --recode vcf \
      --bfile $IMPUT_DIR/input/all_individuals.merged.merge_pars_x \
      --out $IMPUT_DIR/input/all_individuals.merged.for_imputation

# Sort chromosomes
sed 's/^23/X/g; s/^24/Y/g; s/^25/PARS/g; s/^26/MT/g; s/ID=23/ID=X/g; s/ID=24/ID=Y/g; s/ID=25/ID=PARS/g; s/ID=26/ID=MT/g' \
    $IMPUT_DIR/input/all_individuals.merged.for_imputation.vcf \
    > /tmp/all_individuals.merged.for_imputation.vcf && \
    mv /tmp/all_individuals.merged.for_imputation.vcf \
    $IMPUT_DIR/input/all_individuals.merged.for_imputation.vcf

# Run VCF checks
conda activate python2
python $CHECKVCF_SCRIPT \
       -r $FASTA_FILE \
       -o $IMPUT_DIR/input/vcfcheck \
       $IMPUT_DIR/input/all_individuals.merged.for_imputation.vcf

# Zip
conda activate plink
bcftools view $IMPUT_DIR/input/all_individuals.merged.for_imputation.vcf \
         -Oz -o $IMPUT_DIR/input/all_individuals.merged.for_imputation.vcf.gz

# Split by chromosomes
bash split_all_vcfs.sh \
    *for_imputation*.vcf.gz \
    $IMPUT_DIR/input/ $IMPUT_DIR/input/
