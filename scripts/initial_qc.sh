bash vcf2plink.sh *numeric*.vcf ../rawdata/vcf_files/ ../qced_data/plink/

bash edit_ids.sh *.bed ../qced_data/plink/ ../rawdata/E181906_SarahRice_GSA_290323/E181906_SarahRice_GSA_290323_SamplesTable.txt

bash merge_plinks.sh *.bed all_individuals.merged ../qced_data/plink/ ../qced_data/merged/

plink --keep-allele-order --bfile ../qced_data/merged/all_individuals.merged --update-sex ../rawdata/sample_metadata.for_plink.tsv --covar ../rawdata/sample_metadata.for_plink.tsv --make-bed --out ../qced_data/merged/all_individuals.merged.with_covar

plink --bfile ../qced_data/merged/all_individuals.merged.with_covar --freq --out ../qced_data/merged/qc_stats/
plink --bfile ../qced_data/merged/all_individuals.merged.with_covar --missing --out ../qced_data/merged/qc_stats/
plink --bfile ../qced_data/merged/all_individuals.merged.with_covar --hardy --out ../qced_data/merged/qc_stats/

plink --keep-allele-order --nonfounders --split-x b37 --maf 0.05 --geno 0.025 --mind 0.025 --hwe 0.000001 --bfile ../qced_data/merged/all_individuals.merged.with_covar --make-bed --out ../qced_data/merged/all_individuals.merged.q_filt



plink  --keep-allele-order --check-sex --bfile ../qced_data/merged/all_individuals.merged.q_filt --make-bed --out ../qced_data/merged/all_individuals.merged.q_filt

plink  --keep-allele-order --genome --bfile ../qced_data/merged/all_individuals.merged.q_filt --out ../qced_data/IBD/ibd_report

plink  --keep-allele-order --pca --bfile ../qced_data/merged/all_individuals.merged.q_filt --out ../qced_data/merged/pca/all_individuals.pca
plink  --keep-allele-order --cluster --mds-plot 20 --bfile ../qced_data/merged/all_individuals.merged.q_filt --out ../qced_data/merged/pca/all_individuals.mds

plink --keep-allele-order --bfile ../qced_data/merged/all_individuals.merged.q_filt --recode vcf --out ../qced_data/merged/all_individuals.merged.q_filt

sed 's/^23/X/g; s/^24/Y/g; s/^25/PARS/g; s/^26/MT/g' ../qced_data/merged/all_individuals.merged.q_filt.vcf > /tmp/all_individuals.merged.q_filt.vcf && mv /tmp/all_individuals.merged.q_filt.vcf ../qced_data/merged/all_individuals.merged.q_filt.vcf





python ../utilities/genotyping/imputation/checkVCF/checkVCF.py -r ../utilities/genome/GRCh37/human_g1k_v37.fasta -o vcfcheck ../qced_data/merged/all_individuals.merged.q_filt.vcf

bcftools view ../qced_data/merged/all_individuals.merged.q_filt.vcf -Oz -o ../qced_data/merged/all_individuals.merged.q_filt.vcf.gz
bcftools index ./qced_data/merged/all_individuals.merged.q_filt.vcf.gz
bash split_all_vcfs.sh *individual*.vcf.gz ../qced_data/merged/ ../results/imputation/input/



bash flip_strand.sh all_individuals.merged.q_filt.bim ../utilities/genome/GRCh37/human_g1k_v37.xymt_numeric.fasta ../qced_data/merged/ ../qced_data/strandflip/

plink --keep-allele-order --indep-pairwise 1000 50 0.2 --bfile ../qced_data/strandflip/all_individuals.merged.q_filt.strand_flip --make-bed --out ../qced_data/merged/all_individuals.merged.ld_filt



