bash sidestep_genomestudio.grch37.sh

bash replace_xymt.sh 20*.vcf ../rawdata/vcf_files/

bash vcf2plink.sh 20*numeric*.vcf ../rawdata/vcf_files/ ../qced_data/plink/

bash edit_ids.sh *.bed ../qced_data/plink/ ../rawdata/E181906_SarahRice_GSA_290323/E181906_SarahRice_GSA_290323_SamplesTable.txt

bash merge_plinks.sh *.bed all_individuals.merged ../qced_data/plink/ ../qced_data/merged/

plink --keep-allele-order --bfile ../qced_data/merged/all_individuals.merged --update-sex ../rawdata/sample_metadata.for_plink.tsv --covar ../rawdata/sample_metadata.for_plink.tsv --make-bed --out ../qced_data/merged/all_individuals.merged.with_covar

plink --bfile ../qced_data/merged/all_individuals.merged.with_covar --freq --out ../qced_data/merged/qc_stats/
plink --bfile ../qced_data/merged/all_individuals.merged.with_covar --missing --out ../qced_data/merged/qc_stats/
plink --bfile ../qced_data/merged/all_individuals.merged.with_covar --hardy --out ../qced_data/merged/qc_stats/

plink --keep-allele-order --nonfounders --split-x b37 --maf 0.05 --geno 0.025 --mind 0.025 --hwe 0.000001 --bfile ../qced_data/merged/all_individuals.merged.with_covar --make-bed --out ../qced_data/merged/all_individuals.merged.q_filt