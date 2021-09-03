module spider Foo/11.1
ml Anaconda2/4.4.0
cd /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/pnes_saige_ldsc/ldsc
source activate ldsc
conda info --envs
module spider Anaconda
conda env create --file environment.yml
source activate ldsc

##running FS phenotype
/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/pnes_saige_ldsc/ldsc/munge_sumstats.py \
  --out /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/METAL_BioVU_UKBBR56_CC/ldsc_meta_pnes \
  --merge-alleles /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/pnes_saige_ldsc/ldsc/w_hm3.snplist \
  --N 494008 \
  --N-cas 7897 \
  --N-con 486111 \
  --sumstats /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/metal_fs_gwas_filtered_20210831.tbl \
  --snp snpid \
  --N-col Weight
  
##calculate heritability for phenotype
/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/pnes_saige_ldsc/ldsc/ldsc.py \
  --h2 /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/METAL_BioVU_UKBBR56_CC/ldsc_meta_pnes.sumstats.gz \
  --ref-ld-chr /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/pnes_saige_ldsc/ldsc/eur_w_ld_chr/ \
  --out /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/METAL_BioVU_UKBBR56_CC/ldsc_meta_pnes_h2 \
  --w-ld-chr /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/pnes_saige_ldsc/ldsc/eur_w_ld_chr/ \
  --samp-prev 0.015985571 \
  --pop-prev 0.0014
  
##calculate rg for related phenotypes (must munge sumstats separately for each)
/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/pnes_saige_ldsc/ldsc/ldsc.py \
--rg /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/METAL_BioVU_UKBBR56_CC/ldsc_FE.sumstats.gz,/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/METAL_BioVU_UKBBR56_CC/ldsc_meta_pnes.sumstats.gz \
--ref-ld-chr /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/pnes_saige_ldsc/ldsc/eur_w_ld_chr/ \
--w-ld-chr /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/pnes_saige_ldsc/ldsc/eur_w_ld_chr/ \
--out /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/METAL_BioVU_UKBBR56_CC/pnes_FE_corr 
