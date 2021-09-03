#module --ignore-cache load "LLVM/.7.0.1" GCC/8.2.0 OpenMPI/3.1.4 R/3.6.0
#R
.libPaths("~/R/rlib-3.6.0")
library(qqman)


meta_gwas = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/metal_fs_gwas_1.tbl", stringsAsFactors = FALSE, header=TRUE, comment.char="", sep="\t")

snp_ref = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/CC/CC_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="")
snp_ref2 = snp_ref[c("chr_pos", "CHR", "BP","dbSNPid")]
colnames(snp_ref2)=c("MarkerName","CHR","BP","snpid")

snp_ref_ipsych = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/iPSYCH/iPSYCH_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="")
snp_ref_ipsych2 = snp_ref_ipsych[c("chr_pos", "CHR", "BP","SNP")]
colnames(snp_ref_ipsych2)=c("MarkerName","CHR","BP","snpid")

snp_ref_biovu = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/BioVU/BioVU_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="", quote="")
snp_ref_biovu2 = snp_ref_biovu[c("chr_pos","CHR","POS","SNPID")]
colnames(snp_ref_biovu2)=c("MarkerName","CHR","BP","snpid")

snp_ref_mvp = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/MVP/MVP_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="", quote="")
snp_ref_mvp2 = snp_ref_mvp[c("chr_pos","CHR","POS","SNP")]
colnames(snp_ref_mvp2)=c("MarkerName","CHR","BP","snpid")

snp_ref_partners = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/MGBB/MGBB_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="")
snp_ref_partners2 = snp_ref_partners[c("chr_pos","CHR","POS","SNPID_rs")]
colnames(snp_ref_partners2)=c("MarkerName","CHR","BP","snpid")

snp_ref_mtsinai = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/BioMe/BioMe_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="")
snp_ref_mtsinai2 = snp_ref_mtsinai[c("chr_pos","CHR","BP","SNP")]
colnames(snp_ref_mtsinai2) = c("MarkerName","CHR","BP","snpid")

snp_ref_3 = rbind(snp_ref2, snp_ref_ipsych2, snp_ref_mtsinai2, snp_ref_biovu2,snp_ref_mvp2,snp_ref_partners2)
snp_ref_3=unique(snp_ref_3)


meta_snp_chr_bp = merge(meta_gwas, snp_ref_3, by.x='MarkerName', by.y="MarkerName", all.x=TRUE, all.y=FALSE)
nrow(meta_snp_chr_bp)


meta_snp_chr_bp$CHR=as.numeric(meta_snp_chr_bp$CHR)
nrow(meta_snp_chr_bp)



###exclude single SNPs
meta_snp_chr_bp=meta_snp_chr_bp[as.numeric(meta_snp_chr_bp$HetDf)!=0,]
nrow(meta_snp_chr_bp)

write.table(meta_snp_chr_bp, "/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/metal_fs_gwas_filtered_20210831.tbl", sep="\t",row.names=FALSE, quote=FALSE)


manhattan_name = "/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/metal_fs_gwas_filtered_20210831_manhattan.png"
png(manhattan_name, width=5000, height=1440, res=360, type="cairo")
manhattan(meta_snp_chr_bp, snp = 'snpid', chr = "CHR", bp = "BP", p="P.value", cex=0.6, col=c("#000000","#a9a8a9"), chrlabs=c(1:22))
dev.off()

qq_name = "/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/metal_fs_gwas_filtered_20210831_qqplot.png"
png(qq_name, width=2880, height=1440, res=360, type="cairo")
qq(meta_snp_chr_bp$P.value)
dev.off()
