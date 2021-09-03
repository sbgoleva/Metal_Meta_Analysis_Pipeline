#module --ignore-cache load "LLVM/.7.0.1" GCC/8.2.0 OpenMPI/3.1.4 R/3.6.0
#R
.libPaths("~/R/rlib-3.6.0")
library(qqman)

#load in meta analysis
meta_gwas = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/metal_fs_gwas_1.tbl", stringsAsFactors = FALSE, header=TRUE, comment.char="", sep="\t")

#load in each individual sum stat
snp_ref_cc = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/CC/CC_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="")
snp_ref_ipsych = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/iPSYCH/iPSYCH_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="")
snp_ref_biovu = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/BioVU/BioVU_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="", quote="")
snp_ref_mvp = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/MVP/MVP_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="", quote="")
snp_ref_partners = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/MGBB/MGBB_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="")
snp_ref_mtsinai = read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/BioMe/BioMe_autosome_sumstats_for_meta_analysis_cleaned.txt", stringsAsFactors = FALSE, header=TRUE, comment.char="")

#rbind all individual sumstats and remove duplicates
snp_ref_all = rbind(snp_ref_cc[c("RS_ID", "CHR", "BP")], snp_ref_ipsych[c("RS_ID", "CHR", "BP")], snp_ref_mtsinai[c("RS_ID", "CHR", "BP")], snp_ref_biovu[c("RS_ID", "CHR", "BP")],snp_ref_mvp[c("RS_ID", "CHR", "BP")],snp_ref_mgbb[c("RS_ID", "CHR", "BP")])
snp_ref_all=unique(snp_ref_all)

#merge meta and individual sumstats (for CHR and BP info)
meta_snp_chr_bp = merge(meta_gwas, snp_ref_all, by.x='MarkerName', by.y="RS_ID", all.x=TRUE, all.y=FALSE)
meta_snp_chr_bp$CHR=as.numeric(meta_snp_chr_bp$CHR)

###exclude single SNPs
meta_snp_chr_bp=meta_snp_chr_bp[as.numeric(meta_snp_chr_bp$HetDf)!=0,]

##write out sumstats (this will exclude any SNPs from just 1 biobank, so use this for follow up analyses
write.table(meta_snp_chr_bp, "/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/metal_fs_gwas_filtered_20210831.tbl", sep="\t",row.names=FALSE, quote=FALSE)

#plot Manhattan
manhattan_name = "/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/metal_fs_gwas_filtered_20210831_manhattan.png"
png(manhattan_name, width=5000, height=1440, res=360, type="cairo")
manhattan(meta_snp_chr_bp, snp = 'snpid', chr = "CHR", bp = "BP", p="P.value", cex=0.6, col=c("#000000","#a9a8a9"), chrlabs=c(1:22))
dev.off()

#plot QQ plot
qq_name = "/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/metal_fs_gwas_filtered_20210831_qqplot.png"
png(qq_name, width=2880, height=1440, res=360, type="cairo")
qq(meta_snp_chr_bp$P.value)
dev.off()
