#module --ignore-cache load "LLVM/.7.0.1" GCC/8.2.0 OpenMPI/3.1.4 R/3.6.0
#R
.libPaths("~/R/rlib-3.6.0")

####QCing GWAS SumStats prior to meta analysis###
##To Be Done with each individual summary statistic file##
#file info to fill out#
ss_file="/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/20210412_PNES_GWAS_Redo_Matched_EUR_step2_chr.all.txt"
sumstat_name="BioVU"
sep_type=" "
output_dir="/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/"

#column info to fill out#
is_n_col=TRUE #is there a sample size column present
n_col="N"
sample_size=NA #only need this if is_n_col=FALSE
EAF_col="AF_Allele2" #EA=effect allele
EA_col="Allele2"
OA_col="Allele1"
P_col="p.value"
effect_col="BETA"
is_effect_BETA=TRUE #TRUE if Beta, FALSE if OR
se_col="SE"
#if using BioVU sumstats and not yet added Rsq info to the file, set add_imp_data_for_biovu to TRUE, is_rsq to TRUE, and imp_column to "Rsq"
add_imp_data_for_biovu=TRUE #will be added as 'Rsq' column to sumstats
is_rsq=TRUE #TRUE if Rsq, FALSE if INFO score
imp_column="Rsq"
#only need chr_col and bp_col if add_imp_data_for_biovu==TRUE
chr_col="CHR"
bp_col="POS"


lower_EAF_bound=0.05 #desired EAF filters
upper_EAF_bound=0.95

##EAF reference file info. Below is using 1kg snp list in European Ancestry
ref_eaf_file="/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Allele_Freqs/1000GP_p3v5_legends_rbind.noDup.noMono.noCnv.noCnAll.afref.EUR.txt.gz"
ref_eaf_col="eaf"
ref_chrbp_col="cptid"
ref_ea_col="ea"
ref_oa_col="oa"

is_first_sumstat_analysed=TRUE

is_x_present=FALSE
x_col=""
####################################################################



library(ggplot2)

if(!dir.exists(paste0(output_dir,sumstat_name))){
	dir.create(paste0(output_dir,sumstat_name))
}
setwd(paste0(output_dir,sumstat_name))

#create or add to summary document
if(is_first_sumstat_analysed){
	summary.dat<-data.frame("Filter"=c("Initial # SNPs",
		paste0("EAF <",upper_EAF_bound,", >",lower_EAF_bound),
		"Complete cases",
		"Keep only A,C,T,G alleles",
		"Remove non-numeric Betas",
		"Remove negative, non-numeric StErr",
		"# BioVU snps after merging Rsq info",
		"Rsq>0.3,INFO>0.8",
		"Rsq<1,INFO<1",
		"# SNPs after merging with Ref EAF file",
		"Initial # SNPs where SS EA == Ref EA",
		"Initial # SNPs where SS EA == Ref OA",
		"Switched EA and OA cols",
		"Final # SNPs where SS EA == Ref EA",
		"Exclude AEF outside 10% of ref AF",
		"# autosomal SNPs",
		"# xchr SNPs"))
}else{
	summary.dat=read.table(paste0(output_dir,"pre_meta_analysis_qc_summary.txt"),stringsAsFactors=FALSE,header=TRUE,comment.char="",sep="\t")
}
summary.dat[[sumstat_name]]=NA

#Read each set of summary statistics into R for processing. Record the number of SNPs left at each stage of QC in a table
ss.dat=read.table(ss_file, stringsAsFactors=FALSE, quote="", comment.char="", header=TRUE, sep=sep_type)
summary.dat[summary.dat$Filter=="Initial # SNPs",sumstat_name]=nrow(ss.dat)

#Filter EAF <0.05 and >0.95
ss.dat=ss.dat[ss.dat[,EAF_col]>lower_EAF_bound & ss.dat[,EAF_col]<upper_EAF_bound,]
summary.dat[summary.dat$Filter==paste0("EAF <",upper_EAF_bound,", >",lower_EAF_bound),sumstat_name]=nrow(ss.dat)

#Remove missing data (use complete() function)
if(is_n_col){
	ss.dat=ss.dat[complete.cases(ss.dat[c(EAF_col,n_col,OA_col,EA_col,P_col,effect_col,se_col)]),]
} else{
	ss.dat=ss.dat[complete.cases(ss.dat[c(EAF_col,OA_col,EA_col,P_col,effect_col,se_col)]),]
}
summary.dat[summary.dat$Filter=="Complete cases",sumstat_name]=nrow(ss.dat)


#Remove alleles that are not A, C, T, or G
ss.dat=ss.dat[(ss.dat[OA_col]=="A" | ss.dat[OA_col]=="T" | ss.dat[OA_col]=="C" | ss.dat[OA_col]=="G") & (ss.dat[EA_col]=="A" | ss.dat[EA_col]=="T" | ss.dat[EA_col]=="C" | ss.dat[EA_col]=="G"),]
summary.dat[summary.dat$Filter=="Keep only A,C,T,G alleles",sumstat_name]=nrow(ss.dat)


#Remove non-numeric Betas
ss.dat=ss.dat[!(is.na(as.numeric(ss.dat[[effect_col]]))),]
summary.dat[summary.dat$Filter=="Remove non-numeric Betas",sumstat_name]=nrow(ss.dat)

#Remove negative or non-numeric St Err
ss.dat=ss.dat[!(is.na(as.numeric(ss.dat[[se_col]]))),]
ss.dat=ss.dat[as.numeric(ss.dat[[se_col]]) >= 0,]
summary.dat[summary.dat$Filter=="Remove negative, non-numeric StErr",sumstat_name]=nrow(ss.dat)


#Filter on imputation quality, R2>0.3, INFO>0.8 (Use Rsq GSA for CC), and filter out any scores above 1
if(add_imp_data_for_biovu){
	#imp_data=read.table("/fs0/straubp/MEGA_recalled_full_set/post_imputation_qc_rerun_AIMs/snp_filter_r2/all.chr.imputation.info.dedup.txt", stringsAsFactors=FALSE, quote="",comment.char="", header=TRUE)
	#imp_data2=aggregate(Rsq ~ SNP, data=imp_data, mean)
	#write.table(imp_data2, "/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/all.chr.imputation.info.aggregated.mean.txt", row.names=FALSE, quote=FALSE)
	imp_data<-read.table("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/all.chr.imputation.info.aggregated.mean.txt", stringsAsFactors=FALSE, quote="",comment.char="", header=TRUE)
	ss.dat$ChrPosA1A2=paste(ss.dat[[chr_col]],ss.dat[[bp_col]],ss.dat$Allele1,ss.dat$Allele2, sep=":")
	ss.dat=merge(ss.dat,imp_data,by.x="ChrPosA1A2",by.y="SNP",all.x=FALSE,all.y=FALSE)
	ss.dat=ss.dat[!is.na(ss.dat$Rsq),]
	summary.dat[summary.dat$Filter=="# BioVU snps after merging Rsq info",sumstat_name]=nrow(ss.dat)
}

if(is_rsq){
	ss.dat=ss.dat[ss.dat[,imp_column]>0.3,]
	summary.dat[summary.dat$Filter=="Rsq>0.3,INFO>0.8",sumstat_name]=nrow(ss.dat)
} else{
	ss.dat=ss.dat[ss.dat[,imp_column]>0.8,]
	summary.dat[summary.dat$Filter=="Rsq>0.3,INFO>0.8",sumstat_name]=nrow(ss.dat)
}
ss.dat=ss.dat[ss.dat[,imp_column]<=1,]
summary.dat[summary.dat$Filter=="Rsq<1,INFO<1",sumstat_name]=nrow(ss.dat)


#Create ChrPosID column to merge files with (must ensure all sum stats are in the same build for this)
ss.dat$chr_pos=paste(ss.dat[[chr_col]],ss.dat[[bp_col]], sep=":")

#Plot imputation for each summary statistics, make sure for each one there’s a peak at high imputation scores with trailing tail to the left, and that distribution seems similar across sum-stats
pdf(paste0("imputation_histogram_",sumstat_name,".pdf"))
hist(ss.dat[[imp_column]], xlab=imp_column, main=paste0("Imputation values of ",sumstat_name," sum stats"),breaks=100)
dev.off()


#Compare p-val calculated from z-stats (Z=beta/se(beta)) vs p-values calculated by the software and being plotted in the meta analysis
if(!is.numeric(ss.dat[[effect_col]])){
	ss.dat[[effect_col]]=as.numeric(ss.dat[[effect_col]])
}
if(!is.numeric(ss.dat[[se_col]])){
	ss.dat[[se_col]]=as.numeric(ss.dat[[se_col]])
}
if(is_effect_BETA){
	ss.dat$manual.p=2*pnorm(-abs(ss.dat[[effect_col]]/ss.dat[[se_col]]))
}else{
	ss.dat$manual.p=2*pnorm(-abs(log(ss.dat[[effect_col]])/ss.dat[[se_col]]))
}

png(paste0(sumstat_name,"_provided_p_values_vs_manually_calc_p_val.png"), width=2880, height=1440, res=360, type="cairo")
smoothScatter(ss.dat[c("manual.p",P_col)], ylim=c(0,1), xlim=c(0,1), ylab=paste0(sumstat_name," Provided P values"),xlab="Beta-calculated P values",main=paste0(sumstat_name," provided p values vs manually calculated p values"),nrpoints = Inf)
dev.off()



#Plot each GWAS sumstats EAF vs. Reference population EAF (e.g. 1k genomes) before QC
ref_eaf.dat=read.table(ref_eaf_file, stringsAsFactors=FALSE, quote="", comment.char="", header=TRUE)
ref_ss_eaf.dat=merge(ss.dat,ref_eaf.dat[c(ref_chrbp_col,ref_eaf_col,ref_ea_col,ref_oa_col)],by.x="chr_pos",by.y=ref_chrbp_col,all.x=FALSE,all.y=FALSE)
ref_ss_eaf.dat=ref_ss_eaf.dat[ref_ss_eaf.dat[[EA_col]]==ref_ss_eaf.dat[[ref_oa_col]] | ref_ss_eaf.dat[[EA_col]]==ref_ss_eaf.dat[[ref_ea_col]],]

png(paste0("ref_eaf_vs_",sumstat_name,"_eaf_before_qc.png"), width=2880, height=1440, res=360, type="cairo")
smoothScatter(ref_ss_eaf.dat[c(ref_eaf_col,EAF_col)], ylim=c(0,1), xlim=c(0,1), ylab=paste0(sumstat_name," EAF"),xlab="Reference EAF",main=paste0("Reference EAF vs ",sumstat_name," SNPs EAF before QC"),nrpoints = Inf)
dev.off()

names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == ref_eaf_col] <- "ref_eaf"
names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == ref_ea_col] <- "ref_ea"
names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == ref_oa_col] <- "ref_oa"

##calculate the number of SNPs once merged with reference file, and # of SNPs that are aligned vs. misaligned w EA vs. OA
summary.dat[summary.dat$Filter=="# SNPs after merging with Ref EAF file",sumstat_name]=nrow(ref_ss_eaf.dat)
summary.dat[summary.dat$Filter=="Initial # SNPs where SS EA == Ref EA",sumstat_name]=nrow(ref_ss_eaf.dat[ref_ss_eaf.dat[[EA_col]]==ref_ss_eaf.dat$ref_ea,])
summary.dat[summary.dat$Filter=="Initial # SNPs where SS EA == Ref OA",sumstat_name]=nrow(ref_ss_eaf.dat[ref_ss_eaf.dat[[EA_col]]==ref_ss_eaf.dat$ref_oa,])
switch_ea_and_oa_cols=summary.dat[summary.dat$Filter=="Initial # SNPs where SS EA == Ref OA",sumstat_name]

##switch any EA vs OA SNPs that are misaligned, including options for if all snps are misaligned and if no snps are misaligned
if(switch_ea_and_oa_cols==summary.dat[summary.dat$Filter=="# SNPs after merging with Ref EAF file",sumstat_name]){
	names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == EA_col] <- "ss_OA"
	names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == OA_col] <- "ss_EA"
	names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == EAF_col] <- "ss_EAF"
	ref_ss_eaf.dat$ss_EAF=1-ref_ss_eaf.dat$ss_EAF
}else if(switch_ea_and_oa_cols!=0){
	cols_to_switch=ref_ss_eaf.dat[[EA_col]]==ref_ss_eaf.dat$ref_oa
	ref_ss_eaf.dat$ss_EA[cols_to_switch]=ref_ss_eaf.dat[[OA_col]][cols_to_switch]
	ref_ss_eaf.dat$ss_OA[cols_to_switch]=ref_ss_eaf.dat[[EA_col]][cols_to_switch]	
	ref_ss_eaf.dat$ss_EA[!cols_to_switch]=ref_ss_eaf.dat[[EA_col]][!cols_to_switch]
	ref_ss_eaf.dat$ss_OA[!cols_to_switch]=ref_ss_eaf.dat[[OA_col]][!cols_to_switch]
	ref_ss_eaf.dat$ss_EAF[cols_to_switch]=1-ref_ss_eaf.dat[[EAF_col]][cols_to_switch]
	ref_ss_eaf.dat$ss_EAF[!cols_to_switch]=ref_ss_eaf.dat[[EAF_col]][!cols_to_switch]
	ref_ss_eaf.dat=ref_ss_eaf.dat[,!(names(ref_ss_eaf.dat) %in% c(EAF_col,OA_col,EA_col))]
}else{
	names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == EA_col] <- "ss_EA"
	names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == OA_col] <- "ss_OA"
	names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == EAF_col] <- "ss_EAF"
}
#calculate final number of switched SNPs and final number of properly aligned SNPs (should be all of them)
summary.dat[summary.dat$Filter=="Switched EA and OA cols",sumstat_name]=switch_ea_and_oa_cols
summary.dat[summary.dat$Filter=="Final # SNPs where SS EA == Ref EA",sumstat_name]=nrow(ref_ss_eaf.dat[ref_ss_eaf.dat$ss_EA==ref_ss_eaf.dat$ref_ea,])

#if EAF is still actually MAF, then fix that
is_ss_eaf_actually_maf=max(ref_ss_eaf.dat$ss_EAF)<=0.5
if(is_ss_eaf_actually_maf){
	#create column for minor allele in reference data
	ref_ea_minor <- ref_ss_eaf.dat$ref_eaf<=0.5
	ref_ss_eaf.dat$ref_minor_allele=NA
	ref_ss_eaf.dat[ref_ea_minor, "ref_minor_allele"]=ref_ss_eaf.dat$ref_ea[ref_ea_minor]
	ref_ss_eaf.dat[!ref_ea_minor, "ref_minor_allele"]=ref_ss_eaf.dat$ref_oa[!ref_ea_minor]

	##create adjusted ss_eaf column
	#if ss EA==ref MA, no change. if ss EA != ref MAF, then compute 1-ss AF
	ss_aligned <- ref_ss_eaf.dat$EA==ref_ss_eaf.dat$ref_minor_allele
	ref_ss_eaf.dat$ss_eaf_aligned=NA
	ref_ss_eaf.dat[ss_aligned, "ss_eaf_aligned"]=ref_ss_eaf.dat$EAF[ss_aligned]
	ref_ss_eaf.dat[!ss_aligned, "ss_eaf_aligned"]=1-ref_ss_eaf.dat$EAF[!ss_aligned]

	#remove any SNPs that fall outside of 10% of the reference AF
	ref_ss_eaf.dat=subset(ref_ss_eaf.dat, select=-(ss_EAF))
	names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == ss_eaf_aligned] <- "ss_EAF"
}

#remove any SNPs that are outside of 10% of reference MAF (must convert EAF to MAF for this)
ss_is_eaf_maf = ref_ss_eaf.dat$ss_EAF<=0.5
ref_ss_eaf.dat$ss_eaf_maf[ss_is_eaf_maf]=ref_ss_eaf.dat$ss_EAF[ss_is_eaf_maf]
ref_ss_eaf.dat$ss_eaf_maf[!ss_is_eaf_maf]=1-ref_ss_eaf.dat$ss_EAF[!ss_is_eaf_maf]
ref_ss_eaf.dat$ref_eaf_maf[ss_is_eaf_maf]=ref_ss_eaf.dat$ref_eaf[ss_is_eaf_maf]
ref_ss_eaf.dat$ref_eaf_maf[!ss_is_eaf_maf]=1-ref_ss_eaf.dat$ref_eaf[!ss_is_eaf_maf]
ref_ss_eaf.dat=ref_ss_eaf.dat[0.9*ref_ss_eaf.dat$ref_eaf_maf<=ref_ss_eaf.dat$ss_eaf_maf & ref_ss_eaf.dat$ss_eaf_maf<=1.1*ref_ss_eaf.dat$ref_eaf_maf,]
summary.dat[summary.dat$Filter=="Exclude AEF outside 10% of ref AF",sumstat_name]=nrow(ref_ss_eaf.dat)

png(paste0("ref_eaf_vs_",sumstat_name,"_eaf_after_qc.png"), width=2880, height=1440, res=360, type="cairo")
smoothScatter(ref_ss_eaf.dat[c("ref_eaf","ss_EAF")], ylim=c(0,1), xlim=c(0,1), ylab=paste0(sumstat_name," EAF"),xlab="Reference EAF",main=paste0("Reference EAF vs ",sumstat_name," SNPs EAF after QC"),nrpoints = Inf)
dev.off()


#Calculate median standard error and max sample size and produce a plot of c(proportionality constant)/median(SE) vs. Sqrt(max(N))
#	Equation to calculate C ~ median(1/(sqrt(2*MAFj(1-MAFj))
ref_ss_eaf.dat$C_const=1/sqrt(2*ref_ss_eaf.dat$ss_EAF*(1-ref_ss_eaf.dat$ss_EAF)) ##change EAF column to MAF
median_c_const=median(ref_ss_eaf.dat$C_const)
median_se=median(ref_ss_eaf.dat[[se_col]])
c_median_se=median_c_const/median_se
if(is_n_col){
	sqrt_max_n=sqrt(max(ref_ss_eaf.dat[[n_col]]))
}else{
	sqrt_max_n=sqrt(sample_size)
}

if(is_first_sumstat_analysed){
	c_constant.dat=data.frame("Sumstat.name"=sumstat_name,"C.median.se"=round(c_median_se,4),"sqrt.max.N"=round(sqrt_max_n,4))
}else{
	c_constant.dat=read.table(paste0(output_dir,"c_div_by_median_se_vs_sqrt_n_across_biobanks.txt"), stringsAsFactors=FALSE, quote="",comment.char="",header=TRUE)
	c_constant.dat=rbind(c_constant.dat,c(sumstat_name,round(c_median_se,4),round(sqrt_max_n,4)))
}
write.table(c_constant.dat,paste0(output_dir,"c_div_by_median_se_vs_sqrt_n_across_biobanks.txt"),row.names=FALSE,quote=FALSE,sep="\t")
png(paste0(output_dir,"c_div_by_median_se_vs_sqrt_n_across_biobanks.png"), width=2880, height=1440, res=360, type="cairo")
ggplot(c_constant.dat, aes(x=C.median.se, y=sqrt.max.N)) +
	ylim(0,max(max(as.numeric(c_constant.dat$C.median.se)),max(as.numeric(c_constant.dat$sqrt.max.N)))) +
	xlim(0,max(max(as.numeric(c_constant.dat$C.median.se)),max(as.numeric(c_constant.dat$sqrt.max.N)))) +
	geom_abline(intercept = 0,slope=1) +
	labs(title="C/median(SE) vs. sqrt(max(N)) across biobanks", x="C/median(SE)", y="sqrt(max(N))") +
	geom_point() +
	geom_text(
		label=c_constant.dat$Sumstat.name,
		nudge_x=0.1, nudge_y=0.1,
		check_overlap=T
	)
dev.off()



#Reformat summary statistics so Effect allele is the first allele column, and is labeled EA, freq is labeled EAF, and then non-effect allele is labeled OA and is second column
if(add_imp_data_for_biovu){
	ref_ss_eaf.dat[c("ref_eaf","ref_ea","ref_oa","manual.p","ChrPosA1A2","C_const")]=NULL
}else{
	ref_ss_eaf.dat[c(E"ref_eaf","ref_ea","ref_oa","manual.p","C_const")]=NULL
	}
names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == "ss_EA"] <- "EA"
names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == "ss_OA"] <- "OA"
names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == "ss_EAF"] <- "EAF"
names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == se_col] <- "SE"
names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == P_col] <- "P"
if(is_n_col){
	names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == n_col] <- "N"
}
if(is_effect_BETA){
	names(ref_ss_eaf.dat)[names(ref_ss_eaf.dat) == effect_col] <- "BETA"
}




write.table(ref_ss_eaf.dat, paste0(sumstat_name,"_sumstats_for_meta_analysis_cleaned.txt"),quote=FALSE, sep="\t",row.names=FALSE)
if(is_x_present){
	ss_auto.dat=ref_ss_eaf.dat[ref_ss_eaf.dat[[chr_col]]!=x_col,]
	ss_x.dat=ref_ss_eaf.dat[ref_ss_eaf.dat[[chr_col]]==x_col,]
	summary.dat[summary.dat$Filter=="# autosomal SNPs",sumstat_name]=nrow(ss_auto.dat)
	summary.dat[summary.dat$Filter=="# xchr SNPs",sumstat_name]=nrow(ss_x.dat)
	write.table(ss_auto.dat, paste0(sumstat_name,"_autosome_sumstats_for_meta_analysis_cleaned.txt"),quote=FALSE, sep="\t",row.names=FALSE)
	write.table(ss_x.dat, paste0(sumstat_name,"_xchr_sumstats_for_meta_analysis_cleaned.txt"),quote=FALSE, sep="\t",row.names=FALSE)
}else{
	write.table(ref_ss_eaf.dat, paste0(sumstat_name,"_autosome_sumstats_for_meta_analysis_cleaned.txt"),quote=FALSE, sep="\t",row.names=FALSE)
	summary.dat[summary.dat$Filter=="# autosomal SNPs",sumstat_name]=nrow(ref_ss_eaf.dat)
	summary.dat[summary.dat$Filter=="# xchr SNPs",sumstat_name]=0
}

write.table(summary.dat,paste0(output_dir,"pre_meta_analysis_qc_summary.txt"),quote=FALSE,row.names=FALSE,sep="\t")

