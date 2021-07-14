# Metal_Meta_Analysis_Pipeline

First run each set of summary statistics through the pre_meta_analysis_qc.R script in R:

Adapted from Winkler et. Al meta analysis Easy QC steps by Slavi Goleva and Lea Davis

##To Be Done with each individual summary statistic file##

#typical run time per summary stat: ~1hr

Pre-meta-analysis QC (steps 1-16) is done in:

    https://github.com/sbgoleva/Metal_Meta_Analysis_Pipeline/blob/main/pre_meta_analysis_qc.R
    
For this script, just change the parameters commented out at the top for each set of summary statistics and run independently. 

    #file info to fill out#

    ss_file="/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/PNES_GWAS_Saige/20210412_PNES_GWAS_Redo_Matched_EUR_step2_chr.all.txt"
    sumstat_name="BioVU" ##make sure name is unique for each set of sumstats as the output files are based on this
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
    is_x_present=TRUE
    x_col="22"
    #only need chr_col and bp_col if add_imp_data_for_biovu==TRUE
    chr_col="CHR"
    bp_col="POS"


    lower_EAF_bound=0.05 #desired EAF filters
    upper_EAF_bound=0.95

    ##EAF reference file info. Below is using 1kg snp list in European Ancestry
    ref_eaf_file="/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Allele_Freqs/1000GP_p3v5_legends_rbind.noDup.noMono.noCnv.noCnAll.afref.EUR.txt.gz"
    ref_eaf_col="eaf"
    ref_chrbp_col="cptid"

    is_first_sumstat_analysed=TRUE
    ####################################################################

running pre_meta_analysis QC script in R:

    source("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/pre_meta_analysis_qc.R")
    
Then, make sure all the plots produced for each sumstat look ok.

If so, proceed to meta-analysis using outputted summary statistics


Run meta analysis in screen, this only takes around 15 minutes:

##For each set of summary statistics, put this information in individually. 

#If you have used the pre-QC script, 

#then the only thing that will change are the sample size column presence column

#and the file name

    module load GCCcore/.8.2.0
    module load Metal

    metal

    GENOMICCONTROL OFF

    AVERAGEFREQ ON
    MINMAXFREQ ON

    SCHEME SAMPLESIZE
    
    SEPARATOR TAB
    MARKER chr_pos
    ALLELE EA OA
    PVALUELABEL P
    EFFECTLABEL BETA #change to log(OR) if OR instead of BETA
    WEIGHTLABEL N #if no N column, remove and uncomment out two lines beneath this
    #WEIGHTLABEL DONTUSECOLUMN
    #DEFAULTWEIGHT 10000
    STDERR SE
    FREQLABEL EAF

    PROCESS paste0(output_dir,sumstat_name,"_autosome_sumstats_for_meta_analysis_cleaned.txt") #repeat this line and just keep changing ss name


After meta-analysis QC steps:

    17.	Calculate lambda GC for each study file – can do this in LDSC (make sure it doesn’t exceed 1.1)
    18.	Inverse variance weighted meta analysis using a fixed effects model with METAL, as follows
    19.	If you have duplicate SNPs, flag and figure out which alleles to use based on their frequencies. Print to another file and go through them together
    20.	Print min and max difference in allele freqs
    21.	If all summary statistics don’t use the same covariates, then this will in effect act as different transformations to the effect estimates, so use the sample size scheme in Metal to meta-analyze
    22.	Do not include x-chrom for heritability analysis
    
 
Post-meta-analysis Analyses:

    Read summary statistics in to R, and using qqman package, plot Manhattan and QQ Plots of data
    Using LDSC, calculate heritability
    Using LDSC, calculate genetic correlation between related phenotypes
    Using MultiXcan, conduct TWAS
