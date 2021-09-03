Metal_Meta_Analysis_Pipeline

Example can be found in the folder:

    /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss/

First run each set of summary statistics through the pre_meta_analysis_qc.R script in R:

Adapted from Winkler et. Al meta analysis Easy QC steps by Slavi Goleva and Lea Davis

##To be run with each individual summary statistic file##
#comment out variables at top of script to reflect summary statistics or put separately into R and then source() the script#

#typical run time per summary stat: ~15min

Pre-meta-analysis QC (steps 1-16) is done in:

    https://github.com/sbgoleva/Metal_Meta_Analysis_Pipeline/blob/main/pre_meta_analysis_qc.R
    
For this script, just change the parameters commented out at the top for each set of summary statistics and run independently. 

    #To Be Done with each individual summary statistic file##
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
    chr_col="CHR"
    bp_col="POS"
    rs_id_col="SNPID"


    lower_EAF_bound=0.05 #desired EAF filters
    upper_EAF_bound=0.95

    ##EAF reference file info. Below is using 1kg snp list in European Ancestry
    ref_eaf_file="/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Allele_Freqs/1000GP_p3v5_legends_rbind.noDup.noMono.noCnv.noCnAll.afref.EUR.txt.gz"
    ref_eaf_col="eaf"
    ref_chrbp_col="cptid"
    ref_ea_col="ea"
    ref_oa_col="oa"
    ref_x_name="X"

    is_first_sumstat_analysed=TRUE

    is_x_present=FALSE
    x_col="" ##will be written out as '23' to standardize across biobanks

    is_ss_reanalysis=TRUE #TRUE if you have already analyzed these sumstats before


running pre_meta_analysis QC script in R:

    source("/data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/pre_meta_analysis_qc.R")
    
To check after each sumstats has run:

    1. 1/Median SE vs Sqrt(N) plot runs along identity line, and if not this can be explained by differences in GWAS analysis type, etc
    2. pre_meta_analysis_qc_summary.txt percents lost and # of SNP look ok
        A. If BioVU file, "# BioVU snps after merging Rsq info" step should be 100%
        B. If any step loses more than ~15%, make sure there is a good reason/it makes sense
        C. Most SNPs remain after "# SNPs after merging with Ref EAF file" step
        D. "Initial # SNPs where SS EA == Ref EA" and "Initial # SNPs where SS EA == Ref OA" should add up to "# SNPs after merging with Ref EAF file" step
        E. "Initial # SNPs where SS EA == Ref OA" step should match "Switched EA and OA cols" step
        F. "Final # SNPs where SS EA == Ref EA" should be 100%
        G. "Exclude AEF outside 10% of ref AF" should not exceed ~15%. If it does, consider using more appropriate reference population file.
        H. "#autosomal SNPs" and "# xchr SNPs" should add up to "Exclude AEF outside 10% of ref AF"
    3. Check the following in the directory for each summary statistics file:
        A. ref_eaf_vs_[biobank_name]_eaf_after_qc.png plot should be along identity line, with no outliers present
        B. imputation_histogram_[biobank_name].pdf should not contain values below 0.3 if Rsq and 0.8 if INFO. 
            Most scores should be clusered around 1. 
            No values should exceed 1. 
        C. [biobank_name]_provided_p_values_vs_manually_calc_p_val.png should be on identity line. There should be no outliers. 
    4.  Calculate lambda GC for each study file – can do this in LDSC (make sure it doesn’t exceed 1.1)



If everything looks good, proceed to meta-analysis using outputted summary statistics

Edit and run metal.sh script screen, this takes around 15 minutes:


If you have used the pre-QC script, then the only thing that will change are the output file names and the summary statistic file names. 
Columns are all standardized as part of the pre-QC script output
To change in run_metal.sh
   
    1. script name:
        /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/METAL_BioVU_UKBBR56_CC/metal_script.sh
    2. log file output:
        /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/metal_fs_gwas_try.log
    3. sum-stat names. Each starts with "PROCESS." The number of lines should equal the number of sum-stats to be meta-analyzed
    4. outfile name ("outfile/prefix .tbl") ***make sure a space separates these two***
        /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/metal_fs_gwas_try .tbl
    5. If all summary statistics don’t use the same covariates, then this will in effect act as different transformations to the effect estimates, 
       so use the sample size scheme in Metal to meta-analyze.
       Else, you can use StErr scheme
    6. ***Do not include x-chrom for heritability analysis or in meta-analysis***
        The x-chr must be meta-analyzed separately and can be concatonated with the rest of the results for follow-up analysis.

To run metal script in Unix screen

    metal /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/METAL_BioVU_UKBBR56_CC/metal_script.sh > /data/davis_lab/golevasb/Projects/PNES_PRS/GWAS/Meta_Analysis_cleaned_ss_absolute_filter_ref_maf/metal_fs_gwas_try.log

Post meta-analysis QC steps:

    1. Check log file for duplicate SNP warnings. For each one, manually figure out which alleles to use based on their frequencies
        e.g. from log file: "WARNING: Bad alleles for marker 'rs16856772', expecting 'a/g' found 'a/t'"
 
Post meta-analysis sumstats filtering and Manhattan/QQ Plots:

    post_meta_plotting.R
    
    
 
    Here you will read summary statistics in to R
    Using the qqman package, this script plots Manhattan and QQ Plots of data
         ***may need to read in individual sumstats to merge CHR and BP info.***
         ***working on getting new metal package installed on accre which can do this automatically***
    This also filters SNPs from only 1 biobank and writes out summary statistics to be used in all future post-meta analysis

Post meta-analysis analyses:

    1. Using running_ldsc.sh, calculate meta-analysis heritability
    2. Using running_ldsc.sh, calculate genetic correlation between related phenotypes
    3. Use FUMA to load sumstats for follow up analyses: https://fuma.ctglab.nl/snp2gene
    4. Using MultiXcan, conduct TWAS (contact Kritika for script)
