############################################################################
#
#		This file defines the QC function for GWAS, Reference SNP alignment, 
#		and LD mismatch function that can be used in the multiple-ancestry
#		fine-mapping
#
############################################################################
##############################################
#	Step 1. QC of GWAS for each ancestry
###############################################

#' GWAS Quality Control Function
#'
#' This function performs several quality control steps on GWAS data:
#' 1. Removes strand ambiguous variants.
#' 2. Removes multi-allelic variants.
#' 3. Excludes the MHC complex region.
#' 4. Filters SNPs based on Minor Allele Frequency (MAF).
#'
#' @param dat A data frame containing the GWAS data. Columns "CHR", "POS", "CHR_POS", 
#' "REF", "ALT", and (optionally) "MAF" are required.
#' @param MAF_threshold The threshold for Minor Allele Frequency filtering. SNPs 
#' with MAF below this threshold or above (1 - threshold) will be removed.
#'
#' @return A filtered data frame after performing the quality control steps.
#' @import dplyr
#' @export

GWAS_QC <- function(dat, MAF_threshold, rm_strand = TRUE, rm_multi = TRUE, rm_indel = TRUE, rm_mhc = TRUE,rm_maf = TRUE) { 
  # Column name to uppers
  colnames(dat) <- toupper(colnames(dat))
  # Check if required columns are present and to upper
  if(!all(c("CHR", "POS", "REF", "ALT") %in% colnames(dat))) {
    stop("Required columns 'CHR', 'POS', 'REF', and 'ALT' are missing.")
  }
  
  dat <- dat %>% 
    mutate(CHR_POS = paste0(CHR, "_", POS))
  # CHR, POS, REF, ALT, and MAF are required for QC
  # QC Step 1: Remove Strand Ambiguous Variants
  # These are SNPs that can't be differentiated based on the strand and may introduce errors.
  if(rm_strand){
    dat <- dat %>%
      filter(!(
        (REF == "G" & ALT == "C") |
          (REF == "C" & ALT == "G") |
          (REF == "T" & ALT == "A") |
          (REF == "A" & ALT == "T")
      ))
  }
  # QC Step 2: Remove Multi-allelic Variants
  if(rm_multi){
    multi_allelic <- dat$CHR_POS[duplicated(dat$CHR_POS)]
    if(length(multi_allelic)>0){
      dat<-dat %>%filter(!(CHR_POS%in%multi_allelic))
    }
  }
  
  # QC Step 3: Remove Indel if necessary
  if(rm_indel){
    dat <- dat %>% 
      filter(!(nchar(REF) + nchar(ALT) > 2))	
  }
  
  # QC Step 4: Exclude MHC Complex Region
  # This region can be particularly complex and may need special handling.
  if(rm_mhc){
    dat <- dat %>% 
      filter(!(CHR == 6 & POS > 25e6 & POS < 34e6))
  }
  
  
  # QC Step 5: Minor Allele Frequency (MAF) Filter
  # Remove SNPs with MAF below the threshold or above (1 - threshold).
  if(rm_maf){
    if("MAF"%in%colnames(dat)){
      dat <- dat %>% 
        filter(MAF > MAF_threshold, MAF < 1 - MAF_threshold)
    }
  }
  return(dat)
}
##############################################
#	Step 2. Find common set of SNPs, and align 
#	the reference allele, switch the sign of beta
#	or z score based on the allele
###############################################

#' Identify Common SNPs Between Two Datasets
#'
#' This function identifies common SNPs (based on "CHR_POS") between a summary statistics 
#' dataset (`sumstat`) and a reference dataset (`ref`). It then ensures that the 
#' alleles match between the two datasets, either in the original or flipped orientation.
#'
#' @param sumstat A data frame containing the summary statistics data with at least 
#' the columns "CHR_POS", "REF", and "ALT".
#' @param ref A data frame serving as the reference, with at least the columns 
#' "CHR_POS", "REF", and "ALT".
#'
#' @return A vector of "CHR_POS" values corresponding to the SNPs common between 
#' `sumstat` and `ref` with matching allele information.
#' @import dplyr
#' @export

find_common_snps <- function(sumstat, ref) {
  
  # Identifying common SNPs between sumstat and ref datasets
  common_SNP <- intersect(sumstat$CHR_POS, ref$CHR_POS)
  
  # Subsetting data frames based on common SNPs and arranging by order of common_SNP
  sumstat <- sumstat %>% 
    filter(CHR_POS %in% common_SNP) %>% 
    arrange(match(CHR_POS, common_SNP))
  
  ref <- ref %>% 
    filter(CHR_POS %in% common_SNP) %>% 
    arrange(match(CHR_POS, common_SNP))
  
  matched_pos<-which((sumstat$REF==ref$REF&sumstat$ALT==ref$ALT)|(sumstat$REF==ref$ALT&sumstat$ALT==ref$REF))
  
  common_SNP<-sumstat[matched_pos,]%>%pull(CHR_POS)
  
  return(common_SNP)
}
# Match the reference allele with the sumstat used, 
# flip the beta and zscore if the reference allele is flipped.

#' Allele Flip Function for SNP Matching
#'
#' This function modifies the reference dataset (`ref`) based on the summary 
#' statistics dataset (`sumstat`) to ensure that the alleles match. If the alleles 
#' do not match, the function will flip the alleles and corresponding statistics 
#' in the reference dataset.
#'
#' @param sumstat A data frame containing the summary statistics data with at least 
#' the columns "REF".
#' @param ref A data frame serving as the reference with at least the columns 
#' "REF", "ALT", "BETA", and "SE".
#'
#' @return A modified version of `ref` where alleles have been flipped to match 
#' `sumstat` and corresponding statistics (BETA and Z) have been updated accordingly.
#' @import dplyr
#' @export

allele_flip<-function(sumstat,ref){
  ref<-ref %>% mutate(flip_index = ifelse(ref$REF==sumstat$REF,1,-1))%>%
    mutate(BETA = BETA * flip_index)%>%
    mutate(Z = BETA/SE, REF = sumstat$REF, ALT = sumstat$ALT)%>%
    dplyr::select(-flip_index)
  return(ref)
}
#################################################
#	Step 3. Match up the reference panel with
#	the GWAS reference and alternative
#	alleles
###############################################
# Here we use plink to do the referernce allele flip
# Write out the reference SNP for the genotype data
# plink --bfile ref_geno --ref-allele ref.txt --make-bed --out ref_geno_flip
