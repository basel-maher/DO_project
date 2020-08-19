# GWAS-eQTL colocalization

* Essentially, for each GWAS, we took the lead variants and identified all genes within +/- 1 Mbp of a variant. Then, we converted BAN genes to their human homologs using the MGI homolgy table. For each BAN within +/- 1 Mbp of a GWAS lead variant, we took all GTEx V7 variants within +/- 200 Kbp of the lead variant, and colocalized their eQTL p-values with the GWAS p-values. This was performed for "gwas" and "gwas". Much of this was performed on our high-performance computing servers.



1. Download GTEx V7 all eQTL and uncompress.
  * Download all GTEx V7 SNP-gene associations for all 48 tissues.
    * From https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_all_associations.tar.gz
    * Also download GTEx V7 lookup table to convert RSIDs to chr_pos_ref_alt_build format (https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz)
    * Also download sample counts by tissue summary. This data isn't available anymore for GTEx v7, but can be found here:
    
  * Uncompress using gunzip
  
  * Get tissue list
  
  ```bash
  ls | cut -d . -f1 > tissue_list_v7
  ```
  
  * split into smaller files (100 M) to facilitate grabbing relevant variants, and get the split tissue list.
  
  ```bash
  for f in *.txt; do out=$(echo $f | cut -d. -f1); split -C 100M -d $f ./split/$out; done
  ```

  * get split tissue list
  
  ```bash
  ls > ../../tis_supersplit_v7
  ```
  
2. Convert BANs to homologs, and get homologs within 1 Mbp of GWAS variants.

  * This is performed in [src/coloc_pre.R](../src/coloc_pre.R), and uses the following data:
  
    * MGI human-mouse homology table (http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt) (Mar-14-2020)
    * BANs from the KDA analyses ([src/KDA_working.R](../src/KDA_working.R) and [src/KDA_sex_based.R](../src/KDA_sex_based.R))
    * Morris et al. 2018 eBMD GWAS data (http://www.gefos.org/?q=content/data-release-2018)
    * Estrada et al. 2012 FNBMD and LSBMD GWAS data (http://www.gefos.org/?q=content/data-release-2012)
    * For each GWAS, I also extracted the GWAS lead variants.
  
3. Get all variants within 200 Kbp of genes for analysis.

  * Script for getting the relevant variants from GTEx for eBMD GWAS: [src/get_200k_morris.sh](../src/get_200k_morris.sh). This was run using a SLURM job array. Needs to be run for each of the GTEx tissues (i=1-48)
  
  * Convert RSIDs for the Estrada GWAS and the lead SNPs to GTEx format (chr_pos_ref_alt_build) using the GTEx lookup table (from step 1 above).
    * This was done in R using merge(), to merge by RSID.
    * Output files are named lead_merged, lsbmd_merged and fnbmd_merged.
  
  * Script for getting the relevant variants from GTEx for Estrada FN/LSBMD GWAS: [src/get_200k_estrada.sh](../src/get_200k_estrada.sh)
  
  * For each resulting set of variant files for the two GWAS's, merge into one file per tissue:
  
  ```bash 
  while read line; do
  cat $(find ../eqtl_200k -name *${line}* -print | sort) > ../merged_eqtl_200k/$line
  echo $line
  done < ./tissue_list_v7
  ```

4. Colocalization

  * Colocalize:
    * [src/coloc_morris_v7.R](../src/coloc_morris_v7.R)
    * [src/coloc_estrada_v7.R](../src/coloc_estrada_v7.R)
  
  * Get the results:
    * [src/get_coloc_results.R](../src/get_coloc_results.R)

5. Miscellaneous analyses:

  * [src/coloc_post.R](../src/coloc_post.R)