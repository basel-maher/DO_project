# DO project pipeline

1. Munge phenotypic files into one nice file.
	* use [src/munge_raw_pheno_data.R](../src/munge_raw_pheno_data.R)
    * This takes pheno files from data/pheno_data
	* outputs:
		 * results/flat/full_pheno_table.csv
         	 * results/Rdata/full_pheno_table.Rdata

2. Fix GeneSeek files
	* There was a sample confusion. Re-genotyped samples have a .1 appended to their name, except 371.
    	* In ./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20190228_FinalReport, change Sample ID 371 to 371.1
        * This is done in [src/fix_geneseek.R](../src/fix_geneseek.R)
        	* Add header back manually. See R file for instructions
    * Fix/create ./data/GIGAMUGA/merged/Merged_Sample_Map.txt manually. Change the second 371 entry to 371.1
    	* create Merged_Sample_Map.txt by merging all sample maps from geneseek experiments and editing indices appropriately. Done in TextEdit.
        * Generate data/GIGAMUGA/merged/Merged_FinalReport.txt. 
        	* This is done in [src/fix_geneseek.R](../src/fix_geneseek.R)

3. Convert Geneseek FinalReports files to a format r/qtl2 understands. This encodes the DO genotypes.

    * From https://kbroman.org/qtl2/pages/prep_do_data.html
    * use [src/geneseek2qtl2_mod.R](../src/geneseek2qtl2_mod.R)
    * output is in results/GIGAMUGA/qtl_batches1-4

4. Remove bad samples and markers using Argyle.

    * use [src/GIGAMUGA_QC.Rmd](../src/GIGAMUGA_QC.Rmd)
    * output:
		  * results/GIGAMUGA/geno.final_merged.RDS


5. Create cross file for QTL mapping. 

    * [src/make_crossfile.R](../src/make_crossfile.R)

6. calculate genotype and allele probs, as well as kinship matrices 

	* [src/calc_probs.R](../src/calc_probs.R)

7. more QC - Broman et al.

    * script adapted from https://github.com/kbroman/Paper_MPPdiag/blob/master/R/diagnostics.Rmd
    * This is largely buggy, but everything looked fine. Led to the removal of several hundred markers (479)that had greater than 5% genotyping errors.
    * output is results/Rdata/cross_basic_cleaned.Rdata
    	* 109,427 markers remaining

8. recalculate genotype and allele probs using the new cross_basic_cleaned

    * [src/calc_probs.R](../src/calc_probs.R)




9. QTL mapping / analysis

    * [src/map_qtl.R](../src/map_qtl.R)
    * can also get qtl_loc from this script

		* GWAS-QTL overlap
		  * [src/GWAS_qtl_overlap.R](../src/GWAS_qtl_overlap.R)
        
    
   
10. eQTL mapping
    * [doc/rna_seq_pipeline.md](./rna_seq_pipeline.md)

    * Make the annotation file that links between Stringtie gene IDs and gene names:
      * [src/make_annot.R](../src/make_annot.R)
        * requires gene abundance files (Stringtie output)
        * In: results/flat/RNA-seq/abund/
        * Output:
          * results/flat/annot_file.csv
    
    * after RNA-seq preprocessing, normalize the counts
    	* [src/normalize_RNAseq.R](../src/normalize_RNAseq.R)

    * make an eqtl specific cross file
        * [src/make_cross_eQTL.R](../src/make_cross_eQTL.R)

    * eQTL permutations
      * [src/calc_eqtl_perms.R](../src/calc_eqtl_perms.R)
    	* based on mapping local varying covariates, we chose to map local
            eqtl with sex and all 48 PEER covars.
        	* output in ./results/Rdata/eqtl_perms/  
                * use maximum LODs at alpha=0.05 for mapping:
                    * Autosomal: 10.89
                    * chrX: 11.55

    * Map eQTL

    	* [src/map_local_eqtl.R](../src/map_local_eqtl.R)
    	* outputs: 
        	* ./results/Rdata/local_eqtl.Rdata
        
 
11. Merge Analysis 
  * These two create the merge analysis objects:
    * [src/merge_analysis_QTL.R](../src/merge_analysis_QTL.R)
    * [src/merge_analysis_eQTL.R](../src/merge_analysis_eQTL.R)
    
  * This colocalizes merge analysis qtl/eqtl and identifies common missense variants:
  
    * [src/merge_analysis_misc.R](../src/merge_analysis_misc.R)
 
  * This performs actual analyses on CHR1 using the merge objects:
  
    * [src/merge_analysis_chr1.R](../src/merge_analysis_chr1.R)


12. QSOX1 analysis
  * [src/qsox1_analysis.R](../src/qsox1_analysis.R)
  
  
13. IMPC Phenstat analysis
  * This analysis was performed by Samuel Haddox (sh7qe@virginia.edu). The R script "IMPCgetDataTable.R" includes the commands to install and load the necessary packages and an example of how to find datasets with the "IMPC_DXA_004_001" parameter. The "IMPCgetDataTable.R" file also contains the getIMPCtable functions used to build a single master table, "DEXATable"", that contains all the necessary functions to retrieve the data for analysis with Phenstats. "DEXATable" is output as a comma separated values file "IMPC_DEXA_Table.csv", which is then input for the "IMPCPhenStats.R" script.
  The “IMPCPhenStats.R” script iterates through the DEXATable line by line retrieving the
dataset for each knockout mouse strain. PhenStats mixed model regression analysis is then
carried out on datasets with enough samples and the appropriate controls. Results that could
not be obtained either return “Sample N” error, if there are not enough samples for analysis, or
“Not 2 Genotypes” error, if there are not 2 genotypes to test. Results or Errors are appended to
the variable “ResultsTable”, which is also written to a csv file every 10 iterations to prevent loss
of data when run on a personal computer. The “FinalResultsTable” variable is created by
filtering the “ResultsTable” of lines with errors. The “FinalResultsTable” is then ordered by P
value and significant results are subset into three seperate tables based on sexual dimorphism,
“Not_Dimorphic”, “Female_Dimorphs”, and “Male_Dimorphs”.

    * [src/IMPCgetDataTable.R](../src/IMPCgetDataTable.R)
    
    * [src/IMPCPhenStats.R](../src/IMPCPhenStats.R)



14. Networks
    
    * Make bone "superset" [src/make_bone_geneset.R](../src/make_bone_geneset.R)
		* used code from https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
        * Get gene ontology identifiers from AmiGO2.
            * Use "bone" and "osteo* and "ossif*" as terms. stored in ./data/GO_term_bone.txt
        * Also from MGI
        	* human and mouse "osteoporosis" and "bone mineral density" and "osteoblast" and "osteoclast" and "osteocyte"
        
  
  
    * Full WGCNA networks constructed in [src/WGCNA_working.R](../src/WGCNA_working.R)
    	* Gene ontology for WGCNA network modules also performed here
    

        
    * Sex-specific WGCNA networks constructed in [src/WGCNA_sex_specific.R](../src/WGCNA_sex_specific.R)
    * Gene ontology for sex-specific WGCNA network modules also performed here
    

	* Bayesian networks learned with [src/learn_bn.R](../src/learn_bn.R)
      * This was run on Rivanna supercomputing cluster

    * Key Driver Analysis performed in [src/KDA_working.R](../src/KDA_working.R) and [src/KDA_sex_based.R](../src/KDA_sex_based.R)
   
    * Miscellaneous KDA analyses, such as converting to human homologs and Gene Ontology analyses in [src/KDA_analyses.R](../src/KDA_analyses.R)
 
   * Annotate with coloc results from GTEx
      * [src/annotate_KDA_GTEx.R](../src/annotate_KDA_GTEx.R)
    

15. Single Cell RNA-seq
    * Key Driver Analysis performed in [src/seurat_analysis.R](../src/seurat_analysis.R)
      
16. IMPC data for Glt8d2

    * [src/impc_glt8d2.R](../src/impc_glt8d2.R)

#mouse_human_syntenic regions supplementary fig
#plots for paper
#coloc

