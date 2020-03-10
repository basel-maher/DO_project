# DO_project pipeline

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




9. Map QTL

    * [src/map_qtl.R](../src/map_qtl.R)
    * will also perform permutation analysis
        * output: 
		* results/Rdata/qtl_perms/
		
		* Ellipticity calculation and mapping in [src/ellipticity.R](../src/ellipticity.R)
        

10. eQTL mapping
    * [put pipeline here](rna_seq_pipeline.md)
    * PEER stuff
        -used VST + qnormed counts
        - 48 PEER factors, no covars and no intercept


    * Make the annotation file that links between Stringtie gene IDs and gene names:
      * [src/make_annot.R](../src/make_annot.R)
        * requires gene abundance files (Stringtie output)
        * In: results/flat/RNA-seq/abund/
        * Output:
          * results/flat/annot_file.csv
    
    * after RNA-seq preprocessing, normalize the counts
    	* [src/normalize_RNAseq.R](../src/normalize_RNAseq.R)

    * make an eqtl specific cross file
        * [src/make_cross_eqtl.R](../src/make_cross_eqtl.R)

    - calc eqtl perms
    	* based on mapping local varying covariates, we chose to map local
            eqtl with sex and all 48 PEER covars.
        	* output in ./results/Rdata/eqtl_perms/  
                * use maximum LODs at alpha=0.05 for mapping:
                    * Autosomal: 10.89
                    * chrX: 11.55

    * Map eQTL

    	* [src/map_local_and_distal_eqtl.R](../src/map_local_and_distal_eqtl.R)
    	* outputs: 
        	* ./results/flat/local_eqtl_peaks.csv 
        	* ./results/Rdata/local_eqtl.Rdata
			* CHECK OUTS, SOME FROM RIVANNA (1-4)
        
        

11. Networks
    
    * Make bone "superset" [src/make_bone_geneset.R](../src/make_bone_geneset.R)
		* used code from https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/
        * Get gene ontology identifiers from AmiGO2.
            * Use "bone" and "osteo* and "ossif*" as terms. stored in ./data/GO_term_bone.txt
        * Also from MGI
        	* human and mouse "osteoporosis" and "bone mineral density" and "osteoblast" and "osteoclast" and "osteocyte"
  
  
    * Full WGCNA networks constructed in [src/WGCNA_working.R](../src/WGCNA_working.R)
    	* Gene ontology for WGCNA network modules also performed here
    
    * outputs:
        * results/Rdata/networks/edata_full.Rdata
        * results/Rdata/networks/wgcna_4.RDS
        * results/Rdata/networks/moduleTraitPvalue_full_4.RData
        * results/Rdata/networks/moduleTraitCor_full_4.RData
        * results/Rdata/networks/geneModMemAnnot_power4.RData
        * results/Rdata/networks/GO_sft4.RData
        
    * Sex-specific WGCNA networks constructed in [src/WGCNA_sex_specific.R](../src/WGCNA_sex_specific.R)
    * Gene ontology for sex-specific WGCNA network modules also performed here
    
     * outputs:
        * results/Rdata/networks/edata_m.Rdata
        * results/Rdata/networks/edata_f.Rdata
        * results/Rdata/networks/wgcna_m_5.RDS
        * results/Rdata/networks/wgcna_f_4.RDS
        * results/Rdata/networks/geneModMemAnnot_m_power5.RData
        * results/Rdata/networks/geneModMemAnnot_f_power4.RData
        * results/Rdata/networks/moduleTraitPvalue_f.RData
        * results/Rdata/networks/moduleTraitPvalue_m.RData
        * results/Rdata/networks/moduleTraitCor_f.RData
        * results/Rdata/networks/moduleTraitCor_m.RData
        * results/Rdata/networks/GO_Females_sft4.RData
		* results/Rdata/networks/GO_Males_sft5.RData


	* Bayesian networks learned with [src/learn_bn.R](../src/learn_bn.R)
    
      * This was run on Rivanna supercomputing cluster
      * output in ./results/Rdata/networks/bn_4/

    * Key Driver Analysis performed in [src/KDA_working.R](../src/KDA_working.R)
      * output:
      	* results/flat/key_driver_analysis.csv
    
    - munge to add coloc data
    

        


