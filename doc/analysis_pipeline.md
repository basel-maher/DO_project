# DO_project pipeline

1. Munge phenotypic files into one nice file.
  -use ./src/munge_raw_pheno_data.R
    - This takes pheno files from ./data/pheno_data and outputs ./results/flat/full_pheno_table.csv and results/Rdata/full_pheno_table.Rdata

2. Fix GeneSeek files

    - There was a sample confusion. Re-genotyped samples have a .1 appended to their name, except 371.
        - In ./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20190228_FinalReport, change Sample ID 371 to 371.1
            - This is done in ./src/fix_geneseek.R
                -Add header back manually. See R file for instructions
        - Fix/create ./data/GIGAMUGA/merged/Merged_Sample_Map.txt manually. Change the second 371 entry to 371.1
            - create Merged_Sample_Map.txt by merging all sample maps from geneseek experiments and editing indices appropriately. Done in TextEdit.
        - Generate ./data/GIGAMUGA/merged/Merged_FinalReport.txt. 
            - This is done in ./src/fix_geneseek.R

3. Convert Geneseek FinalReports files to a format r/qtl2 understands. This encodes the DO genotypes.

    - From https://kbroman.org/qtl2/pages/prep_do_data.html
    - use ./src/geneseek2qtl2_mod.R
    - output is in results/GIGAMUGA/qtl_batches1-4

4. Remove bad samples and markers using Argyle.

    - use ./src/GIGAMUGA_QC.Rmd
    - output ./results/GIGAMUGA/geno.final_merged.RDS


5. Create cross file for QTL mapping. 

    - ./src/make_crossfile.R

6. calculate genotype and allele probs, as well as kinship matrices 

    - how is kinship calculated?
    - ./src/calc_probs.R

7. more QC - Broman et al.

    - script adapted from https://github.com/kbroman/Paper_MPPdiag/blob/master/R/diagnostics.Rmd
    - cite paper
    - This is largely buggy, but everything looked fine. Led to the removal of several hundred markers (479)that had greater than 5% genotyping errors.
    - output is cross_basic_cleaned (./results/Rdata/cross_basic_cleaned.Rdata)
    - 109,427 markers remaining

8. recalculate genotype and allele probs using the new cross_basic_cleaned

    - ./src/calc_probs.R

9. Calc QTL

    - ./src/map_qtl.R
    - will also perform permutation analysis
        - output in ./results/Rdata/eqtl_perms/   
        

10. Calc eQTL

    - separate instructions for sequencing data preprocessing?
        -probably. go to rna_seq_pipeline.md
    - PEER stuff
        -used VST + qnormed counts
        - 48 PEER factors, no covars and no intercept

    - after RNA-seq preprocessing, normalize the counts
        - ./src/normalize_RNAseq.R

    - make an eqtl specific cross file
        - ./src/make_cross_eqtl.R

    - calc eqtl perms
        -  based on mapping local and distal eqtl with varying covariates, we chose to map local
            eqtl with sex and the first 35 PEER covars, and distal eqtl with all 48 PEER covariates, and not sex.
            - ./src/calc_eqtl_perms.R
                - output in ./results/Rdata/eqtl_perms/   
                - use maximum LODs at 0.05 alpha=0.05 for mapping
                    - local_autosomal: 9.9
                    - local_X: 10.9
                    - distal_autosomal: 10.7
                    - distal_X: 11
                

11. Map eQTL

    - ./src/map_local_and_distal_eqtl.R
    -   outputs: 
        - ./results/flat/local_eqtl_peaks.csv 
        - ./results/flat/distal_eqtl_peaks.csv 
        - ./results/flat/gene_annot_file.csv
        - ./results/Rdata/local_eqtl.Rdata
        - ./results/Rdata/distal_eqtl.Rdata
        
        
-make_annot.R


12. eQTL analysis: trans eQTL hotspots, mediation analysis, NEO structural equation  modeling

13. Networks
    
    - Make bone "superset" ./src/make_bone_geneset.R
        - Get gene ontology identifiers from AmiGO2.
            - Use "bone" and "osteo* and "ossif*" as terms. stored in ./data/
        - also from MGI
         - human and mouse "osteoporosis" and "bone mineral density" and "osteoblast" and "osteoclast" and "osteocyte"
  
  
    - Full WGCNA networks constructed in ./src/WGCNA_working.R
      -Gene ontology for WGCNA network modules also performed here
    
     -outputs:
        - ./results/Rdata/networks/edata_full.Rdata
        - ./results/Rdata/networks/wgcna_4.RDS
        -./results/Rdata/networks/moduleTraitPvalue_full_4.RData
        -./results/Rdata/networks/moduleTraitCor_full_4.RData
        -./results/Rdata/networks/geneModMemAnnot_power4.RData
        -./results/Rdata/networks/GO_sft4.RData
        
    - Sex-specific WGCNA networks constructed in ./src/WGCNA_sex_specific.R
      -Gene ontology for sex-specific WGCNA network modules also performed here
    
     -outputs:
        -./results/Rdata/networks/edata_m.Rdata
        -./results/Rdata/networks/edata_f.Rdata
        -./results/Rdata/networks/wgcna_m_5.RDS
        -./results/Rdata/networks/wgcna_f_4.RDS
        -./results/Rdata/networks/geneModMemAnnot_m_power5.RData
        -./results/Rdata/networks/geneModMemAnnot_f_power4.RData
        -./results/Rdata/networks/moduleTraitPvalue_f.RData
        -./results/Rdata/networks/moduleTraitPvalue_m.RData
        -./results/Rdata/networks/moduleTraitCor_f.RData
        -./results/Rdata/networks/moduleTraitCor_m.RData
        -./results/Rdata/networks/GO_Females_sft4.RData





    
    -Bayesian networks learned with ./src/learn_bn.R
    
      -This was run on Rivanna supercomputing cluster
      -output in ./results/Rdata/networks/bn_4/

    -Key Driver Analysis performed in ./src/KDA_working.R
      -output:
        ./results/flat/key_driver_analysis.csv
    
    -munge to add coloc data
    

        


