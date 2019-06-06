# DO_project pipeline

1. Munge phenotypic files into one nice file.
    - use ./src/munge_raw_pheno_data.R
        - This takes pheno files from ./data/pheno_data and outputs ./results/flat/full_pheno_table.csv and results/Rdata/full_pheno_table.Rdata

2. Fix GeneSeek files
    There was a sample confusion. Re-genotyped samples have a .1 appended to their name, except 371.
        - In ./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20190228_FinalReport, change Sample ID 371 to 371.1
            - This is done in ./src/fix_geneseek.R
                -Add header back manually. See R file for instructions
        - Fix/create ./data/GIGAMUGA/merged/Merged_Sample_Map.txt manually. Change the second 371 entry to 371.1
            - create Merged_Sample_Map.txt by merging all sample maps from geneseek experiments and editing indices appropriately. Done in TextEdit.
        - Generate ./data/GIGAMUGA/merged/Merged_FinalReport.txt. 
            - This is done in ./src/fix_geneseek.R

2. Convert Geneseek FinalReports files to a format r/qtl2 understands. This encodes the DO genotypes.
    - From https://kbroman.org/qtl2/pages/prep_do_data.html
    - use ./src/geneseek2qtl2_mod.R
    - output is in results/GIGAMUGA/qtl_batches1-4

3. Remove bad samples and markers using Argyle.
    - use ./src/GIGAMUGA_QC.Rmd
    - output ./results/GIGAMUGA/geno.final_merged.RDS


4. Create cross file for QTL mapping. 
    - ./src/make_crossfile.R

5. calculate genotype and allele probs, as well as kinship matrices 
    - how is kinship calculated?
    - ./src/calc_probs.R

6. more QC - Broman et al.
    - script adapted from https://github.com/kbroman/Paper_MPPdiag/blob/master/R/diagnostics.Rmd
    - cite paper
    - This is largely buggy, but everything looked fine. Led to the removal of several hundred markers that had greater than 5% genotyping errors.
    - output is cross_basic_cleaned (./results/Rdata/cross_basic_cleaned.Rdata)
    - 109,427 markers remaining

7. recalculate genotype and allele probs using the new cross_basic_cleaned
    - ./src/calc_probs.R

8. Calc QTL
    - ./src/map_qtl.R
    - will also perform permutation analysis

9. Calc eQTL

    - separate instructions for sequencing data preprocessing?
    - PEER stuff

    - after RNA-seq preprocessing, normalize the counts
        - ./src/normalize_RNAseq.R

    - make an eqtl specific cross file
        - ./src/make_cross_eqtl.R

    - map Cis and trans differently, using different covariates
10. Calc eQTL perms

    - Also separately for cis and trans

11. eQTL analysis: cis and trans eqtl, trans eQTL hotspots, mediation analysis, NEO structural equation  modeling

12. Networks
