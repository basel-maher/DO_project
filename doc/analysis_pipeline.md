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

5. more QC - Broman et al.

5. calculate genotype and allele probs, as well as kinship matrices 
    - ./src/calc_probs.R