# DO_project pipeline

1. Munge phenotypic files into one nice file.
    - use ./src/munge_raw_pheno_data.R
        - This takes pheno files from ./data/pheno_data and outputs ./results/flat/full_pheno_table.csv and results/Rdata/full_pheno_table.Rdata

2. Convert Geneseek FinalReports files to a format r/qtl2 understands. This encodes the DO genotypes.
    - From https://kbroman.org/qtl2/pages/prep_do_data.html
    - use ./src/geneseek2qtl2_mod.R
    - output is in results/GIGAMUGA/qtl_batches1-4
    - Also create Merged_Sample_Map.txt by merging all sample maps from geneseek experiments and editing indices appropriately. Done in TextEdit. ./data/GIGAMUGA/merged/
    - Also create merged_FinalReport.txt by merging all finalReport files from geneseek experiments . Done in TextEdit. In ./data/GIGAMUGA/merged/


3. Remove bad samples and markers using Argyle.
    - use ./src/GIGAMUGA_QC.Rmd
    - output ./results/GIGAMUGA/geno.final_merged.RDS"