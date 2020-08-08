# Bone Strength in the Diversity Outbred

This project contains scripts used in our Diversity Outbred project, currently published as a [preprint](https://www.biorxiv.org/content/10.1101/2020.06.24.169839v1).


Check out the data on our [QTL viewer](http://qtlviewer.uvadcos.io). You can also download an RData file here which contains Mouse phenotypes and allele probabilities.


Sequencing data is available on GEO:
    - [RNA-seq data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152708)<br/>
    - [Single-cell RNA-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152806)


Go to [analysis_pipeline](/doc/analysis_pipeline.md) for the analysis pipeline.


Please contact Basel Al-Barghouthi (bma8ne AT virginia DOT edu) or Charles Farber (crf2s AT virginia DOT edu) with any questions, or for any data.



Folders are as follows:
    - /bin/ contains scripts that I got from other people

    - /data/contains raw data. For the most part, this is a read-only file (some exceptions)

        - /data/GIGAMUGA/ contains raw data pertaining to the GigaMUGA genotyping. Some writing has been done to this

        - /data/pheno_data/ contains phenotypic data files

    - /doc/ contains documents.

    - /results/ contains output from the analyses described. A bit messy but:

        - /results/flat/ contains flat file outputs (.txt, .csv, etc)

        - /results/GIGAMUGA/ contains output pertaining to the genotyping

        - /results/Rdata contains .Rdata and .Rds output







