#RNA-seq pipeline

1. Fastq files were downloaded from Illumina Basespace. 
  * Paired-end, 2x76 bp
  * QC performed by FASTQC and multiQC.
  * Each sample was sequenced across multiple lanes and multiple runs. All samples were sequenced across 4 runs, while some (low aligners) were resequenced in a 5th run
    * For each run, al lanes for a sample were merged. For example, all 4 lanes for sample 1 R1 for run 1 were merged.
    
2. Each samples merged lanes (per run) were aligned using Hisat2. Alignment was done to the GRCm38_snp index.

  * last base was trimmed (3' end)
  * 
  ```bash 
  i=$SLURM_ARRAY_TASK_ID
  line=$(sed -n "${i}{p;}" < ../hisat_aligned_lanes/samples)

  hisat2 --dta --trim3 1 2>../hisat_aligned_lanes/aln_sums/${line}.sum -x /grcm38_snp/genome_snp -1 $(grep -e "/${line}_R1" ../hisat_aligned_lanes/merged_files) -2 $(grep -e "/${line}_R2" ../hisat_aligned_lanes/merged_files) | samtools view -bS > ../hisat_aligned_lanes/aligned/${line}_hisat_aligned_lanes.bam
  
  ```
  * Above is a sample. Basically, a file with sample names is used to find all sample_R1 and sample_R2 files, and align to GRCm38_snp. Output is a BAM file and an alignment summary file.

    