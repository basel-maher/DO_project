#RNA-seq pipeline

1. Fastq files were downloaded from Illumina Basespace. 
  * Paired-end, 2x76 bp
  * QC performed by FASTQC and multiQC.
  * Each sample was sequenced across multiple lanes and multiple runs. All samples were sequenced across 4 runs, while some (low aligners) were resequenced in a 5th run
    * For each run, al lanes for a sample were merged. For example, all 4 lanes for sample 1 R1 for run 1 were merged.
    
2. Each samples merged lanes (per run) were aligned using Hisat2. Alignment was done to the GRCm38_snp index.

  * last base was trimmed (3' end)

  ```bash 
  i=$SLURM_ARRAY_TASK_ID
  line=$(sed -n "${i}{p;}" < ../hisat_aligned_lanes/samples)

  hisat2 --dta --trim3 1 2>../hisat_aligned_lanes/aln_sums/${line}.sum -x /grcm38_snp/genome_snp -1 $(grep -e "/${line}_R1" ../hisat_aligned_lanes/merged_files) -2 $(grep -e "/${line}_R2" ../hisat_aligned_lanes/merged_files) | samtools view -bS > ../hisat_aligned_lanes/aligned/${line}_hisat_aligned_lanes.bam
  
  ```
  * Above is a sample. Basically, a file with sample names is used to find all sample_R1 and sample_R2 files, and align to GRCm38_snp. Output is a BAM file and an alignment summary file.
  
  
3. BAM files were sorted, merged across runs (per sample), and sorted again.

  ```bash
  samtools sort sample.bam -o sample.sorted
  
  samtools merge -r sample.mergedAll.bam $(grep _${sample} sorted_aligned_merged_files | sort)
  
  samtools sort ${sample} -o ${out_name}.sorted
  ```
  
4. Assembly was performed with Stringtie, to GRCm38.98.gtf

  ```bash
  stringtie -e -G ../../Mus_musculus.GRCm38.98.gtf -o ../assembled/${name}_assembled.gtf -A ../gene_abund/${name}.gene_abund.tab ${sample}
  ```

5. A single GTF was with common transcripts was generated, and assembly was performed again, but to this GTF instead.

  ```bash
  stringtie --merge -G GRCm38/Mus_musculus.GRCm38.98.gtf -o mus_stringtie_merged.gtf assembled_files
  ```

  ```bash
  stringtie -eB -G mus_stringtie_merged.gtf -o ../assembled_again/${name}/${name}_assembled.gtf -A ../gene_abund_again/${name}.gene_abund.tab ${sample}
  ```
  
6. Take all gene abundance files (1 per sample) and concatenate into a single file. Then, identify genes that had more than 0.1 transcripts per million in at least 38 samples (20% of samples)

  ```bash
  #Took gene abundances for each file, concatenated into all_gene_abundances
  cat ../gene_abund_again/*.tab > all_gene_abundances

  #used AWK to filter on more than 0.1 TPM (all_gene_abund_filt_0.1tpm)
  awk '($9 > 0.1) {print $0 >> "all_gene_abund_filt_0.1tpm"}' all_gene_abundances

  #counted how many pass filter (counts_pass_0.1_filt), 
  awk '{a[$1]++;} END{for(i in a) print a[i]"  "i >> "counts_pass_0.1_filt"}' all_gene_abund_filt_0.1tpm

  #took genes that pass in more than 38 samples (~20%)
  awk '($1 > 38) {print $2 >> "pass_0.1_filt_over_38samps"}' counts_pass_0.1_filt
  ```
  
7. Create a file with the re-assembled GTF names, then use prepDE.py to generate count matrices.

  ```bash
  find ../assembled_again/ -type f | grep -E "*gtf" > prepDE_list.txt

  awk -F'[/-]' '{print $9, $0 > "prepDE_list.txt"}' prepDE_list.txt

  #generate gene_count_matrix and transcript_count_matrix (length 75)
  python prepDE.py -i prepDE_list.txt -l 75
  ```
  
8. Use R script to find number of samples per gene with more than 6 reads. Keep those more than 6 reads in more than 38 samples (20%).
Then, keep genes that pass the TPM filter in step 6 above.

```R
  x = read.csv("gene_count_matrix.csv", stringsAsFactors=F, header=T, check.names=F)

  genes_over6reads_over38samps = c()

  for(i in 1:nrow(x)){
  if (length(which(x[i,2:ncol(x)] > 6))>38){
  genes_over6reads_over38samps= append(genes_over6reads_over38samps, x[i,1])
  }
  }

  tpm = read.csv("pass_0.1_filt_over_38samps", stringsAsFactors=F)

  z=genes_over6reads_over38samps[which(genes_over6reads_over38samps %in% tpm[,1])]

  write(z,"genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps.txt", sep="\t")

  final=x[which(x[,1] %in% z),]

  colnames(final)=colnames(x)

  write.csv(final, "genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS.csv",row.names=F, quote=F)
```


  

    