#!/bin/bash

###############################################
#  Shell script for RRBS data analysis
###############################################

## Prerequisite
## Ensure the following software is available in the current environment:
## fastqc, multiqc, trim_galore, bismark, samtools, parallel

## Set directories
fastq_dir="fastq"
fastqc_dir="fastqc"; mkdir -p ${fastqc_dir}
trimmed_fastq_dir="trimmed_fastq"; mkdir -p ${trimmed_fastq_dir}
bam_dir="bam_dir"; mkdir -p ${bam_dir}
multiqc_dir="multiqc"; mkdir -p ${multiqc_dir}
genome_dir="genome_index"

## Number of parallel threads
threads=36

# 1. Quality control with FastQC
fastqc -f fastq -t ${threads} -o ${fastqc_dir} -q ${fastq_dir}/*.fastq.gz

# 2. Create a multiqc report (pre-trimming)
multiqc -o ${multiqc_dir} -f -n pre_trimming_report.html ${fastqc_dir}

# 3. Adapter trimming and quality filtering using Trim Galore
ls -1 ${fastq_dir}/*.fastq.gz | cut -d'R' -f1 | sort | uniq | parallel -j 6 --tmpdir tmp/ \
    'trim_galore --rrbs --illumina --gzip --fastqc -o ${trimmed_fastq_dir} {}R1_001.fastq.gz'

# 4. Renaming trimmed files for consistency
for file in ${trimmed_fastq_dir}/*_trimmed.fq.gz; do
    mv $file ${file/_trimmed.fq.gz/.fastq.gz}
done

# 5. Quality control with FastQC
fastqc -f fastq -t ${threads} -o ${fastqc_dir} -q ${fastq_dir}/*.fastq.gz

# 6. Create a multiqc report (pre-trimming)
multiqc -o ${multiqc_dir} -f -n after_trimming_report.html ${fastqc_dir}

# 7. Concatenate paired fastq files for each sample
for file in ${trimmed_fastq_dir}/*L001_R1_001.fastq.gz; do
    base=$(basename $file _L001_R1_001.fastq.gz)
    cat ${trimmed_fastq_dir}/${base}_L001_R1_001.fastq.gz \
        ${trimmed_fastq_dir}/${base}_L002_R1_001.fastq.gz > ${trimmed_fastq_dir}/${base}.merged.fastq.gz
done

# 8. Mapping with Bismark
ls -1 ${trimmed_fastq_dir}/*.merged.fastq.gz | cut -d'/' -f2 | parallel -j 3 --tmpdir tmp/ \
    'bismark -q --unmapped --ambiguous --genome_folder ${genome_dir} -p 4 --non_directional --multicore 3 -se ${trimmed_fastq_dir}/{}'

# 9. Sort and index BAM files
for bam in ${bam_dir}/*_bismark_bt2.bam; do
    base=$(basename $bam _bismark_bt2.bam)
    samtools sort -@ ${threads} -m 4G -o ${bam_dir}/${base}.sorted.bam ${bam_dir}/${base}_bismark_bt2.bam
    samtools index -@ ${threads} ${bam_dir}/${base}.sorted.bam
done

# 10. Create a final multiqc report (post-processing)
multiqc -o ${multiqc_dir} -f -n post_processing_report.html ${trimmed_fastq_dir} ${bam_dir}

echo "RRBS data preprocessing completed!"
