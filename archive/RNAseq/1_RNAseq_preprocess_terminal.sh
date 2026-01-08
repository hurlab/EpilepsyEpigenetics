###############################################
#  Shell script for EAT RNA-Seq data analysis
###############################################

## Prerequisite
## The folloing software should be available in the current environment
## fastqc, multiqc, samtools, hisat2, trimmomatic

## Create subdirectories
fastq_dir="fastq"
fastqc_dir="fastqc"; mkdir -p fastqc
trimmed_fastq_dir="trimmed_fastq"; mkdir -p trimmed_fastq
trimmed_fastqc_dir="trimmed_fastqc"; mkdir -p trimmed_fastqc
multiqc_dir="multiqc"; mkdir -p multiqc
bam_dir="bam_dir"; mkdir -p bam_dir

## Number of parallel processing (multi-threads)
threads=36

## Run fastqc
fastqc -f fastq -t ${threads} -o ${fastqc_dir} -q ${fastq_dir}/*.fastq.gz

## Create a multiqc report
multiqc -o ${multiqc_dir} -f -n before_trimming_report.html ${fastqc_dir}

## Ensure adapter sequence file exists in the current directory
## TruSeq3_adapters.fa contains the sequences from TruSeq3-PE-2.fa and TruSeq3-SE.fa 

## Perform trimming using trimmomatic
for i in `ls ${fastq_dir}/*.fastq.gz | cut -d'/' -f2`
do 
    base=`echo ${i} | cut -f1 -d'.'`
	trimmomatic SE -threads ${threads} -phred33 ${fastq_dir}/${base}.fastq.gz ${trimmed_fastq_dir}/${base}.clean.fastq.gz ILLUMINACLIP:TruSeq3_adapters.fa:2:30:10 HEADCROP:12 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
done

## Run fastqc for trimmed fastq files
fastqc -f fastq -t ${threads} -o ${trimmed_fastqc_dir} -q ${trimmed_fastq_dir}/*.fastq.gz

## Create a multiqc report
multiqc -o ${multiqc_dir} -f -n after_trimming_report.html ${trimmed_fastqc_dir}

## the filenames need to be slightly modifed
## '-' has some issues with batch processing below
cd ${trimmed_fastq_dir}
rename s/\-/\_/ *.gz
cd ..

## Map using hisat2
## Ensure hisat2 index are available
hisat2_index="/data/genomes/rn6/genome"

for i in `ls ${trimmed_fastq_dir}/*_R1_001.clean.fastq.gz | cut -f2 -d'/'` 
do 
	base=`echo ${i} | cut -f1 -d'.'`
    hisat2 -x ${hisat2_index} -U ${trimmed_fastq_dir}/${base}.clean.fastq.gz -p ${threads} | samtools view -Sbh > ${trimmed_fastq_dir}/${base}.bam
done

## Merge bam files per sample. There are two files from two lanes for each sample
for i in `ls ${trimmed_fastq_dir}/*L001*.bam | cut -f2 -d'/'` 
do 
	base=`echo ${i} | cut -f1 -d'_'`
    samtools merge -@ ${threads} ${bam_dir}/${base}_S${base}.merged.bam ${trimmed_fastq_dir}/${base}_S${base}_L001_R1_001.bam ${trimmed_fastq_dir}/${base}_S${base}_L002_R1_001.bam
	samtools sort -@ ${threads} -o ${bam_dir}/${base}_S${base}.sorted.bam ${bam_dir}/${base}_S${base}.merged.bam 
	samtools index -@ ${threads} ${bam_dir}/${base}_S${base}.sorted.bam
done

## Now the sorted.bam files in bam_dir directory are ready for further processing in R
