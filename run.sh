#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --partition=normal_q
#SBATCH -n 32
#SBATCH --mem=100G
#SBATCH --account=aipmm

THREADS=32

CUTADAPT="3.cutadapt"
STAR="5.star"
FC="7.features"

REF_GEN="data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
REF_RNA="data/GCF_000001405.40_GRCh38.p14_rna.fna.gz"
ANNO="data/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
REF_GEN_BUILD="${STAR}/build"

FEATURES_COUNT="${FC}/features_count.txt"
FEATURES_MATRIX="${FC}/features_matrix.txt"

# STEP 0: GET DATA: from SRR_Acc_List.txt
for i in $(cat data/SRR_Acc_List.txt)
do
    parallel-fastq-dump --outdir data --skip-technical --split-3 --gzip --threads $THREADS --sra-id $i
done


# STEP 1. QUALITY CONTROL: remove reads being shorter than 20 nucleotides, reads having quality score smaller than 20
mkdir 1.fastqc
fastqc -t 32 data/*.fastq* -o 1.fastqc
multiqc 1.fastqc -f -o 2.multiqc
mkdir $CUTADAPT
for i in data/*_1.fastq.gz
do
    SAMPLE=$(echo ${i} | sed "s/_1\.fastq\.gz//")
    in_1=${SAMPLE}_1.fastq.gz
    in_2=${SAMPLE}_2.fastq.gz
    name="$(basename $SAMPLE)"
    out_1=$CUTADAPT/${name}_1.fastq.gz
    out_2=$CUTADAPT/${name}_2.fastq.gz
    report=$CUTADAPT/${name}.txt
    cutadapt --cores=0 --minimum-length 20 --quality-cutoff 20 -o $out_1 -p $out_2 $in_1 $in_2 > $report
done
multiqc $CUTADAPT -f -o 4.multiqc

# STEP 2. MAPPING: minimum length of segments of read is 18, disable coverage search for junctions,
mkdir $STAR
STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir $REF_GEN_BUILD --genomeFastaFiles $REF_GEN --sjdbGTFfile $ANNO --sjdbOverhang 36
for i in $CUTADAPT/*_1.fastq.gz
do
    SAMPLE=$(echo ${i} | sed "s/_1\.fastq\.gz//")
    IN_1=${SAMPLE}_1.fastq.gz
    IN_2=${SAMPLE}_2.fastq.gz
    NAME="$(basename $SAMPLE)"
    BAM=$STAR/${NAME}Aligned.sortedByCoord.out.bam
    STAR --runThreadN $THREADS --readFilesCommand zcat --genomeDir $REF_GEN_BUILD --readFilesIn $IN_1,$IN_2 --outFileNamePrefix $STAR/$NAME --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN $THREADS
done
multiqc $STAR -f -o 6.multiqc

# STEP 3. COUNTING: Count reads for each gene based on the sorted bam files, minimum mapping quality is 10, reverse strand
mkdir $FC
featureCounts -s 2 -Q 10 -t exon -g gene_id -T $THREADS -a $ANNO -o $FEATURES_COUNT $STAR/*.bam
cut -d$'\t' -f2,3,4,5,6 --complement $FEATURES_COUNT > $FEATURES_MATRIX
awk -F'\t' '{ print $1,"\t",$7,"\t",$8,"\t",$9,"\t",$10}' $FEATURES_COUNT > $FEATURES_MATRIX
multiqc $FC -f -o 8.multiqc

# STEP 4. DIFFERENTIAL EXPRESSION ANALYSIS: the input is features_matrix.txt and conditions.txt
Rscript deseq2.r

# STEP 5. PATHWAY ANALYSIS (https://david.ncifcrf.gov/tools.jsp) and find possible biological pathways that related to the DE genes.
# Attached images