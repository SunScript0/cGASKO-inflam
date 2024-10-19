#!/bin/bash

# Required modules
module load fastqc
module load trimgalore
module load cutadapt
module load rnastar
module load samtools

# Initialize options
r1f=""
r2f=""
fname=""
genome=""
outpath=""
threads=1

# Collect options
while getopts "1:2:n:g:c:o:" option; do
	case $option in
		1) r1f="$OPTARG";;
		2) r2f="$OPTARG";;
		n) fname="$OPTARG";;
		g) genome="$OPTARG";;
		c) threads=$OPTARG;;
		o) outpath="$OPTARG";;
		\?) # Invalid option
			echo "Error: Invalid option"
			exit;;
	esac
done

# Check that required options were passed
if [[ -z $r1f ]]; then
	echo "R1 file unspecified"
	exit 22
elif [[ -z $r2f ]]; then
	echo "R2 file unspecified"
	exit 22
elif [[ -z $fname ]]; then
	echo "File name unspecified"
	exit 22
elif [[ -z $genome ]]; then
	echo "Genome unspecified"
	exit 22
elif [[ -z $outpath ]]; then
	echo "Output path unspecified"
	exit 22
fi

# Actual execution
cd $outpath

echo "R1: $r1f"
echo "R2: $r2f"
echo "Name: $fname"
echo "Genome: $genome"
date

mkdir -p 01_raw_fastqc_reports
fastqc -v
fastqc $r1f --outdir=./01_raw_fastqc_reports --threads=$threads
fastqc $r2f --outdir=./01_raw_fastqc_reports --threads=$threads

mkdir -p 02_trimmed_fastq
echo "TrimGalore version:"
trim_galore -v
trim_galore --output_dir ./02_trimmed_fastq --paired $r1f $r2f --basename "${fname}"

mkdir -p 03_trimmed_fastqc_reports
fastqc "./02_trimmed_fastq/${fname}_R1_val_1.fq.gz" --outdir=./03_trimmed_fastqc_reports --threads=$threads
fastqc "./02_trimmed_fastq/${fname}_R2_val_2.fq.gz" --outdir=./03_trimmed_fastqc_reports --threads=$threads

mkdir -p 04_mapped
echo "STAR " $(STAR --version)
STAR --readFilesIn "./02_trimmed_fastq/${fname}_R1_val_1.fq.gz" "./02_trimmed_fastq/${fname}_R2_val_2.fq.gz" \
	 --genomeDir "$genome" --outFileNamePrefix "./04_mapped/${fname}_" \
	 --runThreadN $threads --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
	 --outFilterMultimapNmax  100 \
	 --winAnchorMultimapNmax 200 \
	 --outFilterMismatchNoverLmax 0.04 

date
