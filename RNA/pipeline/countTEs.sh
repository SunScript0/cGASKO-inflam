#!/bin/bash

# Required modules
source /software/miniconda3/4.12.0/bin/activate TEtranscripts

# Initialize options
file=""
fname=""
genes=""
repeats=""
outpath=""
threads=1

# Collect options
while getopts "f:n:g:r:c:o:" option; do
	case $option in
		f) file="$OPTARG";;
		n) fname="$OPTARG";;
		g) genes="$OPTARG";;
		r) repeats="$OPTARG";;
		c) threads="$OPTARG";;
		o) outpath="$OPTARG";;
		\?) # Invalid option
			echo "Error: Invalid option"
			exit;;
	esac
done

# Check that required options were passed
if [[ -z $file ]]; then
	echo "File unspecified"
	exit 22
elif [[ -z $fname ]]; then
	echo "File name unspecified"
	exit 22
elif [[ -z $genes ]]; then
	echo "Genes GTF unspecified"
	exit 22
elif [[ -z $repeats ]]; then
	echo "Repeats GTF unspecified"
	exit 22
elif [[ -z $outpath ]]; then
	echo "Output path unspecified"
	exit 22
fi

# Actual execution
cd $outpath

echo "File: $file"
echo "Name: $fname"
echo "Genes GTF: $genes"
echo "Repeats GTF: $repeats"
date

mkdir -p 05_counts
TEcount --sortByPos -b $file \
	--GTF $genes \
	--TE $repeats \
	--project ./05_counts/$fname #Was just $fname
	
date