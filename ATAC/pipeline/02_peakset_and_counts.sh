#!/bin/bash
#SBATCH --partition=standard --time=10:00:00 -c 4 --mem=40G --output=peakset_clock.log

module load bedtools
module load subread

basepath=/scratch/fmorandi/internal/John/cGAS_KO/PAPER/ATAC/pipeline_out
threads=4

cd $basepath

date

mkdir -p 07_peaksets_and_tables

# Concatenate all narrowPeak files and sort
cat ./06_macs2_outputs/*_peaks.narrowPeak | cut -f1,2,3,4 | bedtools sort > ./07_peaksets_and_tables/peaks_cat.bed
# Bedtools merge on concatenated peak file
bedtools merge -i ./07_peaksets_and_tables/peaks_cat.bed > ./07_peaksets_and_tables/peaks_merged.bed

# Multiinter on all narrowPeak files
bedtools multiinter -i ./06_macs2_outputs/*.narrowPeak > ./07_peaksets_and_tables/multinter.txt
# Select regions which were called as peaks in at least 2 samples (careful: not greater or equal so 1)
cat ./07_peaksets_and_tables/multinter.txt | awk '$4>1' > ./07_peaksets_and_tables/multinter_thresholded.txt

# Only keep peaks which overlap regions called as peak in at least n samples
bedtools intersect -wa -u -sorted \
	-a ./07_peaksets_and_tables/peaks_merged.bed \
	-b ./07_peaksets_and_tables/multinter_thresholded.txt \
	> ./07_peaksets_and_tables/peaks_merged_filt.bed
	
# Make SAF for featureCounts
awk 'OFS="\t" {print $1":"$2"-"$3, $1, $2, $3, "."}' ./07_peaksets_and_tables/peaks_merged_filt.bed > ./07_peaksets_and_tables/peaks_filt.saf
# Count 5' ends of reads over peaks
featureCounts -F "SAF" -p -B --read2pos 5 -T $threads \
	-a ./07_peaksets_and_tables/peaks_filt.saf \
	-o ./07_peaksets_and_tables/counts.tsv \
	./05_clean_BAMs/*.bam
	
# rm ./07_peaksets_and_tables/peaks_cat.bed
# rm ./07_peaksets_and_tables/multinter.txt
# rm ./07_peaksets_and_tables/multinter_thresholded.txt
# rm ./07_peaksets_and_tables/peaks_merged_filt.bed

date