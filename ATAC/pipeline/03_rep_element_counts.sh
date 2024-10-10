#!/bin/bash
#SBATCH --partition=standard --time=10:00:00 -c 4 --mem=80G --output=re_counts.log

module load bedtools
module load subread

basepath=/scratch/fmorandi/internal/John/cGAS_KO/PAPER/ATAC/pipeline_out
rep_saf=/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/RepeatMaskerOut/GRCm39_repeats.saf
threads=4

cd $basepath

date

# Make a copy of the premade saf
tail -n +2 $rep_saf > ./07_peaksets_and_tables/rtes.saf

# # Count 5' ends of reads over REs
featureCounts -F "SAF" -p -B --read2pos 5 -T $threads \
	-a ./07_peaksets_and_tables/rtes.saf \
	-o ./07_peaksets_and_tables/counts_rtes.tsv \
	./05_clean_BAMs/*.bam
	
# Concatenate peakset and repeat regions
cat ./07_peaksets_and_tables/peaks_filt.saf ./07_peaksets_and_tables/rtes.saf > ./07_peaksets_and_tables/peaks_and_rtes.saf
	
# Count 5' ends of reads over concatenated peaks and REs
# This time I will need to account for overlapping features
featureCounts -F "SAF" -p -B --read2pos 5 -T $threads \
	-O --fraction \
	-a ./07_peaksets_and_tables/peaks_and_rtes.saf \
	-o ./07_peaksets_and_tables/counts_ocrs_and_rtes.tsv \
	./05_clean_BAMs/*.bam

date