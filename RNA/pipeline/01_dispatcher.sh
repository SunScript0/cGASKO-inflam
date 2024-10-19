#!/bin/bash

# Job settings
partition="standard"
max_time="8:00:00"
cores=4
mem="40G"

# File list
# files="$(find /gpfs/fs2/scratch/fmorandi/internal/RE_clock/data/rna_seq/00_fastq -name "*fastq.gz")"
# files="$(find /gpfs/fs2/scratch/fmorandi/extras/ChromAcc-clock/data/fastq_ucar_rna -name "*fastq.gz")"
# files="$(find /gpfs/fs2/scratch/fmorandi/internal/John/L1KD/rna_seq2/00_fastq -name "*fastq.gz")"
files="$(find /gpfs/fs2/scratch/fmorandi/internal/Huiru/00_fastq -name "*fastq.gz")"
echo "$files"
# exit

# Script settings
single_end=false
script="/home/fmorandi/Documents/scripts/pipelines/rte_rna/trim_and_map_STAR.sh"
script_se="/home/fmorandi/Documents/scripts/pipelines/rte_rna/trim_and_map_STAR_se.sh"
# genome="/scratch/fmorandi/external/references/GRCh37-hg19-Ensembl/STAR_with_SJs"
# genome="/scratch/fmorandi/external/references/GRCh38-hg38-Ensembl/STAR_with_SJs"
# genome="/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/STAR_with_SJs"
genome="/scratch/fmorandi/external/references/Degu/STAR_with_SJs"

outpath="/scratch/fmorandi/internal/Huiru"
r1_pattern="_1\.fastq\.gz"
r2_pattern="_2\.fastq\.gz"
# r1_pattern="_R1\.fastq\.gz"
# r2_pattern="_R2\.fastq\.gz"

cd $outpath
mkdir -p A_mapping_logs

if [ "$single_end" = false ]; then
	echo "Submitting in paired-end mode"
	r1_files=$(echo "$files" | grep "$r1_pattern")
	for r1f in $r1_files; do
		r2f=$(echo $r1f | sed "s/$r1_pattern/$r2_pattern/")
		if [ $(echo $files | grep -c $r2f) -ne 1 ]; then
			echo "R2 file of $r1f not found"
			continue
		fi
		fname=$(basename $r1f | sed "s/$r1_pattern//")
		logf="./A_mapping_logs/$(date +%Y-%m-%d)_$fname.txt"
		sbatch --partition=$partition --time=$max_time -c $cores --mem=$mem --output=$logf $script -1 "$r1f" -2 "$r2f" -n "$fname" -g "$genome" -c $cores -o "$outpath"
		sleep 0.2
	done
elif [ "$single_end" = true ]; then
	echo "Submitting in single-end mode"
	for r1f in $files; do
		fname=$(basename $r1f | sed "s/.fastq.gz//")
		logf="./A_mapping_logs/$(date +%Y-%m-%d)_$fname.txt"
		sbatch --partition=$partition --time=$max_time -c $cores --mem=$mem --output=$logf $script_se -1 "$r1f" -n "$fname" -g "$genome" -c $cores -o "$outpath"
		sleep 0.2
	done
fi

