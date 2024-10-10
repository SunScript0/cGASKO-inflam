#!/bin/bash
# Job settings
partition="standard"
max_time="48:00:00" # 24h
cores=6 #2
mem="40G" #30G

# File list
files="$(find /scratch/fmorandi/internal/John/cGAS_KO/PAPER/ATAC/pipeline_out/00_fastq -name "*.fq.gz")"

# Script settings
script="/home/fmorandi/Documents/scripts/pipelines/atac_seq/map_and_clean.sh"
genome="/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/Bowtie2/GRCm39"
outpath="/scratch/fmorandi/internal/John/cGAS_KO/PAPER/ATAC/pipeline_out"
r1_pattern="_1\.fq.gz"
r2_pattern="_2\.fq.gz"

cd $outpath
mkdir -p A_mapping_logs

r1_files=$(echo "$files" | grep "$r1_pattern")

for r1f in $r1_files
do
	r2f=$(echo $r1f | sed "s/$r1_pattern/$r2_pattern/")
	if [ $(echo $files | grep -c $r2f) -ne 1 ]; then
		echo "R2 file of $r1f not found"
		continue
	fi
	fname=$(basename $r1f | sed "s/$r1_pattern//")
	logf="${outpath}/A_mapping_logs/$(date +%Y-%m-%d)_$fname.txt"
	sbatch --partition=$partition --time=$max_time -c $cores --mem=$mem --output=$logf $script -1 "$r1f" -2 "$r2f" -n "$fname" -g "$genome" -c $cores -o "$outpath"
	sleep 0.2
done