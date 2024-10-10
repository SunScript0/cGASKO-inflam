#!/bin/bash
# Job settings
partition="preempt"
max_time="6:00:00"
cores=1
mem="5G"

files="$(find "/scratch/fmorandi/internal/John/cGAS_KO/PAPER/ATAC/pipeline_out/05_clean_BAMs" -name "*_clean.bam" | sed 's/_clean.bam//')"

echo "$files"

# Script settings
script="/home/fmorandi/Documents/scripts/pipelines/atac_seq/make_bigwig.sh"
basepath=/scratch/fmorandi/internal/John/cGAS_KO/PAPER/ATAC/pipeline_out
genome_index=/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/Mus_musculus.GRCm39.dna.primary_assembly.fa.fai
qc_summary=/scratch/fmorandi/internal/John/cGAS_KO/PAPER/ATAC/results/qc_summary.tsv

cd $basepath
mkdir -p B_bigwig_logs

for f in $files
do
	fname=$(basename $f)
	logf="${basepath}/B_bigwig_logs/$(date +%Y-%m-%d)_$fname.txt"
	sbatch --partition=$partition --time=$max_time -c $cores --mem=$mem --output=$logf $script -f "$f" -g "$genome_index" -q "$qc_summary" -o "$basepath"
	sleep 0.2
done