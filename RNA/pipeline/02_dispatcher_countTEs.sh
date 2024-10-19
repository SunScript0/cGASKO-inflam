#!/bin/bash
# Job settings
partition="standard"
max_time="6:00:00"
cores=8
mem="40G"

# File list
# files="$(find /gpfs/fs2/scratch/fmorandi/internal/RE_clock/data/rna_seq_ucar/04_mapped -name "*out.bam")"
files="$(find /gpfs/fs2/scratch/fmorandi/internal/John/L1KD/rna_seq2/04_mapped -name "*out.bam")"

echo "$files"

# Script settings
script="/home/fmorandi/Documents/scripts/pipelines/rte_rna/countTEs.sh"
genes="/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/Mus_musculus.GRCm39.108.gtf"
# genes="/scratch/fmorandi/external/references/GRCh38-hg38-Ensembl/Homo_sapiens.GRCh38.109.gtf"
# genes="/scratch/fmorandi/external/references/GRCh37-hg19-Ensembl/Homo_sapiens.GRCh37.87.gtf"
repeats="/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/GRCm39_Ensembl_rmsk_TE.gtf"
# repeats="/scratch/fmorandi/external/references/GRCh38-hg38-Ensembl/GRCh38_Ensembl_rmsk_TE.gtf"
# repeats="/scratch/fmorandi/external/references/GRCh37-hg19-Ensembl/GRCh37_Ensembl_rmsk_TE.gtf"
# outpath="/scratch/fmorandi/internal/RE_clock/data/rna_seq_ucar"
outpath="/scratch/fmorandi/internal/John/L1KD/rna_seq2"

cd $outpath
mkdir -p B_TEcounts_logs

for f in $files
do
	fname=$(basename $f | sed "s/_Aligned.sortedByCoord.out.bam//")
	logf="./B_TEcounts_logs/$(date +%Y-%m-%d)_$fname.txt"
	sbatch --partition=$partition --time=$max_time -c $cores --mem=$mem --output=$logf $script -f "$f" -n "$fname" -g "$genes" -r "$repeats" -c $cores -o "$outpath"
	sleep 0.2
done