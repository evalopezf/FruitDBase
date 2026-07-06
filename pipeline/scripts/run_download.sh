#!/bin/bash
#
#SBATCH -p generic
#SBATCH -J download
#SBATCH --chdir=/lustre/home/cebas/emlopez/rnaSeqAnalysis/test-run-snakemake/
#SBATCH -o /lustre/home/cebas/emlopez/rnaSeqAnalysis/test-run-snakemake/output_%j.log
#SBATCH --error=/lustre/home/cebas/emlopez/rnaSeqAnalysis/test-run-snakemake/error_%j.log
#SBATCH --cpus-per-task=6
#SBATCH --mail-user=eva.lopezf00@gmail.com
#SBATCH --mail-type=END


bash "/lustre/home/cebas/emlopez/rnaSeqAnalysis/scripts-paper/workflow/scripts/download.sh" "/lustre/home/cebas/emlopez/rnaSeqAnalysis/test-run-snakemake/newMetadata20260701.tsv" "/lustre/home/cebas/emlopez/rnaSeqAnalysis/test-run-snakemake/fastq/"
