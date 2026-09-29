#!/bin/bash
#
#SBATCH -p generic
#SBATCH --chdir=/lustre/home/cebas/emlopez/rnaSeqAnalysis/scripts-paper/workflow/
#SBATCH -J run_snakemake
#SBATCH --cpus-per-task=5
#SBATCH --mail-type=END
#SBATCH --mail-user=eva.lopezf00@gmail.com
#SBATCH -o /lustre/home/cebas/emlopez/rnaSeqAnalysis/updated_new_samples/output_%j.out
#SBATCH --error=/lustre/home/cebas/emlopez/rnaSeqAnalysis/updated_new_samples/my_job_%j.err


module load   rama0.4  GCC/12.2.0  OpenMPI/4.1.4
module load snakemake/7.32.3 


snakemake --rerun-incomplete --cores 1 --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem {resources.mem_mb}M --cpus-per-task {cluster.cpus} -o {cluster.output} -e {cluster.error}" --jobs 4 --use-conda

