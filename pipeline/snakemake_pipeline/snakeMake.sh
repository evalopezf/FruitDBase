#!/bin/bash
#
#SBATCH -p generic
#SBATCH --chdir=
#SBATCH -J run_snakemake
#SBATCH --cpus-per-task=5
#SBATCH -o 
#SBATCH --error=

module load   rama0.4  GCC/12.2.0  OpenMPI/4.1.4
module load snakemake/7.32.3 


snakemake --rerun-incomplete --cores 1 --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem {resources.mem_mb}M --cpus-per-task {cluster.cpus} -o {cluster.output} -e {cluster.error}" --jobs 4 --use-conda

