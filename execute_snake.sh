#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH -t 168:00:00
#SBATCH --mem-per-cpu=10G # megabytes
#SBATCH -J SNAKE
#SBATCH -o SNAKE.o%j
#SBATCH -e SNAKE.e%j
#SBATCH --mail-user=Rita.Tam@anu.edu.au
#SBATCH --mail-type=ALL


source /opt/conda/etc/profile.d/conda.sh
conda activate /mnt/data/wright/home/groups/schwessinger/condaEnvs/snakemake

cd /mnt/data/wright/home/scratch/groups/schwessinger/Pst198E_funannotate/funannotate_snakemake
snakemake -c32 -j32 