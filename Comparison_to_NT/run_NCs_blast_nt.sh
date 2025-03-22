#!/bin/bash

#SBATCH --job-name=NCs_blast_nt
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=/scratch/hb-llnext/VLP_public_data/nc_project/for_upload/Negativome/fragmented/nt_%A_%a.out
#SBATCH --mem=32GB

QUERY_LIST=$1

echo "QUERY list: ${QUERY_LIST}"

QUERY=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${QUERY_LIST})

module purge
module load BLAST+/2.13.0-gompi-2022a

blastn \
	-query /scratch/hb-llnext/VLP_public_data/nc_project/for_upload/Negativome/fragmented/${QUERY}.fa \
        -task blastn \
        -db /scratch/p282752/databases/NT_v1.1_nov_23/nt \
        -out /scratch/hb-llnext/VLP_public_data/nc_project/for_upload/Negativome/fragmented/${QUERY}.outfmt6.txt \
        -evalue 0.001 \
	-num_threads 2 \
	-outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle qcovs'
