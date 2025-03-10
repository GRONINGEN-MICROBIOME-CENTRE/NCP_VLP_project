#!/bin/bash
#SBATCH --job-name=inStrain_profile
#SBATCH --output=./out/ALL/IS_profile_%A_%a.out
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=96GB

SAMPLE_LIST=$1

echo ${SAMPLE_LIST}

SAMPLE_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${SAMPLE_LIST})

echo "SAMPLE_ID=${SAMPLE_ID}"

# COPYING THE READ ALIGNMENT MAPS:
cp ../bam_files/${SAMPLE_ID}_w_neg_der95_NCP.sorted.bam* ./bam

# FILTERING READ ALIGNMENT MAPS TO CONTAIN ONLY NC-IDENTIFIED VOTUS PRESENT PER FILE:
ml SAMtools
samtools view \
	-b \
	-o ./bam/${SAMPLE_ID}_filtered.sorted.bam \
	./bam/${SAMPLE_ID}_w_neg_der95_NCP.sorted.bam $(cat ./keep_vOTUs_per_sample/${SAMPLE_ID}_vOTUs_to_keep)

rm ./bam/${SAMPLE_ID}_w_neg_der95_NCP.sorted.bam
rm ./bam/${SAMPLE_ID}_w_neg_der95_NCP.sorted.bam.bai

samtools index \
	./bam/${SAMPLE_ID}_filtered.sorted.bam \
	-@ $((${SLURM_CPUS_PER_TASK}-1))

module purge

# --- LOADING MODULES ---
ml Python/3.10.8-GCCcore-12.2.0
source /scratch/p282752/tools/python_envs/instrain/bin/activate
ml SAMtools

inStrain profile \
	./bam/${SAMPLE_ID}_filtered.sorted.bam \
	../all_extended_pruned_viral_renamed_NCP.fasta \
	-o ./ALL_PROFILE/${SAMPLE_ID}.inStrain \
	--database_mode \
	-s ./genomes.stb \
	-p ${SLURM_CPUS_PER_TASK} \
	--scaffolds_to_profile ./vOTUs_shared_by_samples_and_ncs \
	--min_cov 1 \
	--skip_plot_generation

deactivate

module purge
