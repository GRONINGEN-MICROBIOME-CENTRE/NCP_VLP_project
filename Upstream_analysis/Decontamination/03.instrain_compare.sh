cd ./ALL_PROFILE

ml Python/3.10.8-GCCcore-12.2.0
source /scratch/p282752/tools/python_envs/instrain/bin/activate
ml SAMtools

inStrain compare \
	-cov 0.5 \
	--min_cov 1 \
	-i ./* \
	-o ../instrain_compare_ALL \
	-sc ../vOTUs_shared_by_samples_and_ncs \
	-s ../genomes.stb \
	--database_mode 
deactivate
