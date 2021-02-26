#!/bin/sh

# JMIB [jmi.barth@gmail.com], 2020-06-10
# Using RELATE to infere populations size and demography

#SBATCH --job-name=rel_demog
#SBATCH --account=XXX
#SBATCH --time=23:59:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset
	module load BCFtools/1.9-foss-2018b
	module load R/3.5.1-foss-2018b
module list    # For easier debugging

## Set variables, and define output directory.
LG=$1 #LG01 > specify in command line
n=`echo ${LG} | sed 's/LG0//g' | sed 's/LG//g'` #1
poplabels=samples.poplabels.txt
out_dir=`readlink -f ../res/relate/output/inversion_ancestor_ne`
mkdir -p ${out_dir}

## Copy input files to the work directory:
echo -n "Copying infiles..."
	cp ../res/relate/output/inversion_ancestor_anc_mut/${LG}_inversion.bial.renamed.preped_chr${n}.anc.gz $SCRATCH
	cp ../res/relate/output/inversion_ancestor_anc_mut/${LG}_inversion.bial.renamed.preped_chr${n}.mut.gz $SCRATCH
	cp ../data/tables/${poplabels} $SCRATCH
	cp -r ../bin/relate/ $SCRATCH
echo " done."

## Change to temporary directory.
cd $SCRATCH

## Run RELATE to infer estimate population sizes (note, LG1, 2, 7, 12 excluded, so only 19 chromosomes!)
echo "Running RELATE to infer estimate population sizes..."
echo "./relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i ${LG}_inversion.bial.renamed.preped --first_chr ${n} --last_chr ${n} -m 1.63922e-8 --poplabels ${poplabels} --years_per_gen 10 --threshold 0.5 --num_iter 5 --seed 1 -o ${LG}_inversion.bial.renamed.preped.ne"
	  ./relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i ${LG}_inversion.bial.renamed.preped --first_chr ${n} --last_chr ${n} -m 1.63922e-8 --poplabels ${poplabels} --years_per_gen 10 --threshold 0.5 --num_iter 5 --seed 1 -o ${LG}_inversion.bial.renamed.preped.ne
echo "...done."

## Copy files back
	echo -n "Copying result files back to the submission directory..."
		cp *inversion.bial.renamed.preped.ne* ${out_dir}
	echo " done."
