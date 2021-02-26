#!/bin/sh

# JMIB [jmi.barth@gmail.com], 2020-08-17
# Using RELATE to produce .anc (trees per region) and .mut files (SNP information) and estimate population sizes.

# Load modules.
module purge
module load BCFtools/1.9-foss-2018b
module load R/3.5.1-foss-2018b

# Run this script using:
# bash run_relate.sh [relative_path_to_chromosome] [no._of_chr.] [relative_path_to_map] [relative_path_to_poplabels_file] 
#[relative_path_to_output directory] [mutation rate] [estimate_of_ne] [seed] [total_number_of_chromosomes_or_last_chr._analyzed] [years_per_generation]
#[fraction of trees to be dropped]

## Get command line arguments, set variables, and define output directory.
chr=$1 #../data/vcf/NC_031965.f5.2.vcf.gz
	tmp=${chr%.vcf.gz}
	chr_name=${tmp##*/}
n=$2 #1
map=$3 #../data/map/NC_031965.flat.map
poplabels=$4 #../data/ids/samples.poplabels.txt
out_dir=$5
mutation_rate=$6 #7.1496e-9
ne=$7 #100000
seed=$8 #877
n_total=$9
gen=${10}
threshold=${11}

## Feedback
echo "run_relate.sh was called with the following arguments:"
echo "Chromosome": $chr
echo "Chromosome name:" $chr_name
echo "Number of chromosome": $n
echo "Path to map": $map
echo "Path to poplables": $poplabels
echo "Path to output dir": $out_dir
echo "Mutation rate": $mutation_rate
echo "Ne": $ne
echo "Seed": $seed
echo "Total number of chromosomes/last chromosome": $n_total
echo "Generation time in years": $gen
echo "Threshold": $threshold

## Convert to SHAPEIT's haplotype format
echo "Converting to SHAPEIT's haplotype format..."
	echo "../bin/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps ${out_dir}${chr_name}.haps --sample ${out_dir}${chr_name}.sample -i ${chr%.vcf.gz}"
		  ../bin/relate/bin/RelateFileFormats --mode ConvertFromVcf --haps ${out_dir}${chr_name}.haps --sample ${out_dir}${chr_name}.sample -i ${chr%.vcf.gz}
echo "...done."

## Gzip haps file
echo "Gzip haps file..."
	 bgzip -c ${out_dir}${chr_name}.haps > ${out_dir}${chr_name}.haps.gz
		rm -f ${out_dir}${chr_name}.haps
echo "...done."

## Run RELATE to infer genome-wide genealogies and mutations (Note: output needs to be in working directory)
echo "Running RELATE to infer genome-wide genealogies and mutations..."
	echo "../bin/relate/bin/Relate --mode All -m ${mutation_rate} -N ${ne} --seed ${seed} --haps ${out_dir}${chr_name}.haps.gz --sample ${out_dir}${chr_name}.sample --map ${map} -o ${chr_name}"
		  ../bin/relate/bin/Relate --mode All -m ${mutation_rate} -N ${ne} --seed ${seed} --haps ${out_dir}${chr_name}.haps.gz --sample ${out_dir}${chr_name}.sample --map ${map} -o ${chr_name}
echo "...done."

## Bgzip and re-name anc and mut output files
echo "gzip and re-name anc and mut output files..."
	 bgzip -c ${chr_name}.mut > ${out_dir}${chr_name}_chr${n}.mut.gz
		rm -f ${chr_name}.mut
	 bgzip -c ${chr_name}.anc > ${out_dir}${chr_name}_chr${n}.anc.gz
		rm -f ${chr_name}.anc
echo "...done."

## Run RELATE to infer estimate population sizes
echo "Running RELATE to infer estimate population sizes..."
	echo "../bin/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i ${out_dir}${chr_name} --first_chr ${n} --last_chr ${n_total} -m ${mutation_rate} --poplabels ${poplabels} --years_per_gen ${gen} --threshold ${threshold} --num_iter 5 --seed ${seed} -o ${out_dir}${chr_name}.ne"
	  	  ../bin/relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i ${out_dir}${chr_name} --first_chr ${n} --last_chr ${n_total} -m ${mutation_rate} --poplabels ${poplabels} --years_per_gen ${gen} --threshold ${threshold} --num_iter 5 --seed ${seed} -o ${out_dir}${chr_name}.ne
echo "...done."



