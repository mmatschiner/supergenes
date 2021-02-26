#!/bin/sh

# JMIB [jmi.barth@gmail.com], 2020-08-17

## Load modules.
module --quiet purge
module load BCFtools/1.9-foss-2018b

# Get the recombination map.
wget -O tmp.map.txt https://raw.githubusercontent.com/johnbowes/CRAFT-GP/master/source_data/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr20.txt

# Compress the vcf.
rm -f ../res/msprime/sim_lg02_with_recmap.vcf.gz
bgzip -c ../res/msprime/lg02_with_recmap.vcf > ../res/msprime/sim_lg02_with_recmap.vcf.gz

# Prepare directories.
mkdir -p ../res/relate/output/simulated_recmap
rm -f ../res/relate/output/simulated_recmap/*
mkdir -p ../res/relate/output/simulated_flatmap
rm -f ../res/relate/output/simulated_flatmap/*
rm -rf sim_lg02_with_recmap

# Run RELATE with the recombination map.
bash run_relate_simulated.sh ../res/msprime/sim_lg02_with_recmap.vcf.gz 1 tmp.map.txt ../data/tables/samples.poplabels.simulated.txt ../res/relate/output/simulated_recmap/ 1.63922e-8 100000 1 1 10 0.5

# Run RELATE with a flat map.
bash run_relate_simulated.sh ../res/msprime/sim_lg02_with_recmap.vcf.gz 1 ../data/maps/flat_map_GRCh37_chr20.map ../data/tables/samples.poplabels.simulated.txt ../res/relate/output/simulated_flatmap/ 1.63922e-8 100000 1 1 10 0.5

# Clean up.
rm -f tmp.map.txt
