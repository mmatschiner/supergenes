# m_matschiner Tue Mar 10 10:48:07 CET 2020

# Load modules.
module --quiet purge
module load Python/3.8.2-GCCcore-9.3.0
module load GSL/2.6-GCC-9.3.0
module load R/4.0.0-foss-2020a

# Make the result directory.
mkdir -p ../res/msprime

# Make the log directory.
mkdir -p ../log/msprime

# Simulate a bottleneck and get the proportion of the genome coalescing in it.
recombination_rate=1e-8
mutation_rate_per_myr=`cat ../res/evo/rate.txt | tail -n 1 | cut -d ":" -f 2`
mutation_rate_per_yr=`echo "${mutation_rate_per_myr} / 1000000" | bc -l`
generation_time=10
bottleneck_population_size=1
chr_length=20000000
out=../log/msprime/pi.out
res=../res/msprime/pi.txt
rm -f ${out}
sbatch -o ${out} simulate_split_with_bottleneck.slurm ${recombination_rate} ${mutation_rate_per_yr} ${generation_time} ${chr_length} ${res}

# Simulate a bottleneck and get the proportion of the genome coalescing in it.
echo -e "bottleneck_time\tpop_size\tn_trees\tmean_height\tp_bottleneck" > ../res/msprime/p_bottleneck.txt
recombination_rate=1e-8
mutation_rate_per_myr=`cat ../res/evo/rate.txt | tail -n 1 | cut -d ":" -f 2`
mutation_rate_per_yr=`echo "${mutation_rate_per_myr} / 1000000" | bc -l`
generation_time=10
bottleneck_population_size=1
chr_length=20000000
for bottleneck_time in 1000 3000 10000 30000 100000 300000 1000000 3000000
do
    for population_size in 1000 3000 10000 30000 100000
    do
        for rep in {1..20}
        do
            python simulate_bottleneck.py ${recombination_rate} ${mutation_rate_per_yr} ${generation_time} ${bottleneck_time} ${population_size} ${chr_length} >> ../res/msprime/p_bottleneck.txt
        done
    done
done

# Get the recombination map.
wget -O tmp.map.txt https://raw.githubusercontent.com/johnbowes/CRAFT-GP/master/source_data/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr20.txt

# Convert the nexus-format tree of the lg02 inversion to newick format.
Rscript convert_nexus_tree_to_newick.r ../data/trees/inversion_lg02.tre tmp.nwk.tre

# Make a temporary list of species to keep in trees.
cat ../data/tables/color_code.txt | cut -f 1 > tmp.keep_list.txt

# Prune the tree to remove outgroups.
Rscript prune_newick.r tmp.nwk.tre tmp.keep_list.txt tmp.pruned.nwk.tre

# Simulate evolution along the tree for lg02 with an extreme bottleneck, and using a recombination map (later to be analysed with a flat map to test the effect of map misspecification)
recombination_map=tmp.map.txt
chr_length=`tail -n 1 ${recombination_map} | cut -f 2`
mutation_rate_per_myr=`cat ../res/relate/input/rate.txt | tail -n 1 | cut -d ":" -f 2`
mutation_rate_per_yr=`echo "${mutation_rate_per_myr} / 1000000" | bc -l`
echo "Using ${mutation_rate_per_yr} as the mutation rate per bp and per year."
generation_time=10
population_size=30000
vcf=../res/msprime/lg02_with_recmap.vcf
python simulate_lg02_with_recmap.py tmp.pruned.nwk.tre ${recombination_map} ${mutation_rate_per_yr} ${generation_time} ${population_size} tmp.vcf ${chr_length}

# Sort samples in vcf alphabetically.
module --quiet purge
module load BCFtools/1.10.2-GCC-8.3.0
bcftools query -l tmp.vcf | sort > tmp.sample_order.txt
bcftools view -S tmp.sample_order.txt tmp.vcf > ${vcf}

# Clean up.
rm -f tmp.map.txt
rm -f tmp.nwk.tre
rm tmp.keep_list.txt
rm tmp.pruned.nwk.tre
rm -f tmp.vcf
rm -f tmp.sample_order.txt
