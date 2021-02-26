# m_matschiner Thu Oct 31 09:44:57 CET 2019

# Set the list of samples.
samples=../data/tables/samples.txt

# Run dsuite separately for all snps, snps outside of inversion regions, and snps in inversion regions.
vcfgz=../res/gatk/concatenated.masked.variable.vcf.gz

# Prepare a file that specifies which trios should be investigated.
echo -e "gadmac\tgadoga\tgadmon" > tmp.trios.txt
echo -e "gadmac\tgadoga\tgadmoc" >> tmp.trios.txt
echo -e "borsai\tarcgla\tgadmon" >> tmp.trios.txt
echo -e "borsai\tarcgla\tgadmoc" >> tmp.trios.txt

# Run dsuite.
out=../log/misc/dsuite_dinvestigate.out
rm -f ${out}
sbatch -o ${out} run_dsuite_dinvestigate.slurm ${vcfgz} ${samples} tmp.trios.txt ../res/tables
