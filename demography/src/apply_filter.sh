# michaelm Thu Mar 5 18:41:38 CET 2020

# Set the mask for repeat regions.
cd ../data
tar -xzf masks.tgz
rm -r masks.tgz
cd -
mask=../data/masks/Atlantic_cod_repeats.tab

# Set the reference.
mkdir -p ../data/assemblies
cp ../../cod_phylogenomics/data/assemblies/gadMor2.fasta ../data/assemblies
ref=../data/assemblies/gadMor2.fasta

# Generate filtered versions of all vcf files.
for vcfgz in ../res/gatk/LG??.vcf.gz
do
    lg=`basename ${vcfgz%.vcf.gz}`
    filtered_vcfgz=${vcfgz%.vcf.gz}.filtered.vcf.gz
    if [ ! -f ${filtered_vcfgz} ]
    then
        out=../log/gatk/filter.${lg}.out
        log=../log/gatk/filter.${lg}.log
        # Remove log files if they exist already.
        rm -f ${out}
        rm -f ${log}
        sbatch -o ${out} apply_filter.slurm ${vcfgz} ${ref} ${mask} ${filtered_vcfgz} ${log}
    fi
done
