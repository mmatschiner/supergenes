# m_matschiner Sun Sep 3 19:49:36 CEST 2017

# Specify the name of the reference file.
ref=../data/assemblies/gadMor2.fasta

# Use GATKs CatVariants to concatenate all vcfs of all lgs.
concatenated_vcfgz=../res/gatk/concatenated.masked.variable.vcf.gz
if [ ! -f ${concatenated_vcfgz} ]
then
    out="../log/misc/concatenate_vcfs.out"
    rm -f ${out}
    sbatch -o ${out} concatenate_vcfs.slurm ${ref} ${concatenated_vcfgz}
fi