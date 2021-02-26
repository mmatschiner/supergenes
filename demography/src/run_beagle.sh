# m_matschiner Fri Jun 5 21:25:09 CEST 2020

# Make the output directory.
mkdir -p ../res/beagle

# Make the log directory.
mkdir -p ../log/beagle

# Make the bin directory.
mkdir -p ../bin

# Download beagle.
if [ ! -f ../bin/beagle.18May20.d20.jar ]
then
    wget https://faculty.washington.edu/browning/beagle/beagle.18May20.d20.jar
    mv beagle.18May20.d20.jar ../bin
fi

# Run beagle per linkage group.
for gzvcf in ../res/gatk/LG??.filtered.vcf.gz
do
    gzvcf_id=`basename ${gzvcf%.vcf.gz}`
    phased_gzvcf=../res/beagle/${gzvcf_id}.vcf.gz
    out="../log/beagle/${gzvcf_id}.out"
    rm -f ${out}
    sbatch -o ${out} run_beagle.slurm ${gzvcf} ${phased_gzvcf}
done
