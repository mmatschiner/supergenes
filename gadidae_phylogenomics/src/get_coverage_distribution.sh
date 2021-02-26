# m_matschiner Tue Sep 25 09:23:29 CEST 2018

# Make the output directory if it doesn't exist yet.
mkdir -p ../log/mapping

# Get the coverage distribution of all bam files.
for bam in ../res/mapping/??????.bam
do
    id=`basename ${bam%.bam}`
    out_file="../log/mapping/coverage.${id}.out"
    rm -f {out_file}
    sbatch -o ${out_file} get_coverage_distribution.slurm ${bam} ../data/assemblies/gadMor2.fasta
done