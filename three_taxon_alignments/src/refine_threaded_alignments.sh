# m_matschiner Fri Jun 29 23:09:50 CEST 2018

# Refine the threaded alignments.
for alignment in ../res/alignments/merged/lg??.threaded.fasta
do
    refined_alignment=${alignment%.fasta}.refined.fasta
    chromosome_id=`basename ${alignment%.threaded.fasta}`
    out=../log/misc/refine_alignment.${chromosome_id}.out
    if [ ! -f ${refined_alignment} ]
    then
		rm -f ${out}
		sbatch -o ${out} refine_threaded_alignment.slurm ${alignment} ${refined_alignment}
    fi
done
