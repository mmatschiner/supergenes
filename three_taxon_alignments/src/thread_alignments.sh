# m_matschiner Tue Jun 26 23:53:27 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/alignments/merged

# Make the log directory if it doesn't exist yet.
mkdir -p ../log/misc

# Merge and thread alignments for each chromosome.
for number in `seq -w 1 23`
do
    # Set file names.
    chromosome_id=lg${number}
    alignment_without_ns_gadMor_Stat=../res/alignments/masa_cudalign/pairwise_without_ns/gadMor_Stat/${chromosome_id}.fasta
    alignment_with_ns_gadMor_Stat=../res/alignments/masa_cudalign/pairwise_with_ns/gadMor_Stat/${chromosome_id}.fasta
    alignment_lastz_gadMor_Stat=../res/lastz/pairwise_alignments/gadMor_Stat/${chromosome_id}.fasta
    alignment_without_ns_melAeg=../res/alignments/masa_cudalign/pairwise_without_ns/melAeg/${chromosome_id}.fasta
    alignment_with_ns_melAeg=../res/alignments/masa_cudalign/pairwise_with_ns/melAeg/${chromosome_id}.fasta
    alignment_lastz_melAeg=../res/lastz/pairwise_alignments/melAeg/${chromosome_id}.fasta
    alignment_merged=../res/alignments/merged/${chromosome_id}.fasta
    alignent_threaded=../res/alignments/merged/${chromosome_id}.threaded.fasta

    # Merge alignments.
    if [ ! -f ${alignment_merged} ]
    then
	    cat ${alignment_without_ns_gadMor_Stat} ${alignment_with_ns_gadMor_Stat} ${alignment_lastz_gadMor_Stat} ${alignment_without_ns_melAeg} ${alignment_with_ns_melAeg} ${alignment_lastz_melAeg} > ${alignment_merged}
    fi

    # Thread alignments.
    if [ ! -f ${alignent_threaded} ]
    then
	    out=../log/misc/thread_alignment.${chromosome_id}.out
	    rm -f ${out}
	    sbatch -o ${out} thread_alignment.slurm ${alignment_merged} ${alignent_threaded}
    fi
    
done
