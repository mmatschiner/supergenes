# m_matschiner Mon Jul 2 16:31:17 CEST 2018

# Finish the alignments. This means:
# - make a strict-consensus sequences for gadMor_Stat and melAeg
# - use only sites for which no variation is found between the blast/lastz-based approaches and the mapping-based approach.
# - remove sites that are gaps in the first sequence (which is assumed to be the reference).
for refined_alignment in ../res/alignments/merged/lg??.threaded.refined.fasta
do
    # Set variables.
    finished_alignment=${refined_alignment%.fasta}.finished.fasta
    chromosome_id=`basename ${refined_alignment%.threaded.refined.fasta}`
    out=../log/misc/finish_alignment.${chromosome_id}.out
    mapping_alignment=../res/fasta/converted_from_bam/${chromosome_id}.fasta

    # Prepare an alignment for this chromosome, comprising the three sequences based on mapping.
    rm -f ${mapping_alignment}
    echo ">gadMor2_mapping" > ${mapping_alignment}
    cat ../res/fasta/converted_from_bam/gadMor2.fasta | grep -i -A 1 ${chromosome_id} | tail -n 1 >> ${mapping_alignment}
    echo ">gadMor_Stat_mapping" >> ${mapping_alignment}
    cat ../res/fasta/converted_from_bam/gadMor_Stat.fasta | grep -i -A 1 ${chromosome_id} | tail -n 1 >> ${mapping_alignment}
    echo ">melAeg_mapping" >> ${mapping_alignment}
    cat ../res/fasta/converted_from_bam/melAeg.fasta | grep -i -A 1 ${chromosome_id} | tail -n 1 >> ${mapping_alignment}

    # Finish the alignment.
    if [ ! -f ${finished_alignment} ]
    then
        rm -f ${out}
        sbatch -o ${out} finish_alignment.slurm ${refined_alignment} ${mapping_alignment} ${finished_alignment}
    fi

done
