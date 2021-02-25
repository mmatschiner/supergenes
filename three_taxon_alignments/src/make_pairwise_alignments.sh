# m_matschiner Thu Jun 15 14:23:39 CEST 2017

# Make the output and log directories if they don't exist yet.
mkdir -p ../res/alignments/masa_cudalign/pairwise_with_ns/gadMor_Stat
mkdir -p ../res/alignments/masa_cudalign/pairwise_with_ns/melAeg
mkdir -p ../res/alignments/masa_cudalign/pairwise_without_ns/gadMor_Stat
mkdir -p ../res/alignments/masa_cudalign/pairwise_without_ns/melAeg 
mkdir -p ../log/masa_cudalign/

# Pairwise align sequences with masa-cudalign.
for specimen in gadMor_Stat melAeg
do
    for i in ../res/fasta/merged_contigs_per_lg_with_ns/${specimen}/*.fasta
    do
	    lg=`basename ${i%.fasta}`
	    log_file="../log/masa_cudalign/${specimen}_${lg}_with_ns.out"
	    res_dir="../res/alignments/masa_cudalign/pairwise_with_ns/${specimen}/${lg}"
	    if [ ! -f ${res_dir}/alignment.00.txt ]
	    then
	        rm -f ${log_file}
	        sbatch -o ${log_file} run_masa_cudalign.slurm ${i} ../data/assemblies/gadMor2/gadMor2_${lg}.fasta ${res_dir}
	    fi
    done
    for i in ../res/fasta/merged_contigs_per_lg_without_ns/${specimen}/*.fasta
    do
	    lg=`basename ${i%.fasta}`
	    log_file="../log/masa_cudalign/${specimen}_${lg}_without_ns.out"
	    rm -f ${log_file}
	    sbatch -o ${log_file} run_masa_cudalign.slurm ${i} ../data/assemblies/gadMor2/gadMor2_${lg}.fasta ../res/alignments/masa_cudalign/pairwise_without_ns/${specimen}/${lg}
    done
done
