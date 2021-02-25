# m_matschiner Tue Jun 26 18:07:33 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Make the output directories if they don't exist yet.
mkdir -p ../res/lastz/pairwise_alignments/gadMor_Stat
mkdir -p ../res/lastz/pairwise_alignments/melAeg

# Extract sequences from the pairwise alignments made with masa_cudalign, without added ns.
for input_alignment_file in ../res/alignments/masa_cudalign/pairwise_without_ns/gadMor_Stat/lg??/alignment.00.txt
do
    # Set variables for ids and names.
    chromosome_id=`basename ${input_alignment_file%/alignment.00.txt}`
    output_alignment_file=${input_alignment_file%/alignment.00.txt}.fasta
    reference_seq_file=../data/assemblies/gadMor2/gadMor2_${chromosome_id}.fasta
    sequence_id=gadMor_Stat_${chromosome_id}_pairwise_without_ns

    # Extract sequences.
    if [ ! -f ${output_alignment_file} ]
    then
	   ruby convert_to_fasta.rb ${input_alignment_file} ${output_alignment_file} ${reference_seq_file} ${chromosome_id} ${sequence_id}
    fi
done
for input_alignment_file in ../res/alignments/masa_cudalign/pairwise_without_ns/melAeg/lg??/alignment.00.txt
do
    # Set variables for ids and names.
    chromosome_id=`basename ${input_alignment_file%/alignment.00.txt}`
    output_alignment_file=${input_alignment_file%/alignment.00.txt}.fasta
    reference_seq_file=../data/assemblies/gadMor2/gadMor2_${chromosome_id}.fasta
    sequence_id=melAeg_${chromosome_id}_pairwise_without_ns

    # Extract sequences.
    if [ ! -f ${output_alignment_file} ]
    then
           ruby convert_to_fasta.rb ${input_alignment_file} ${output_alignment_file} ${reference_seq_file} ${chromosome_id} ${sequence_id}
    fi
done

# Extract sequences from the pairwise alignments made with masa_cudalign, with added ns.
for input_alignment_file in ../res/alignments/masa_cudalign/pairwise_with_ns/gadMor_Stat/lg??/alignment.00.txt
do
    # Set variables for ids and names.
    chromosome_id=`basename ${input_alignment_file%/alignment.00.txt}`
    output_alignment_file=${input_alignment_file%/alignment.00.txt}.fasta
    reference_seq_file=../data/assemblies/gadMor2/gadMor2_${chromosome_id}.fasta
    sequence_id=gadMor_Stat_${chromosome_id}_pairwise_with_ns

    # Extract sequences.
    if [ ! -f ${output_alignment_file} ]
    then
        ruby convert_to_fasta.rb ${input_alignment_file} ${output_alignment_file} ${reference_seq_file} ${chromosome_id} ${sequence_id}
    fi
done
for input_alignment_file in ../res/alignments/masa_cudalign/pairwise_with_ns/melAeg/lg??/alignment.00.txt
do
    # Set variables for ids and names.
    chromosome_id=`basename ${input_alignment_file%/alignment.00.txt}`
    output_alignment_file=${input_alignment_file%/alignment.00.txt}.fasta
    reference_seq_file=../data/assemblies/gadMor2/gadMor2_${chromosome_id}.fasta
    sequence_id=melAeg_${chromosome_id}_pairwise_with_ns

    # Extract sequences.
    if [ ! -f ${output_alignment_file} ]
    then
        ruby convert_to_fasta.rb ${input_alignment_file} ${output_alignment_file} ${reference_seq_file} ${chromosome_id} ${sequence_id}
    fi
done

# Extract sequences from the pairwise alignments made with lastz.
input_alignment_file=../res/lastz/pairwise_alignments/n_c_ms1_sorted.maf
for number in `seq -w 1 23`
do
    # Set variables for ids and names.
    chromosome_id=lg${number}
    output_alignment_file=../res/lastz/pairwise_alignments/gadMor_Stat/${chromosome_id}.fasta
    reference_seq_file=../data/assemblies/gadMor2/gadMor2_${chromosome_id}.fasta
    sequence_id=gadMor_Stat_${chromosome_id}_lastz

    # Extract sequences.
    if [ ! -f ${output_alignment_file} ]
    then
	ruby convert_to_fasta.rb ${input_alignment_file} ${output_alignment_file} ${reference_seq_file} ${chromosome_id} ${sequence_id}
    fi
done
input_alignment_file=../res/lastz/pairwise_alignments/n_h_ms1_sorted.maf
for number in `seq -w 1 23`
do
    # Set variables for ids and names.
    chromosome_id=lg${number}
    output_alignment_file=../res/lastz/pairwise_alignments/melAeg/${chromosome_id}.fasta
    reference_seq_file=../data/assemblies/gadMor2/gadMor2_${chromosome_id}.fasta
    sequence_id=melAeg_${chromosome_id}_lastz

    # Extract sequences.
    if [ ! -f ${output_alignment_file} ]
    then
        ruby convert_to_fasta.rb ${input_alignment_file} ${output_alignment_file} ${reference_seq_file} ${chromosome_id} ${sequence_id}
    fi
done
