# m_matschiner Wed Jun 14 14:11:42 CEST 2017

# Make the output directories if they don't exist yet.
mkdir -p ../res/fasta/separate_contigs_per_lg/gadMor_Stat
mkdir -p ../res/fasta/separate_contigs_per_lg/melAeg

# Load the ruby module
module load ruby/2.1.5

# Get the command-line arguments.
specimen=${1}
lg=${2}

# Get the name of the reference file.
if [ ${specimen} == "gadMor_Stat" ]
then
    ref="../data/assemblies/gadMor_Stat/gadMor_Stat.fasta"
elif [ ${specimen} == "melAeg" ]
then
    ref="../data/assemblies/melAeg/melAeg.fasta"
else
    echo "ERROR: Unexpected specimen specified: ${specimen}!"
    exit 1
fi

# Clean up results of previous analyses.
rm -f ../res/fasta/separate_contigs_per_lg/${specimen}/${lg}.fasta

# Extract and order contig sequences.
infile_name="../data/tables/${specimen}_contig_order_${lg}.txt"
for i in `cat ${infile_name} | sed 's/\t/,/g'`
do
    contig_id=`echo ${i} | cut -d "," -f 1`
    orientation=`echo ${i} | cut -d "," -f 2`
    if [ ${orientation} == "for" ]
    then
	    ./fastagrep -t -p ${contig_id} ${ref} >> ../res/fasta/separate_contigs_per_lg/${specimen}/${lg}.fasta
    elif [ ${orientation} == "rev" ]
    then
	    ./fastagrep -t -p ${contig_id} ${ref} > tmp_${specimen}_${lg}.fasta
	    ruby make_rev_comp.rb tmp_${specimen}_${lg}.fasta >> ../res/fasta/separate_contigs_per_lg/${specimen}/${lg}.fasta
	    rm tmp_${specimen}_${lg}.fasta
    else
	    echo "ERROR: Found unexpected contig orientation: ${orientation}"
	exit 1
    fi
done
