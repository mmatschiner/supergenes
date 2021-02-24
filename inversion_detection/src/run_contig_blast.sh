# michaelm Sun Jan 29 17:53:42 CET 2017

# Load the BLAST module.
module load blast+/2.2.29

# Make the output directories if they don't exist yet.
mkdir -p ../res/blast/gadMor_Stat_vs_gadMor2
mkdir -p ../res/blast/melAeg_vs_gadMor2

# For each of the contigs identified manually as potentially crossing inversion boundaries, run BLAST against the respective linkage group.
cat ../res/manual/tables/candidate_contigs.txt | grep -v "#" > tmp.txt
while read i
do
    taxon=`echo $i | cut -d " " -f 1`
    scf_id=`echo $i | cut -d " " -f 2`
    gadMor2_lg=`echo $i | cut -d " " -f 3`
    if [ $taxon == "gadMor_Stat" ]
    then
	    taxon_fasta="../data/assemblies/gadMor_Stat.fasta"
    elif [ $taxon == "melAeg" ]
    then
	    taxon_fasta="../data/assemblies/melAeg.contigs.fasta"
    fi
    fastagrep -t -p ${scf_id} ${taxon_fasta} > ../data/assemblies/${taxon}/${scf_id}.fasta
    blastn -query ../data/assemblies/${taxon}/${scf_id}.fasta -db ../data/assemblies/gadMor2_${gadMor2_lg}.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq" -task blastn -evalue 1e-10 -out ../res/blast/${taxon}_vs_gadMor2/${scf_id}.out -reward 1 -penalty -2 -gapopen 2 -gapextend 1 -num_threads 4 -max_target_seqs 10000
done < tmp.txt

# Cleanup.
rm tmp.txt
