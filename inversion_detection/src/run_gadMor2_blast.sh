# michaelm Mon Jan 16 18:15:37 CET 2017

# Load the BLAST module.
module load blast+/2.2.29

# Make result directories.
mkdir -p ../res/blast/gadMor2_vs_gadMor_Stat
mkdir -p ../res/blast/gadMor2_vs_melAeg

# For each linkage group, run BLAST searches in parallel.
for lg in `seq -w 23`
do

	# Run a BLAST search for each target region, with the gadMor_Stat assembly as database.
	for i in ../data/assemblies/gadMor2_lg${lg}_*.fasta;
	do
	    target_region_id_raw=`basename ${i}`
	    target_region_id_raw2=${target_region_id_raw%.fasta}
	    target_region_id=`echo ${target_region_id_raw2} | sed 's/gadMor2_//g'`
	    blastn -query ${i} -db ../data/assemblies/gadMor_Stat.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq" -task blastn -evalue 1e-10 -out ../res/blast/neac_vs_coastal/${target_region_id}.out -reward 1 -penalty -2 -gapopen 2 -gapextend 1 -num_threads 1 -max_target_seqs 10000 &
	done
	wait

	# Run a BLAST search for each target region, with the haddock assembly as database.
	for i in ../data/assemblies/gadMor2_lg${lg}_*.fasta;
	do
	    target_region_id_raw=`basename ${i}`
	    target_region_id_raw2=${target_region_id_raw%.fasta}
	    target_region_id=`echo ${target_region_id_raw2} | sed 's/gadMor2_//g'`
	    blastn -query ${i} -db ../data/assemblies/melAeg.contigs.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq" -task blastn -evalue 1e-10 -out ../res/blast/neac_vs_haddock/${target_region_id}.out -reward 1 -penalty -2 -gapopen 2 -gapextend 1 -num_threads 1 -max_target_seqs 10000 &
	done
	wait

done
