# m_matschiner Mon Jul 16 00:16:26 CEST 2018

# Keep only the first 23 sequences (without mtdna and unplaced).
for assembly_id in gadMor2 gadMor_Stat melAeg
do
	cat ../res/fasta/converted_from_bam/${assembly_id}.fasta | head -n 46 > tmp.fasta
	mv -f tmp.fasta ../res/fasta/converted_from_bam/${assembly_id}.fasta
done