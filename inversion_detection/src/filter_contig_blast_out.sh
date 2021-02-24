# michaelm Wed Jan 18 15:21:34 CET 2017

# Load modules.
module load ruby/2.1.5

# Remove previous results.
rm -f ../res/blast/gadMor_Stat_vs_gadMor2/*filtered.out
rm -f ../res/blast/melAeg_vs_gadMor2/*filtered.out

# Filter the output of BLAST searches with gadMor_Stat scaffolds as queries.
for i in ../res/blast/gadMor_Stat_vs_gadMor2/scf*.out
do
    ruby filter_blast_out.rb ${i}
    rm $i
done

# Filter the output of BLAST searches with melAeg scaffolds as queries.
for i in ../res/blast/melAeg_vs_gadMor2/MeA*.out
do
    ruby filter_blast_out.rb ${i}
    rm $i
done
