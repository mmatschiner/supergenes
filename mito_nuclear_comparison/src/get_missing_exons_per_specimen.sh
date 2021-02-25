# m_matschiner Tue Nov 5 13:26:25 CET 2019

# Get the specimen ids.
first_align=`ls ../res/alignments/nuclear/03/*.fasta | head -n 1`
specimens=`cat ${first_align} | grep ">" | tr -d ">"`

# Get the number of sequences for which each specimen does have some data.
for specimen in ${specimens}
do
    seq_count=0
    for align in ../res/alignments/nuclear/03/*.fasta
    do
        seq_length=`cat ${align} | grep -A 1 ${specimen} | tail -n 1 | tr -d "-" | wc -m`
        if [[ ${seq_length} > 1 ]]
        then
            seq_count=$((seq_count+1))
        fi
    done
    echo -e "${specimen}\t${seq_count}"
done