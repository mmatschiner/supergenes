# m_matschiner Fri Jul 20 11:43:20 CEST 2018

# Load the blast module.
module load blast+/2.2.29

# Make a blast database for each fasta subject.
for i in ../data/subjects/*.fasta ../res/kollector/*.fasta
do
    if [ ! -f ${i}*.nhr ]
    then
        makeblastdb -in ${i} -dbtype nucl
    fi
done