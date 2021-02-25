# m_matschiner Tue Nov 5 09:38:17 CET 2019

# Load the blast module.
module load blast+/2.2.29

# Make a blast database for each local assembly of nuclear genes made with kollector.
for i in ../data/subjects/*.fasta ../res/kollector/nuclear/*.fasta
do
    if [ ! -f ${i}*.nhr ]
    then
        makeblastdb -in ${i} -dbtype nucl
    fi
done