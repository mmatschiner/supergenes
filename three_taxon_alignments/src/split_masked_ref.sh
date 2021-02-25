# m_matschiner Wed Apr 4 11:49:41 CEST 2018

# Set the masked reference assembly.
ref="../data/assemblies/gadMor2/gadMor2.fasta.masked"

# Split the masked reference into separate files per lg.
for n in `seq -w 23`
do
	lg_id="LG${n}"
    lg_upper=`echo ${lg} | tr a-z A-Z`
    ./fastagrep -t -p ${lg_upper} ${ref} > ../data/assemblies/gadMor2/gadMor2_${lg}.mask.rm.fasta
done