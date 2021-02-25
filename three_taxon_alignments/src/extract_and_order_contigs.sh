# m_matschiner Wed Jun 14 15:53:49 CEST 2017

# Extract and order contigs for all lgs and for both the gadMor_Stat and the melAeg assemblies.
for specimen in melAeg gadMor_Stat
do
    for n in `seq -w 23`
    do
    	lg_id="LG${n}"
		bash extract_and_order_contigs_per_lg.sh ${specimen} ${lg_id}
    done
done