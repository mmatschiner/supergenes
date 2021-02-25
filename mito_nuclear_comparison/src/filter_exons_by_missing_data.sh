# m_matschiner Tue Nov 5 11:54:10 CET 2019

# Make the result directory.
mkdir -p ../res/alignments/nuclear/03

# Copy alignments to the 03 directory if they contain a sequence for the low-coverage arcgla sample.
for align in ../res/alignments/nuclear/02/*_nucl.fasta
do
    align_id=`basename ${align%_nucl.fasta}`
    arcgla1_length=`cat ${align} | grep -A 1 SRR2906188 | tail -n 1 | tr -d "-" | wc -m`
    arcgla1_length=$((arcgla1_length-1))
    seqs_missing=`cat ${align} | grep -v ">" | grep -v A | grep -v C | grep -v G | grep -v T | wc -l`
    if (( ${seqs_missing} < 6 ))
    then
        echo "Copying alignment ${align_id} with ${seqs_missing} missing sequences and ${arcgla1_length} bp for arcgla."
        cp ${align} ../res/alignments/nuclear/03
    elif (( ${seqs_missing} < 11 ))
    then
        if (( ${arcgla1_length} > 0 ))
        then
            echo "Copying alignment ${align_id} with ${seqs_missing} sequences and ${arcgla1_length} bp for arcgla."
            cp ${align} ../res/alignments/nuclear/03
        else
            echo "Excluding alignment ${align_id} with ${seqs_missing} missing sequences and ${arcgla1_length} bp for arcgla."
        fi
    else
        echo "Excluding alignment ${align_id} with ${seqs_missing} missing sequences and ${arcgla1_length} bp for arcgla."
    fi
done