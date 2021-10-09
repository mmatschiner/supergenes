# m_matschiner Thu Sep 9 16:52:43 CEST 2021

# Make the data directory.
mkdir -p ../data/assemblies

# Download each chromosome from ena.
chr=1
for i in {43..65}
do
    echo -n "Downloading chromosome ${chr} sequence..."
    wget https://www.ebi.ac.uk/ena/browser/api/fasta/LR6339${i}.1 &> /dev/null
    echo ">chr${chr}" > tmp.fasta
    tail -n +2 LR6339${i}.1 >> tmp.fasta
    mv -f tmp.fasta LR6339${i}.1
    chr=$((chr+1))
    echo " done."
done

# Concatenate the chromosomes into a single file.
cat LR6339??.1 > ../data/assemblies/gadMor3.fasta

# Clean up.
rm -f LR6339??.1

# Set the account.
acct=nn9244k

# Generate index files for the gadmor3 reference.
fasta_id=gadMor3
out=../log/misc/prepare_${fasta_id}.out
log=../log/misc/prepare_${fasta_id}.log
rm -f ${out}
rm -f ${log}
sbatch --account ${acct} -o ${out} prepare.slurm ../data/assemblies/gadMor3.fasta ${log}
