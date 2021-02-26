# m_matschiner Sun Feb 7 10:44:47 CET 2021

# Download the gadMor2 assembly.
if [ -f ../data/assemblies/gadMor2.fasta ]
cd ../data/assemblies/
    wget https://ndownloader.figshare.com/files/5323414
    mv 5323414 gadMor2.fasta.gz
    gunzip gadMor2.fasta.gz
    cd -
fi

# Set the alternative reference fasta.
melaeg_fasta=../data/assemblies/melAeg_in_gadMor2_coords.fasta

# Copy melAeg sequences from the three-way alignment into a new fasta file.
if [ ! -f ${melaeg_fasta} ]
then
    touch ${melaeg_fasta}
    for n_padded in `seq -w 23`
    do
        echo ">LG${n_padded}" >> ${melaeg_fasta}
        cat ../../three_taxon_alignments/res*/alignments/merged/lg${n_padded}.threaded.refined.finished.fasta | grep -A 2 haddock | tail -n 1 | tr "-" "N" >> ${melaeg_fasta}
    done
fi

# Generate index files for the references.
for fasta in ../data/assemblies/gadMor2.fasta ${melaeg_fasta}
do
    fasta_id=`basename ${fasta%.fasta}`
    out=../log/misc/prepare_${fasta_id}.out
    log=../log/misc/prepare_${fasta_id}.log
    rm -f ${out}
    rm -f ${log}
    sbatch -o ${out} prepare.slurm ${fasta} ${log}
done
