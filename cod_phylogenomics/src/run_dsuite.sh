# michaelm Mon Feb 21 00:59:50 CET 2021

# Make the bin directory.
mkdir -p ../bin

# Make the dsuite results directory.
mkdir -p ../res/dsuite

# Make the log directory.
mkdir -p ../log/dsuite

# Download dsuite.
if [ ! -f ../bin/Dsuite ]
then
    git clone https://github.com/millanek/Dsuite.git
    cd Dsuite
    make
    cd -
    mv Dsuite/Build/Dsuite ../bin
    rm -rf Dsuite
fi

# Repeat for each of the alternative references.
for ref in ../data/assemblies/gadMor2.fasta ../data/assemblies/melAeg_in_gadMor2_coords.fasta
do
    # Set the reference id.
    ref_id=`basename ${ref%.fasta}`

    # Repeat first with the original gadoga specimen and then with the alternative one.
    for gadoga_setting in original_gadoga alternative_gadoga
    do

        # Run dtrios for all regions.
        for gzvcf in ../res/gatk/${ref_id}/${gadoga_setting}/colinear.full.vcf.gz ../res/gatk/${ref_id}/${gadoga_setting}/inversion_lg??.full.vcf.gz
        do
            # Get the vcf id.
            gzvcf_id=`basename ${gzvcf%.full.vcf.gz}`
            
            # Make the results directory for this dsuite analysis.
            res_dir=../res/dsuite/${ref_id}/${gadoga_setting}/${gzvcf_id}
            mkdir -p ${res_dir}

            # Convert tree files to newick format if necessary.
            module --quiet purge
            module load R/3.6.2-foss-2019b
            tre=../res/snapp/regions/${gzvcf_id}/combined/${gzvcf_id}.tre
            nwk=${res_dir}/${gzvcf_id}_nwk.tre
            if [ ! -f ${nwk} ]
            then
                Rscript convert_nexus_trees_to_newick.r ${tre} tmp.tre
                if [ ${gadoga_setting} == "original_gadoga" ]
                then
                    cat tmp.tre | sed 's/_spc//g' > ${nwk}
                elif [ ${gadoga_setting} == "alternative_gadoga" ]
                then
                    cat tmp.tre | sed 's/_spc//g' | sed 's/Gadoga/Gadog2/g' > ${nwk}
                fi
                rm -f tmp.tre
            fi

            # Run the dtrios function of dsuite.
            out=../log/dsuite/run_dsuite.${ref_id}_${gadoga_setting}_${gzvcf_id}.txt
            rm -f ${out}
            sbatch -o ${out} run_dsuite.slurm ${gzvcf} ${nwk} ${res_dir}
        done
    done
done
