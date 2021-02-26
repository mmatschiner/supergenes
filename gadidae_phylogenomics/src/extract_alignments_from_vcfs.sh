# m_matschiner Mon Oct 8 18:02:38 CEST 2018

# Define a function to sleep if too many jobs are queued or running.
function sleep_while_too_busy {
    n_jobs=`squeue -u michaelm | wc -l`
    have_been_waiting=false
    while [ $n_jobs -gt 340 ]
    do
        have_been_waiting=true
        echo -ne "\rWaiting for job capacity..."
        sleep 60
        n_jobs=`squeue -u michaelm | wc -l`
    done
    if [ ${have_been_waiting} == true ]
    then
        echo " done."
    fi
}

# Make the results directories.
for window_size in 5000
do
    mkdir -p ../res/windows/${window_size}bp
done

# Make the log directories.
for window_size in 5000
do
    mkdir -p ../log/windows/${window_size}bp
done

# Set the reference file.
ref=../data/assemblies/gadMor2.fasta

# Set the allowed proportion of missing data and the minimum likelihood difference.
p_missing_allowed=0.8
min_likelihood=20

# Start an analysis for each chromosome and for each window size.
for window_size in 5000
do
    for chromosome_id in LG0{1..9} LG{10..23}
    do

        # Set the compressed vcf file.
        gzvcf=../res/gatk/${chromosome_id}.masked.vcf.gz

        # Set the window start and end positions.
        from=1
        (( to = ${from} + ${window_size} - 1 ))
        
        # Get the chromosome length.
        chr_length=`cat ${ref}.fai | head -n 23 | grep ${chromosome_id} | cut -f 2`

        # Start an analysis for each window on this chromosome.
        while [ ${to} -lt ${chr_length} ]
        do
            # Specify the region.
            region="${chromosome_id}:${from}-${to}"

            # Specify the alignment and info files.
            align=../res/windows/${window_size}bp/${chromosome_id}_${from}_${to}.phy
            info=../res/windows/${window_size}bp/${chromosome_id}_${from}_${to}.info.txt
            #out=../log/windows/${window_size}bp/align.${chromosome_id}_${from}_${to}.out

            # Calculate completeness, and generate an alignment if completeness is sufficient.
            if [ ! -f ${info} ]
            then
            #rm -f ${out}
            echo "Launching alignment extraction for region ${region} (no info file present)."
            bash extract_alignment_from_vcf.sh ${gzvcf} ${region} ${align} ${info} ${p_missing_allowed} ${min_likelihood}
            #sbatch -o ${out} extract_alignment_from_vcf.slurm ${gzvcf} ${region} ${align} ${info} ${p_missing_allowed} ${min_likelihood}
            #sleep_while_too_busy
            else
            n_sites_calculated=`cat ${info} | grep "n_sites:" | wc -l`
            if [[ ${n_sites_calculated} -gt 0 ]]
            then
                n_sites=`cat ${info} | grep "n_sites:" | cut -d ":" -f 2`
                n_sites_required=`echo "${window_size}*${p_missing_allowed}" | bc`
                if (( $(echo "${n_sites} >= ${n_sites_required}" | bc -l ) ))
                then
                if [ ! -f ${align} ]
                then
                    #rm -f ${out}
                    echo "Launching alignment extraction for region ${region} (site number calculated and sufficient)."
                    bash extract_alignment_from_vcf.sh ${gzvcf} ${region} ${align} ${info} ${p_missing_allowed} ${min_likelihood}
                    fi
                fi
            else
                #rm -f ${out}
                echo "Launching alignment extraction for region ${region} (info file present but no site number calculated)."
                bash extract_alignment_from_vcf.sh ${gzvcf} ${region} ${align} ${info} ${p_missing_allowed} ${min_likelihood}
            fi

            fi

            # Shift the window.
            (( from = ${to} + 1 ))
            (( to = ${to} + ${window_size} ))

        done

    done
done
