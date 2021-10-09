# m_matschiner Thu Sep 9 17:42:31 CEST 2021

# Make the results directory.
mkdir -p ../res/tables

# Make the log directory.
mkdir -p ../log/misc

# Set the account.
acct=nn9883k

# Set the window size.
win_size=250000

# Repeat for each vcf file.
for gzvcf in ../data/variants/*.vcf.gz ../res/gatk/gadMor3/original_gadoga/chr{1,2,7,12}.filtered.vcf.gz
do
    # Get the lg id.
    lg_id=`basename ${gzvcf%.filtered.vcf.gz}`

    # Set the population groups based on linkage group.
    if [ ${lg_id} == "LG01" ] || [ ${lg_id} == "chr1" ]
    then
	pop_group1="Gadmor_avo1,Gadmor_avo2,Gadmor_ico1,Gadmor_ico2,Gadmor_lfo1,Gadmor_lfo2,Gadmor_two1,Gadmor_two2"
	pop_group2="Gadmor_avc1,Gadmor_avc2,Gadmor_bat1,Gadmor_bat2,Gadmor_bor1,Gadmor_bor2,Gadmor_icc1,Gadmor_kie1,Gadmor_kie2,Gadmor_lfc1,Gadmor_lfc2,Gadmor_lfc3,Gadmor_low1,Gadmor_low2,Gadmor_twc1,Gadmor_twc2"
    elif [ ${lg_id} == "LG02" ] || [ ${lg_id} == "chr2" ]
    then
	pop_group1="Gadmor_avc1,Gadmor_avc2,Gadmor_kie1,Gadmor_kie2,Gadmor_lfc1,Gadmor_lfc2,Gadmor_lfc3,Gadmor_low1,Gadmor_low2"
	pop_group2="Gadmor_avo1,Gadmor_avo2,Gadmor_bat1,Gadmor_bat2,Gadmor_bor1,Gadmor_bor2,Gadmor_icc1,Gadmor_ico1,Gadmor_ico2,Gadmor_lfo1,Gadmor_lfo2,Gadmor_twc1,Gadmor_twc2,Gadmor_two1,Gadmor_two2"
    elif [ ${lg_id} == "LG07" ] || [ ${lg_id} == "chr7" ]
    then
	pop_group1="Gadmor_avc1,Gadmor_avc2,Gadmor_bor1,Gadmor_bor2,Gadmor_lfc1,Gadmor_lfc2,Gadmor_lfc3,Gadmor_low1,Gadmor_low2,Gadmor_kie1,Gadmor_kie2"
	pop_group2="Gadmor_avo1,Gadmor_avo2,Gadmor_bat1,Gadmor_bat2,Gadmor_icc1,Gadmor_ico1,Gadmor_ico2,Gadmor_lfo1,Gadmor_lfo2,Gadmor_twc1,Gadmor_twc2,Gadmor_two1,Gadmor_two2"
    elif [ ${lg_id} == "LG12" ] || [ ${lg_id} == "chr12" ]
    then
	pop_group1="Gadmor_lfc3,Gadmor_kie1,Gadmor_kie2,Gadmor_low1,Gadmor_low2"
	pop_group2="Gadmor_avc1,Gadmor_avc2,Gadmor_avo1,Gadmor_avo2,Gadmor_bat1,Gadmor_bat2,Gadmor_bor1,Gadmor_bor2,Gadmor_icc1,Gadmor_ico1,Gadmor_ico2,Gadmor_lfo1,Gadmor_lfo2,Gadmor_twc1,Gadmor_twc2,Gadmor_two1,Gadmor_two2"
    else
	echo "ERROR: Linkage group id (${lg_id}) could not be identified!"
	exit 1
    fi
    
    # Submit only if the result doesn't exist yet.
    if [ ! -f ../res/tables/${lg_id}_stats.txt ]
    then

	# Set the name of the log file.
	log=../log/misc/get_window_stats.${lg_id}.txt
	rm -f ${log}

	# Submit the job to calculate stats in sliding windows.
	sbatch --account ${acct} -o ${log} get_window_stats.slurm ${gzvcf} ${lg_id} ${pop_group1} ${pop_group2} ../res/tables/${lg_id}_stats.txt ${win_size}
    fi
done
