# m_matschiner Tue May 1 16:31:41 CEST 2018

# Define a function to sleep if too many jobs are queued or running.
function sleep_while_too_busy {
    n_jobs=`squeue -u michaelm | wc -l`
    have_been_waiting=false
    while [ $n_jobs -gt 300 ]
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

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/kollector/mitochondrial
mkdir -p ../res/kollector/nuclear

# Make the log directory if it doesn't exist yet.
mkdir -p ../log/kollector/mitochondrial
mkdir -p ../log/kollector/nuclear

# Set the account.
acct=nn9244k

# Start a slurm script per specimen.
for query in ../data/queries/mitochondrial_gadmor2.fasta ../data/queries/gadmor_exons_nucl.fasta
do
    if [ ${query} == ../data/queries/gadmor_exons_nucl.fasta ]
    then
	query_type=nuclear
    else
	query_type=mitochondrial
    fi
    for reads2 in `ls ../data/reads/*R2.fastq.gz`
    do
	reads1=${reads2%R2.fastq.gz}R1.fastq.gz
	species_id=`basename ${reads2} | cut -d "_" -f 1`
	specimen_id=`basename ${reads2} | cut -d "_" -f 2`
	out=../log/kollector/${query_type}/${species_id}_${specimen_id}.out
	log=../log/kollector/${query_type}/${species_id}_${specimen_id}.log
	res=../res/kollector/${query_type}/${species_id}_${specimen_id}.kollector.fasta
	if [ ! -f ${res} ]
	then
	    rm -f ${out}
	    rm -f ${log}
	    echo -n "Submitting job for specimen ${specimen_id}..."
	    sbatch --account ${acct} -o ${out} run_kollector.slurm ${reads1} ${reads2} ${query} ${res} ${log} &> /dev/null
	    echo " done."
	    sleep_while_too_busy
	fi
    done
done