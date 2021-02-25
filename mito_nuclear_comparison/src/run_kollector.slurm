#!/bin/bash

# m_matschiner Thu Oct 24 11:13:42 CEST 2019

# Job name:
#SBATCH --job-name=kollct
#
# Wall clock limit:
#SBATCH --time=100:00:00
#
# Processor and memory usage:
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=30G

# Define a function to print the time.
function print_time {
    END=$(date +%s)
    DIFF=$(( $END - $1 ))
    dd=$(echo "$DIFF/86400" | bc)
    dt2=$(echo "$DIFF-86400*$dd" | bc)
    dh=$(echo "$dt2/3600" | bc)
    dt3=$(echo "$dt2-3600*$dh" | bc)
    dm=$(echo "$dt3/60" | bc)
    ds=$(echo "$dt3-60*$dm" | bc)
    if [ $dd -gt 0 ]; then
	echo "done in ${dd} days and ${dh} hours."
    elif [ $dh -gt 0 ]; then
	echo "done in ${dh} hours and ${dm} minutes."
    elif [ $dm -gt 0 ]; then
	echo "done in ${dm} minutes and ${ds} seconds."
    else
	echo "done in ${ds} seconds."
    fi
}

## Set up the job environment.
source /cluster/bin/jobsetup
module load kollector/1.0.1

## Feedback.
echo
echo "run_kollector.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
reads1=`readlink -f ${1}`
reads2=`readlink -f ${2}`
query=`readlink -f ${3}`
res=`readlink -f ${4}`
log=`readlink -f ${5}`

# Change to the temporary directory.
cd ${SCRATCH}

# Run kollector.
START=$(date +%s)
echo -n "Running kollector... "
kollector_multiple.sh -m6 -j10 -r 0.5 -k 32 -K 20 -n 1000000 ${query} ${reads1} ${reads2} &> ${log}
print_time ${START}

# Copy results to the output directory.
cp iteration.6/kollector_assembledtargets.fa ${res}