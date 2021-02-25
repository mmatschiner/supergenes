#!/bin/bash

# m_matschiner Fri Jun 29 23:38:42 CEST 2018

# Job name:
#SBATCH --job-name=rfnthr
#
# Wall clock limit:
#SBATCH --time=24:00:00
#
# Processor and memory usage:
#SBATCH --mem-per-cpu=30G

## Set up the job environment.
source /cluster/bin/jobsetup
module load ruby/2.1.5
module load mafft/7.300

## Feedback.
echo
echo "refine_threaded_alignment.slurm was called with arguments ${@}."
echo

# Get the command-line arguments.
threaded_alignment_with_relative_path=${1}
threaded_alignment_with_absolute_path=`readlink -f ${threaded_alignment_with_relative_path}`
refined_threaded_alignment_with_relative_path=${2}
refined_threaded_alignment_with_absolute_path=`readlink -f ${refined_threaded_alignment_with_relative_path}`
cp refine_threaded_alignment.rb ${SCRATCH}

# Change to the temporary directory.
cd ${SCRATCH}

# Copy the threaded alignment and remove the duplicate reference sequences.
cat ${threaded_alignment_with_absolute_path} | head -n 2 > tmp.fasta
cat ${threaded_alignment_with_absolute_path} | grep --no-group-separator -e gadMor_Stat -e melAeg -A 1 >> tmp.fasta

# Refine the threaded alignments by local realignment with mafft.
ruby refine_threaded_alignment.rb tmp.fasta ${refined_threaded_alignment_with_absolute_path}
