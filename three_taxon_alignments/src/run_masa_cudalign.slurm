#!/bin/bash

# m_matschiner Thu Aug 17 00:05:53 CEST 2017

# This script runs the software MASA-cudalign.
#
# Start with e.g.
# 'sbatch run_masa_cudalign.slurm infile1.fasta infile2.fasta outfile.fasta'

# Job name:
#SBATCH --job-name=mscdln
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
# Project.
#SBATCH --account=nn9244k
#
# Processor and memory usage.
#SBATCH --partition=accel --gres=gpu:1
#SBATCH --mem-per-cpu=60G
#SBATCH --cpus-per-task=1

# Set up the job environment.
source /cluster/bin/jobsetup

# Load the cuda module.
module load cuda/8.0

# Get command line arguments.
infile1_w_rel_path=$1
infile2_w_rel_path=$2
outdir_w_rel_path=$3
infile1=`basename ${infile1_w_rel_path}`
infile2=`basename ${infile2_w_rel_path}`
outdir_w_abs_path=`readlink -f ${outdir_w_rel_path}`

# Copy the input files to the analysis directory.
cp ${infile1_w_rel_path} $SCRATCH
cp ${infile2_w_rel_path} $SCRATCH
 
# Move to the analysis directory.
cd $SCRATCH

# Run masa-cudalign.
/projects/cees/bin/MASA-CUDAlign/3.9.1.1024/cudalign -c --ram-size=55G --alignment-start=2 --alignment-end=2 ${infile1} ${infile2}

# Make the output directory.
mkdir -p ${outdir_w_abs_path}

# Copy the results to the output directory.
cp -r ./work.tmp/* ${outdir_w_abs_path}
