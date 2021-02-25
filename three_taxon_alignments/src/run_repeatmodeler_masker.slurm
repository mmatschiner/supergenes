#!/bin/sh

#SBATCH --job-name=repeatmodeler
#SBATCH --account=nn9244k
#SBATCH --time=168:0:0
#SBATCH --mem-per-cpu=3900M 
#SBATCH --cpus-per-task=16

## Set up the job environment.
source /cluster/bin/jobsetup
module load repeatmodeler/1.0.8

mkdir -p ../data/repeat_libraries
cd ../data/repeat_libraries
mkdir -p $1
cd $1

BuildDatabase -name ${1} -engine ncbi ${2}

RepeatModeler -engine ncbi -pa 16 -database ${1} 1> repeatmodeler.out 2> repeatmodeler.err

RepeatMasker -lib RM*/consensi.fa.classified ${2} -xsmall -pa 16 1> repeatmasker.out 2> repeatmasker.err
