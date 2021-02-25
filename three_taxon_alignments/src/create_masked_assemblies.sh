# ole

#Creating masked assemblies for better/easier whole genome alignment.
#First creating a RepeatModeler library, then masking the genome with RepeatMasker. 
#RepeatModeler is stochastic, so it is not completely reproducable.

# Load modules.
module load repeatmodeler/1.0.8

#create a folder for the repeat libraries and the RepeatModeler run
mkdir -p ../data/repeat_libraries

# Make the log directory.
mkdir -p ../log/repeat_masking

out="../log/repeat_masking/melAeg.txt"
rm -f ${out}
sbatch -o ${out} run_repeatmodeler_masker.slurm melAeg ../data/assemblies/melAeg/melAeg.fasta

out="../log/repeat_masking/gadMor_Stat.txt"
rm -f ${out}
sbatch -o ${out} run_repeatmodeler_masker.slurm gadMor_Stat ../data/assemblies/gadMor_Stat/gadMor_Stat.fasta

out="../log/repeat_masking/gadMor2.txt"
rm -f ${out}
sbatch -o ${out} run_repeatmodeler_masker.slurm gadMor2 ../data/assemblies/gadMor2/gadMor2.fasta

