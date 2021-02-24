# michaelm Tue Jan 24 17:57:13 CET 2017

# Load the plink module.
module load plink/1.90b3b

# Make the results directory.
mkdir -p ../res/plink

# Use plink to calculate ld for the entire chromsome (using the ped files produced with strict filters).
for i in ../data/ped/LG??_strict.ped
do
    id=`basename ${i%.ped}`
    plink --threads 1 --file ../res/vcftools/${id} --aec --r2 --ld-window-kb 50000 --ld-window 50000000 --ld-window-r2 0.5 --out ../res/plink/${id}
done
