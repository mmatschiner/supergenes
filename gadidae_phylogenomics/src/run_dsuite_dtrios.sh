# michaelm Sat May 30 16:16:40 CEST 2020

# Download dsuite.
mkdir -p ../bin
if [ ! -f ../bin/Dsuite ]
then
    module load GCC/8.3.0
    git clone https://github.com/millanek/Dsuite.git
    cd Dsuite
    make
    cd -
    mv Dsuite/Build/Dsuite ../bin
    rm -rf Dsuite
fi

# Set the list of samples and the tree.
samples=../data/tables/samples.txt
tree=../data/trees/aim_species_simple_ladder_1_nwk.tre

# Run dsuite separately for all snps, snps outside of inversion regions, and snps in inversion regions.
for vcfgz in ../res/gatk/concatenated.masked.variable*.vcf.gz
do
    bin_id=`basename ${vcfgz} | cut -d "." -f 4`
    if [ ${bin_id} == "vcf" ]
    then
	    bin_id="full"
    fi

    # Make the result directory.
    mkdir -p ../res/dsuite/${bin_id}

    # Feedback.
    echo "Running dsuite for the ${bin_id} set of snps."

    # Run dsuite.
    ../bin/Dsuite Dtrios -t ${tree}  ${vcfgz} ${samples}

    # Move result files to the result directory.
    mv ${samples%.txt}__* ../res/dsuite/${bin_id}
done
