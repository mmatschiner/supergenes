#!/bin/sh

# JMIB [jmi.barth@gmail.com], 2020-06-10
# Using RELATE to produce .anc (trees per region) and .mut files (SNP information)

#SBATCH --job-name=relate
#SBATCH --account=XXX
#SBATCH --time=23:59:00
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=1

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset modules
	module load BCFtools/1.9-foss-2018b
	module load R/3.5.1-foss-2018b
module list    # For easier debugging

## Feedback
echo
echo "run_relate_prep.slurm was called with arguments ${@}."
echo

## Get command line arguments, set variables, and define output directory.
chr=$1 #NC_031965.f5.2.vcf.gz

n=`echo ${chr%_colinear.vcf.gz} | sed 's/LG0//g' | sed 's/LG//g'`
n_padded=`echo ${chr%_colinear.vcf.gz} | sed 's/LG//g'`

map=${1%_colinear.vcf.gz}.flat.map #NC_031965.flat.map
poplabels=samples.poplabels.txt
ancestor=ancestor.colinear.LG${n_padded}.fasta
mask=mask.colinear.LG${n_padded}.fasta
out_dir=`readlink -f ../res/relate/output/colinear_ancestor_anc_mut`
mkdir -p ${out_dir}

## Copy input files to the work directory:
echo -n "Copying infiles..."
	cp ../res/relate/input/${chr} $SCRATCH
	cp ../data/maps/${map} $SCRATCH
	cp ../data/tables/${poplabels} $SCRATCH
	cp ../res/relate/input/${ancestor} $SCRATCH
	cp ../res/relate/input/${mask} $SCRATCH
## Copy relate to work (no recent module file available)
	cp -r ../bin/relate/ $SCRATCH
echo " done."

## Change to temporary directory.
cd $SCRATCH

## Replace chr names (scripts expect chromosomes to be named only with integers) and filter for bialellic SNPs
echo "Replace chr names and filter for biallelic SNPs..."
	echo "gunzip -c ${chr} | sed 's/${chr%_colinear.vcf.gz}/${n}/g' | bcftools view -m2 -M2 -v snps | gzip > ${chr%.vcf.gz}.bial.renamed.vcf.gz"
		  gunzip -c ${chr} | sed 's/${chr%_colinear.vcf.gz}/${n}/g' | bcftools view -m2 -M2 -v snps | gzip > ${chr%.vcf.gz}.bial.renamed.vcf.gz
echo "...done."

## Convert to SHAPEIT's haplotype format
echo "Converting to SHAPEIT's haplotype format..."
	echo "./relate/bin/RelateFileFormats --mode ConvertFromVcf --haps ${chr%.vcf.gz}.bial.renamed.haps --sample ${chr%.vcf.gz}.bial.renamed.sample -i ${chr%.vcf.gz}.bial.renamed"
		  ./relate/bin/RelateFileFormats --mode ConvertFromVcf --haps ${chr%.vcf.gz}.bial.renamed.haps --sample ${chr%.vcf.gz}.bial.renamed.sample -i ${chr%.vcf.gz}.bial.renamed
echo "...done."

## Prep input files: flip haplotypes according to the ancestral genome, adjust distances between SNPs using a genomic mask
echo "Preparing input files..."
	echo "./relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps ${chr%.vcf.gz}.bial.renamed.haps --sample ${chr%.vcf.gz}.bial.renamed.sample --ancestor ${ancestor} --mask ${mask} --poplabels ${poplabels} -o ${chr%.vcf.gz}.bial.renamed.preped"
		  ./relate/scripts/PrepareInputFiles/PrepareInputFiles.sh --haps ${chr%.vcf.gz}.bial.renamed.haps --sample ${chr%.vcf.gz}.bial.renamed.sample --ancestor ${ancestor} --mask ${mask} --poplabels ${poplabels} -o ${chr%.vcf.gz}.bial.renamed.preped
echo "...done."


## Run RELATE to infer genome-wide genealogies and mutations
echo "Running RELATE to infer genome-wide genealogies and mutations..."
	echo "./relate/bin/Relate --mode All -m 1.63922e-8 -N 100000 --seed 1 --haps ${chr%.vcf.gz}.bial.renamed.preped.haps.gz --sample ${chr%.vcf.gz}.bial.renamed.preped.sample.gz --map ${map} -o ${chr%.vcf.gz}.bial.renamed.preped"
		  ./relate/bin/Relate --mode All -m 1.63922e-8 -N 100000 --seed 1 --haps ${chr%.vcf.gz}.bial.renamed.preped.haps.gz --sample ${chr%.vcf.gz}.bial.renamed.preped.sample.gz --map ${map} -o ${chr%.vcf.gz}.bial.renamed.preped
echo "...done."

## Gzip anc and mut output files
echo "Bgzipping output files..."
	bgzip -c ${chr%.vcf.gz}.bial.renamed.preped.mut > ${chr%.vcf.gz}.bial.renamed.preped.mut.gz
		rm -f ${chr%.vcf.gz}.bial.renamed.preped.mut
	bgzip -c ${chr%.vcf.gz}.bial.renamed.preped.anc > ${chr%.vcf.gz}.bial.renamed.preped.anc.gz
		rm -f ${chr%.vcf.gz}.bial.renamed.preped.anc
echo "...done."

## Copy files back
	echo -n "Copying result files back to the submission directory..."
		cp *.bial.renamed.preped.* ${out_dir}
	echo " done."

