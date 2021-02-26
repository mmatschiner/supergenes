# m_matschiner Thu Jun 21 12:15:31 CEST 2018

# Purge and load modules.
module --quiet purge
module load Boost/1.71.0-GCC-8.3.0
module load Ruby/2.7.1-GCCcore-8.3.0
module load BCFtools/1.10.2-GCC-8.3.0

# Download the evo program.
rm -rf evo
if [ ! -f ../bin/evo ]
then
    git clone https://github.com/millanek/evo.git
    cd evo
    make
    cd -
    mv evo/Build/evo ../bin
    rm -r evo
fi

# Make the result directory.
mkdir -p ../res/evo

# Set the divergence time.
divergence_time=0.6394 # (= 0.0654 + 10 * 0.0574; the population divergence time between canadian and other cod pops in the colinear snapp analysis plus the assumed generation time in years times the pop size inferred in the same snapp analysis; this is the expected genetic divergence time.)

# Get the size of the genome excluding lgs 01, 02, 07, 12, and masked sites..
mask=../res/relate/input/mask.fasta
for lg in LG03 LG04 LG05 LG06 LG08 LG09 LG10 LG11 LG13 LG14 LG15 LG16 LG17 LG18 LG19 LG20 LG21 LG22 LG23
do
    ./fastagrep -t -p ${lg} ../res/relate/input/mask.fasta | grep -v ">"
done | awk -F"P" '{print NF-1}' | awk '{ sum+=$1} END {print sum}' > tmp.n_pass.txt
n_pass=`cat tmp.n_pass.txt`
masked_genome_size=`cat tmp.n_pass.txt`
rm -f tmp.n_pass.txt

# Calculate distance matrices.
gzvcf=../res/relate/input/gadmor.vcf.gz
bcftools view -O z -o tmp.vcf.gz ${gzvcf} LG03,LG04,LG05,LG06,LG08,LG09,LG10,LG11,LG13,LG14,LG15,LG16,LG17,LG18,LG19,LG20,LG21,LG22,LG23
gzvcf_id=`basename ${gzvcf%.vcf.gz}`
../bin/evo stats --diff-matrix tmp.vcf.gz
mv tmp.vcf.diff_me_matrix.txt ../res/evo/${gzvcf_id}.vcf.diff_me_matrix.txt
ruby reduce_distance_matrix.rb ../res/evo/${gzvcf_id}.vcf.diff_me_matrix.txt ../data/tables/groups.txt ../res/evo/${gzvcf_id}.per_group.txt
mean_dist=`cat ../res/evo/${gzvcf_id}.per_group.txt | tail -n 1 | cut -f 2`
rate=`echo "${mean_dist} / (${masked_genome_size} * 2 * ${divergence_time})" | bc -l`

# Clean up.
rm -f ${gzvcf_id}.vcf.*.txt

# Output.
echo "Divergence time: ${divergence_time}" > ../res/relate/input/rate.txt
echo "Genetic distance: ${mean_dist}" >> ../res/relate/input/rate.txt
echo "Genome size (exkl. masked sites): ${masked_genome_size}" >> ../res/relate/input/rate.txt
echo "Mutation rate per bp and per myr: ${rate}" >> ../res/relate/input/rate.txt
