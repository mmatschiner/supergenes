# m_matschiner Sun Aug 30 14:40:12 CEST 2020

# Make the input directory for relate.
mkdir -p ../res/relate/input

# Load modules.
module load BCFtools/1.10.2-iccifort-2019.5.281
module load BEDTools/2.29.2-GCC-8.3.0
module load Ruby/2.7.1-GCCcore-8.3.0

# Set the compressed input vcf.
in_gzvcf=../res/beagle/concatenated.filtered.anc_hets_removed.vcf.gz

# Set the reference.
ref=../data/assemblies/gadMor2.fasta

# Set the ancestral reference.
anc_ref=../res/relate/input/ancestor.fasta

# Set the output vcf without ancestor.
out_gzvcf=../res/relate/input/gadmor.vcf.gz

# Set the ancestor id.
anc_id="Gadcha"

# Make a temporary vcf with only the ancestor, and exluding unplaced contigs
bcftools view -o tmp1.vcf -s "${anc_id}" ${in_gzvcf} LG01,LG02,LG03,LG04,LG05,LG06,LG07,LG08,LG09,LG10,LG11,LG12,LG13,LG14,LG15,LG16,LG17,LG18,LG19,LG20,LG21,LG22,LG23

# Add contig lengths to the vcf.
java -jar $EBROOTPICARD/picard.jar UpdateVcfSequenceDictionary \
    INPUT=tmp1.vcf \
    OUTPUT=tmp2.vcf \
    SEQUENCE_DICTIONARY=${ref%.fasta}.dict
mv -f tmp2.vcf tmp1.vcf

# Make an ancestral reference.
gatk IndexFeatureFile -F tmp1.vcf
gatk FastaAlternateReferenceMaker \
    -R ${ref} \
    -O tmp.fasta \
    -V tmp1.vcf

# Correct the ids in the ancestral reference.
cat tmp.fasta | sed -e 's/:.*//g' | sed -e 's/>.*LG/>LG/g' | awk '/MT_genome/ {exit} {print}' > ${anc_ref}

# Clean up.
rm -f ${anc_ref}.fai
rm -f ${anc_ref%.fasta}.dict
rm -f tmp.fasta

# Make a version of the reference that exludes unplaced contigs.
rm -f tmp.ref.fasta
touch tmp.ref.fasta
for lg in LG0{1..9} LG{10..23}
do
    ./fastagrep -t -p ${lg} ${ref} >> tmp.ref.fasta
done

# Make the ancestral sequence.
ruby make_fasta_from_ref_and_vars.rb tmp1.vcf tmp.ref.fasta Gadcha ${anc_ref}

# Clean up.
rm tmp1.vcf

# Write a vcf file without the ancestor.
bcftools view -s "^${anc_id}" ${in_gzvcf} LG01,LG02,LG03,LG04,LG05,LG06,LG07,LG08,LG09,LG10,LG11,LG12,LG13,LG14,LG15,LG16,LG\
17,LG18,LG19,LG20,LG21,LG22,LG23 | bcftools view -O z -o ${out_gzvcf} -e 'AC==0 || AC==AN'
tabix ${out_gzvcf}

# Make a mask in fasta format combining the indel masks, the general mask, and the mask for sites heterozygous in the ancestor.
for bed in ../res/gatk/LG??.filtered.indels.bed ../res/beagle/anc_hets.bed ../data/masks/Atlantic_cod_repeats.tab
do
    tail -n +2 ${bed} > tmp.bed
    bedtools maskfasta -fi tmp.ref.fasta -bed tmp.bed -fo tmp2.ref.fasta
    ruby remove_line_breaks_from_fasta.rb tmp2.ref.fasta tmp.ref.fasta
    rm -f tmp2.ref.fasta
done
cat tmp.ref.fasta | sed 's/>LG/>LX/g' | sed 's/[ACGT]/P/g' | sed 's/>LX/>LG/g' > ../res/relate/input/mask.fasta

# Clean up.
rm -f tmp.ref.fasta
