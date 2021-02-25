# m_matschiner Wed Jul 18 09:41:57 CEST 2018

# Load the ruby module.
module load ruby/2.1.5

# Make the output directory.
mkdir -p ../res/masks

# Make a mask for each finished alignment.
for alignment in ../res/alignments/merged/lg??.threaded.refined.finished.fasta
do
    chromosome_id=`basename ${alignment%.threaded.refined.finished.fasta}`
    mask=../res/masks/${chromosome_id}.bed
    ruby make_mask.rb ${alignment} ${chromosome_id} ${mask}
done

# Concatenate all masks into a single one.
cat ../res/masks/lg??.bed > ../res/masks/unreliable.bed
