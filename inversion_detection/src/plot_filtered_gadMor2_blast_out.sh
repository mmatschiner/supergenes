# m_matschiner Mon Jan 16 23:15:58 CET 2017

# Load modules.
module load ruby/2.1.5

# Plot the filtered BLAST output for the gadMor_Stat assembly.
for i in ../res/blast/gadMor2_vs_gadMor_Stat/lg*_filtered.out
do
    target_region_id=`basename ${i%_filtered.out}`
    lg=`echo ${target_region_id} | cut -d "_" -f 1 | sed 's/lg/LG/g'`
    ruby plot_filtered_gadMor2_blast_out.rb ${i} ../res/plink/${lg}.ld ../data/assemblies/gadMor2_${target_region_id}_features.txt ../res/breakdancer/summary/links.txt ${i%.out}.svg 1000
done

# Plot the filtered BLAST output for the melAeg assembly.
for i in ../res/blast/gadMor2_vs_melAeg/lg*_filtered.out
do
    target_region_id=`basename ${i%_filtered.out}`
    lg=`echo ${target_region_id} | cut -d "_" -f 1 | sed 's/lg/LG/g'`
    ruby plot_filtered_gadMor2_blast_out.rb ${i} ../res/plink/${lg}.ld ../data/assemblies/gadMor2_${target_region_id}_features.txt ../res/breakdancer/summary/links.txt ${i%.out}.svg 1000
done
