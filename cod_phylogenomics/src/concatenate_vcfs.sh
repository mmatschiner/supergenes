# m_matschiner Sun Sep 3 19:49:36 CEST 2017

# Repeat gatk1 for each of the alternative references.
for ref in ../data/assemblies/gadMor2.fasta ../data/assemblies/melAeg_in_gadMor2_coords.fasta
do
    # Set the reference id.
    ref_id=`basename ${ref%.fasta}`

    # Repeat first with the original gadoga specimen and then with the alternative one.
    for gadoga_setting in original_gadoga alternative_gadoga
    do

        # Use GATKs CatVariants to concatenate all vcfs of all lgs.
        gzvcf_dir=../res/gatk/${ref_id}/${gadoga_setting}
        concatenated_vcfgz=${gzvcf_dir}/concatenated.masked.vcf.gz
        out="../log/misc/concatenate_vcfs.${ref_id}.${gadoga_setting}.out"
        rm -f ${out}
        sbatch -o ${out} concatenate_vcfs.slurm ${ref} ${gzvcf_dir} ${concatenated_vcfgz}
    done
done
