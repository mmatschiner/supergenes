# m_matschiner Thu Mar 5 14:40:16 CET 2020

# Set and make the output directory.
res_dir=../res/snapp/xmls
mkdir -p ${res_dir}

# Make the log directory.
mkdir -p ../log/snapp/make_xmls/

# Set the window size.
window_size=250000

# Set the minimum number of sites for a window to be used.
min_n_sites=500

# Set the maximum number of sites to be used by snapp per window.
max_n_sites=1000

# Set the minimum distance between sites.
min_site_dist=50

# Make sliding-window xml files for each linkage group.
for gzvcf in ../res/gatk/gadMor2/original_gadoga/LG12.filtered.vcf.gz
do
    lg_id=`basename ${gzvcf%.filtered.vcf.gz}`
    out=../log/snapp/make_xmls/${lg_id}.out
    log=../log/snapp/make_xmls/${lg_id}.log
    rm -f ${out}
    rm -f ${log}
    sbatch -o ${out} make_snapp_xmls_windows.slurm ${gzvcf} ${res_dir} ${window_size} ${min_n_sites} ${max_n_sites} ${min_site_dist} ${log}
done
