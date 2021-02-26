# m_matschiner Tue Oct 16 18:00:59 CEST 2018

# Make the result directory.
mkdir -p ../res/aim/age_prior/replicates

# Prepare the analysis directories.
xml=../res/aim/xml/aim_age_prior.xml
xml_base=`basename ${xml}`
for n in `seq -w 1 10`
do
    replicate_dir=../res/aim/age_prior/replicates/r${n}
    mkdir -p ${replicate_dir}
    cp ${xml} ${replicate_dir}/aim.xml
    cp run_aim.slurm ${replicate_dir}
done
