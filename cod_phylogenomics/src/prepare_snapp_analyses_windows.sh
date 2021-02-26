# m_matschiner Thu Dec 5 13:51:19 CET 2019

# Copy the snapp xmls and slurm scripts to 2 replicate directories.
for xml in ../res/snapp/xmls/LG??_*_*.xml
do
    snapp_id=`basename ${xml%.xml}`
    for n in `seq -w 1 2`
    do
        snapp_dir=../res/snapp/windows/${snapp_id}/replicates/r0${n}
        mkdir -p ${snapp_dir}
        if [ ! -f ${snapp_dir}/${snapp_id}.xml ]
        then
            cp ${xml} ${snapp_dir}
            cat run_snapp.slurm | sed "s/QQQQQQ/${snapp_id}/g" > ${snapp_dir}/run_snapp.slurm
        fi
    done
done
