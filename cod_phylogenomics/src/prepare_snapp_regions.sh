# michaelm Sun Jun 7 00:08:43 CEST 2020

# Copy the snapp xmls and slurm scripts to three replicate directories.
for xml in ../res/snapp/xmls/colinear.xml ../res/snapp/xmls/inversion_lg??.xml
do
    snapp_id=`basename ${xml%.xml}`
    for n in `seq -w 1 3`
    do
	    snapp_dir=../res/snapp/regions/${snapp_id}/replicates/r0${n}
	    mkdir -p ${snapp_dir}
	    if [ ! -f ${snapp_dir}/${snapp_id}.xml ]
	    then
	        cp ${xml} ${snapp_dir}
	        cat run_snapp.slurm | sed "s/QQQQQQ/${snapp_id}/g" > ${snapp_dir}/run_snapp.slurm
	    fi
    done
done
