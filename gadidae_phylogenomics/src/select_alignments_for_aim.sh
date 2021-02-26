# m_matschiner Thu Oct 11 12:30:24 CEST 2018

# Set the input directory.
in_dir=../res/windows/5000bp

# Make the output directory.
out_dir=../res/aim/alignments
mkdir -p ${out_dir}
rm -f ${out_dir}/*

# Set the stats table.
stats_table=../res/tables/window_stats.txt

# Set the inversion table.
inversion_table=../data/tables/inversion_limits.txt

# Set thresholds for the selection of alignments.
min_n_variable_sites=100
max_n_hemiplasies=10

# Set the path to bmge.
bmge=../bin/BMGE.jar

# Set the log file.
log=../log/misc/bmge.log

# Select alignments.
while read line
do
    lg=`echo ${line} | tr -s " " | cut -d " "  -f 1`
    if [ ! ${lg} == "lg" ]
    then
        use_window=false
        from=`echo ${line} | tr -s " " | cut -d " "  -f 2`
        to=`echo ${line} | tr -s " " | cut -d " "  -f 3`
        n_variable_sites=`echo ${line} | tr -s " " | cut -d " "  -f 6`
        n_hemiplasies=`echo ${line} | tr -s " " | cut -d " "  -f 8`
        if (( $(echo "${n_variable_sites} < ${min_n_variable_sites}" | bc -l) ))
        then
            echo "Excluding alignment (${lg}_${from}_${to}) with ${n_variable_sites} variable sites."
            continue
        elif (( $(echo "${n_hemiplasies} > ${max_n_hemiplasies}" | bc -l) ))
        then
            echo "Excluding alignment (${lg}_${from}_${to}) with ${n_hemiplasies} hemiplasies."
            continue
        else
            lg_in_inversion_table=`cat ${inversion_table} | grep -v "#" | grep ${lg} | wc -l`
            if [[ ${lg_in_inversion_table} == 0 ]]
            then
                use_window=true
            elif [[ ${lg_in_inversion_table} == 1 ]]
            then
                inversion_begin=`cat ${inversion_table} | grep -v "#" | grep ${lg} | cut -f 2`
                inversion_end=`cat ${inversion_table} | grep -v "#" | grep ${lg} | cut -f 3`
                if [[ ${to} < ${inversion_begin} || ${from} > ${inversion_end} ]]
                then
                    use_window=true
                else
                    echo "Excluding alignment (${lg}_${from}_${to}) in inversion region."
                    continue
                fi
            else
                echo "ERROR: Unexpected number of occurrences of lg ${lg} in file ${inversion_table}!"
                exit 1
            fi
        fi
        if [ ${use_window} == true ]
        then
            # Run bmge to remove sites with large gap proportion.
            echo -n "Running BMGE for alignment (${lg}_${from}_${to})... "
            java -jar ${bmge} -i ${in_dir}/${lg}_${from}_${to}.phy -t DNA -o ${out_dir}/${lg}_${from}_${to}.phy &> ${log}
            echo " done."
        fi
    fi
done < ${stats_table}

# Report how many alignments were copied.
count=`ls ${out_dir}/*.phy | wc -l`
echo "Copied ${count} alignments from ${in_dir} to ${out_dir}."
