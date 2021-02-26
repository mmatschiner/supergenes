# m_matschiner Sun Dec 16 13:37:12 CET 2018

# Load modules.
module load ruby/2.1.5

# Set the inversion table.
inversion_table=../data/tables/inversion_limits.txt

# Plot results for likelihood comparisons with all datasets.
for hypothesis in h01 h02
do
    # Set the split tables. 
    full_table=../res/tables/constrained_likelihoods_${hypothesis}.txt
    no_inversion_table=../res/tables/constrained_likelihoods_${hypothesis}.no_inversion.txt
    lg01_inversion_table=../res/tables/constrained_likelihoods_${hypothesis}.lg01_inversion.txt
    lg02_inversion_table=../res/tables/constrained_likelihoods_${hypothesis}.lg02_inversion.txt
    lg07_inversion_table=../res/tables/constrained_likelihoods_${hypothesis}.lg07_inversion.txt

    # Add the header to each split table.
    head -n 1 ${full_table} > ${no_inversion_table}
    head -n 1 ${full_table} > ${lg01_inversion_table}
    head -n 1 ${full_table} > ${lg02_inversion_table}
    head -n 1 ${full_table} > ${lg07_inversion_table}

    # Make a temporary version of the full table without the header line.
    tail -n +2 ${full_table} > tmp.constrained_likelihoods.txt

    # Bin the lines of the constrained-likelihoods file.
    while read line
    do
        # Get the coordinates of this tree.
        lg=`echo ${line} | cut -d " " -f 1 | cut -d "_" -f 1 | tr '[:lower:]' '[:upper:]'`
        from=`echo ${line} | cut -d " " -f 1 | cut -d "_" -f 2`
        to=`echo ${line} | cut -d " " -f 1 | cut -d "_" -f 3`

        # Check if the lg is one of those with inversions.
        lg_in_inversion_table=`cat ${inversion_table} | grep -v "#" | grep ${lg} | wc -l`

        # If the lg is one of those with inversions, check if the tree is from within the inversion.
        if [[ ${lg_in_inversion_table} == 0 ]]
        then
            echo ${line} >> ${no_inversion_table}
        elif [[ ${lg_in_inversion_table} == 1 ]]
        then
            inversion_begin=`cat ${inversion_table} | grep -v "#" | grep ${lg} | cut -f 2`
            inversion_end=`cat ${inversion_table} | grep -v "#" | grep ${lg} | cut -f 3`
            if (( ${to} < ${inversion_begin} || ${from} > ${inversion_end} ))
            then
                echo ${line} >> ${no_inversion_table}
            else
                if [ ${lg} == "LG01" ]
                then
                    echo ${line} >> ${lg01_inversion_table}
                elif [ ${lg} == "LG02" ]
                then
                    echo ${line} >> ${lg02_inversion_table}
                elif [ ${lg} == "LG07" ]
                then
                    echo ${line} >> ${lg07_inversion_table}
                else
                            echo "ERROR: Unexpected linkage group ${lg}!"
                            exit          
                fi
            fi
        else
            echo "ERROR: Unexpected number of occurrences of lg ${lg} in file ${inversion_table}!"
            exit 1
        fi
    done < tmp.constrained_likelihoods.txt

    # Make separate plots for the binned tree sets.
    for table in ${full_table} ${no_inversion_table} ${lg01_inversion_table} ${lg07_inversion_table}
    do
        table_id=`basename ${table%.txt}`
        plot=../res/plots/${table_id}.svg
        ruby plot_constrained_likelihoods.rb ${table} ${plot}
    done
done
