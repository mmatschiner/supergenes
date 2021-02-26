# m_matschiner Tue Oct 9 22:56:28 CEST 2018

# Load the ruby and python modules.
module load ruby/2.1.5
module load python3/3.5.0

# Get stats for each alignment.
for align in ../res/windows/5000bp/*.phy
do
    # Set the info file name.
    info=${align%.phy}.info.txt

    # Feedback.
    align_base=`basename ${align}`
    echo -n "Analyzing alignment file ${align_base}..."

    # Convert the alignment to a temporary nexus file.
    python3 convert.py -f nexus ${align} tmp.nex

    # Check if the proportion of missing data has already been calculated,
    # and do so if it hasn't.
    p_missing_calculated=`cat ${info} | grep "p_missing:" | wc -l`
    if [[ ${p_missing_calculated} == 0 ]]
    then
        p_missing=`ruby get_proportion_of_missing_data.rb tmp.nex`
        echo "p_missing:${p_missing}" >> ${info}
    fi

    # Check if the number of variable sites has already been calculated,
    # and do so if it hasn't.
    n_variable_sites_calculated=`cat ${info} | grep "n_variable_sites:" | wc -l`
    if [[ ${n_variable_sites_calculated} == 0 ]]
    then
        n_variable_sites=`ruby get_number_of_variable_sites.rb tmp.nex`
        echo "n_variable_sites:${n_variable_sites}" >> ${info}
    fi

    # Check if the number of parsimony-informative sites has already been
    # calculated and do so if it hasn't.
    n_pi_sites_calculated=`cat ${info} | grep "n_pi_sites:" | wc -l`
    if [[ ${n_pi_sites_calculated} == 0 ]]
    then
        n_pi_sites=`ruby get_number_of_pi_sites.rb tmp.nex`
        echo "n_pi_sites:${n_pi_sites}" >> ${info}
    fi

    # Check if the number of hemiplasies has already been calculated, and
    # do so if it hasn't.
    n_hemiplasies_calculated=`cat ${info} | grep "n_hemiplasies:" | wc -l`
    if [[ ${n_hemiplasies_calculated} == 0 ]]
    then
        n_hemiplasies=`bash get_number_of_hemiplasies.sh tmp.nex`
        echo "n_hemiplasies:${n_hemiplasies}" >> ${info}
    fi
    echo " done."

    # Clean up.
    rm -f tmp.nex

done
