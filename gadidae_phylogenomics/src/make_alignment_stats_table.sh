# m_matschiner Thu Oct 11 00:46:45 CEST 2018

# Set the name of the stats table.
stats_table=../res/tables/window_stats.txt

# Write the header row.
echo -e "lg\tfrom\tto\tn_sites\tp_missing\tn_variable_sites\tn_pi_sites\tn_hemiplasies" > ${stats_table}

for align in ../res/windows/5000bp/*.phy
do
    # Set the name of the info file.
    info=${align%.phy}.info.txt

    # Get the window positions.
    lg=`basename ${info} | cut -d "." -f 1 | cut -d "_" -f 1`
    from=`basename ${info} | cut -d "." -f 1 | cut -d "_" -f 2`
    to=`basename ${info} | cut -d "." -f 1 | cut -d "_" -f 3`
    
    # Get the window stats.
    n_sites=`cat ${info} | grep "n_sites:" | cut -d ":" -f 2`
    p_missing=`cat ${info} | grep "p_missing:" | cut -d ":" -f 2`
    n_variable_sites=`cat ${info} | grep "n_variable_sites:" | cut -d ":" -f 2`
    n_pi_sites=`cat ${info} | grep "n_pi_sites:" | cut -d ":" -f 2`
    n_hemiplasies=`cat ${info} | grep "n_hemiplasies:" | cut -d ":" -f 2`

    # Add alignment stats to the stats table.
    echo -e "${lg}\t${from}\t${to}\t${n_sites}\t${p_missing}\t${n_variable_sites}\t${n_pi_sites}\t${n_hemiplasies}" >> ${stats_table}

done