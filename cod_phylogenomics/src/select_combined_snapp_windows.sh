# michaelm Thu Mar 19 22:34:40 CET 2020

# Load modules.
module load R/4.0.0-foss-2020a
module load Python/3.8.2-GCCcore-9.3.0

# Make the results directories.
mkdir -p ../res/snapp/trees
mkdir -p ../res/snapp/logs

# Make a temporary list of species to keep in trees.
cat ../data/tables/color_code.txt | cut -f 1 > tmp.keep_list.txt

# Copy log and trees files if ess values are sufficient.
for combined_log in ../res/snapp/windows/LG??_*/combined/*.log
do
    window_id=`echo ${combined_log} | cut -d "/" -f 5`
    tre=${combined_log%.log}.tre
    tre_base=`basename ${tre}`
    min_ess=`Rscript get_min_ess.r ${combined_log}`
    min_ess=`echo ${min_ess}`
    mean_bpp=`python get_mean_node_support.py ${tre}`
    xml=../res/snapp/windows/${window_id}/replicates/r01/${window_id}.xml
    n_bp=`cat ${xml} | head | grep Arcgla | cut -d "\"" -f 8 | wc -c`
    if (( $(echo "${min_ess} > 100" | bc -l) ))
    then
        if (( $(echo "${mean_bpp} > 0.5" | bc -l) ))
        then
            if (( $(echo "${n_bp} > 300" | bc -l) ))
            then
                trees=${combined_log%.log}.trees
                trees_base=`basename ${trees}`
                cp -f ${combined_log} ../res/snapp/logs
                if [ ! -f ../res/snapp/trees/${tre_base} ]
                then
                    Rscript convert_nexus_trees_to_newick.r ${tre} tmp.tre
                    Rscript prune_and_ladderize_newick.r tmp.tre tmp.keep_list.txt ../res/snapp/trees/${tre_base} &> /dev/null
                    rm -f tmp.tre
                fi
                if [ ! -f ../res/snapp/trees/${trees_base} ]
                then
                    Rscript convert_nexus_trees_to_newick.r ${trees} tmp.trees
                    shuf tmp.trees -n 100 > tmp2.trees; mv -f tmp2.trees tmp.trees
                    Rscript prune_and_ladderize_newick.r tmp.trees tmp.keep_list.txt ../res/snapp/trees/${trees_base} &> /dev/null
                    rm -f tmp.trees
                fi
                echo "Copied log and trees files for window ${window_id} (min. ESS = ${min_ess})."
            else
                echo "Skipped log and trees files for window ${window_id} (\#bp = ${n_bp})."
            fi
        else
            echo "Skipped log and trees files for window ${window_id} (mean BPP = ${mean_bpp})."
        fi
    else
        echo "Skipped log and trees files for window ${window_id} (min. ESS = ${min_ess})."
    fi
done

# Clean up.
rm -f tmp.keep_list.txt
