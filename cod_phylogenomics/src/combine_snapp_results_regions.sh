# m_matschiner Thu Dec 5 14:42:15 CET 2019

# Load modules.
module load Beast/2.5.2-GCC-8.2.0-2.31.1
module load Ruby/2.6.3-GCCcore-8.2.0
module load Python/3.7.2-GCCcore-8.2.0

# Get the script to add theta to log files.
if [ ! -f add_theta_to_log.rb ]
then
    wget https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/add_theta_to_log.rb
fi

# Combine results for each snapp analysis.
for xml in ../res/snapp/xmls/colinear.xml ../res/snapp/xmls/inversion_lg??.xml
do
    # Combine replicate log files for the snapp analyses.
    xml_base=`basename ${xml%.xml}`
    mkdir -p ../res/snapp/regions/${xml_base}/combined
    ls ../res/snapp/regions/${xml_base}/replicates/r??/*.log > ../res/snapp/regions/${xml_base}/combined/logs.txt
    ls ../res/snapp/regions/${xml_base}/replicates/r??/*.trees > ../res/snapp/regions/${xml_base}/combined/trees.txt
    python3 logcombiner.py -n 1000 -b 20 ../res/snapp/regions/${xml_base}/combined/logs.txt ../res/snapp/regions/${xml_base}/combined/${xml_base}.log
    python3 logcombiner.py -n 1000 -b 20 ../res/snapp/regions/${xml_base}/combined/trees.txt ../res/snapp/regions/${xml_base}/combined/${xml_base}.trees

    # Add population sizes to log files.
    ruby add_theta_to_log.rb -l ../res/snapp/regions/${xml_base}/combined/${xml_base}.log -t ../res/snapp/regions/${xml_base}/combined/${xml_base}.trees -g 5 -o ../res/snapp/regions/${xml_base}/combined/${xml_base}_w_theta.log

    # Make maximum-clade-credibility consensenssus trees.
    treeannotator -b 0 -heights mean ../res/snapp/regions/${xml_base}/combined/${xml_base}.trees ../res/snapp/regions/${xml_base}/combined/${xml_base}.tre
done
