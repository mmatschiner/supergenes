# michaelm Wed Mar 18 23:44:48 CET 2020

# Load modules.
module load Beast/2.5.2-GCC-8.2.0-2.31.1
module load Ruby/2.6.3-GCCcore-8.2.0
module load Python/3.7.2-GCCcore-8.2.0

# Get the script to add theta to log files.
if [ ! -f add_theta_to_log.rb ]
then
    wget https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/add_theta_to_log.rb
fi

# Make the log directory.
mkdir -p ../log/snapp
out=../log/snapp/combine.out

# Combine results for each snapp analysis.
for dir in ../res/snapp/windows/LG??_*
do
    # Combine replicate log files for the snapp analyses.  
    window_id=`basename ${dir}`
    echo -n "Combining analyses for window ${window_id}..."
    rm -rf ../res/snapp/windows/${window_id}/combined
    mkdir -p ../res/snapp/windows/${window_id}/combined
    log1=../res/snapp/windows/${window_id}/replicates/r01/${window_id}.log
    trees1=../res/snapp/windows/${window_id}/replicates/r01/${window_id}.trees
    log2=../res/snapp/windows/${window_id}/replicates/r02/${window_id}.log
    trees2=../res/snapp/windows/${window_id}/replicates/r02/${window_id}.trees
    if [ ! -w ${log1%.log}_w_theta.log ]
    then
	    ruby add_theta_to_log.rb -l ${log1} -t ${trees1} -g 5 -o ${log1%.log}_w_theta.log &>> ${out}
	    ruby add_theta_to_log.rb -l ${log2} -t ${trees1} -g 5 -o ${log2%.log}_w_theta.log &>> ${out}
    fi
    if [ ! -f ${log1%.log}_w_theta_2000.log ]
    then
	    python logcombiner.py -b 10 -n 2000 ${log1%.log}_w_theta.log ${log1%.log}_w_theta_2000.log
	    python logcombiner.py -b 10 -n 2000 ${log2%.log}_w_theta.log ${log2%.log}_w_theta_2000.log
    fi
    if [ ! -f ${trees1%.trees}_2000.trees ]
    then
	    python logcombiner.py -b 10 -n 2000 ${trees1} ${trees1%.trees}_2000.trees
	    python logcombiner.py -b 10 -n 2000 ${trees2} ${trees2%.trees}_2000.trees
    fi
    if [ ! -f ../res/snapp/windows/${window_id}/combined/${window_id}.log ]
    then
	    logcombiner -log ${log1%.log}_w_theta_2000.log -log ${log2%.log}_w_theta_2000.log -b 0 -o ../res/snapp/windows/${window_id}/combined/${window_id}.log &>> ${out}
	    logcombiner -log ${trees1%.trees}_2000.trees -log ${trees2%.trees}_2000.trees -b 0 -o ../res/snapp/windows/${window_id}/combined/${window_id}.trees &>> ${out}
    fi

    # Make maximum-clade-credibility consensenssus trees.
    if [ ! -f ../res/snapp/windows/${window_id}/combined/${window_id}.tre ]
    then
	    treeannotator -b 0 -heights mean ../res/snapp/windows/${window_id}/combined/${window_id}.trees ../res/snapp/windows/${window_id}/combined/${window_id}.tre &>> ${out}
    fi
    echo " done."
done
