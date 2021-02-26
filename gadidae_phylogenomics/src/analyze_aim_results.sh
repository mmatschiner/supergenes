# michaelm Mon Feb 8 17:43:01 CET 2021

# Load the modules.
module load R/3.6.0-foss-2019a
module load Python/3.7.2-GCCcore-8.2.0
module load Beast/2.5.2-GCC-8.2.0-2.31.1

# Set the r script.
logcombiner_script=`readlink -f logcombiner.py`

# Make the result directory.
mkdir -p ../res/aim/age_prior/combined
rm -f ../res/aim/age_prior/combined/log_files.txt
rm -f ../res/aim/age_prior/combined/tree_files.txt

# Make a temporary modified version of nicola's r script.
cat analyseAIMrun.R | sed "s/.\/combined\/species.trees/tmp.species.trees/g" | grep -v "this.dir" > tmp.script.r
cat tmp.script.r | sed "s/\/Applications\/BEAST//g" > tmp.script2.r
mv -f tmp.script2.r tmp.script.r
cat tmp.script.r | sed "s/\\\\\\ 2.5.0\/bin\///g" > tmp.script2.r
mv -f tmp.script2.r tmp.script.r
cat tmp.script.r | sed "s/\\\logcombiner/logcombiner/g" > tmp.script2.r
mv -f tmp.script2.r tmp.script.r
cat tmp.script.r | sed "s/rm(mig_rates)/suppressWarnings(rm(mig_rates))/g" > tmp.script2.r
mv -f tmp.script2.r tmp.script.r
cat tmp.script.r | sed "s/rm(arrows.data)/suppressWarnings(rm(arrows.data))/g" > tmp.script2.r
mv -f tmp.script2.r tmp.script.r
cat tmp.script.r | sed "s/paste(\"_\", i, \".pdf\", sep=\"\")/paste(\"_\", i, \"_\", bpp__proc_support[i], \".pdf\", sep=\"\")/g" > tmp.script2.r
mv -f tmp.script2.r tmp.script.r
cat tmp.script.r | sed "s/median(log.data/mean(log.data/g" > tmp.script2.r # This changes from median node heights to mean node heights.
mv -f tmp.script2.r tmp.script.r

# Analyze aim results.
converged_replicates=`echo r01 r02 r03 r04 r06 r07 r08 r09 r10`
for rep in ${converged_replicates}
do
    dir=../res/aim/age_prior/replicates/${rep}
    cp tmp.script.r ${dir}
    cd ${dir}
    python3 ${logcombiner_script} -b 10 aim_age_prior_species.trees tmp.species.trees
    rm -f mig_rates
    rm -f arrows.data
    Rscript tmp.script.r
    rm -f tmp.species.trees
    rm -f tmp.script.r
    rm -f Rplots.pdf
    rm -f mig_rates
    rm -f arrows.data
    for i in tmp.species*.{log,pdf}
    do
        new_name=`echo ${i} | sed 's/tmp.species/aim_species/g'`
        mv ${i} ${new_name}
    done
    cd -
    ls ${dir}/aim_age_prior.log >> ../res/aim/age_prior/combined/log_files.txt
    ls ${dir}/aim_age_prior_species.trees >> ../res/aim/age_prior/combined/tree_files.txt
done

# Combine results of replicate analyses with estimated root ages.
python3 ${logcombiner_script} -n 2000 -b 10 ../res/aim/age_prior/combined/log_files.txt ../res/aim/age_prior/combined/aim.log
python3 ${logcombiner_script} -n 2000 -b 10 ../res/aim/age_prior/combined/tree_files.txt ../res/aim/age_prior/combined/aim_species.trees
cp tmp.script.r ../res/aim/age_prior/combined
cd ../res/aim/age_prior/combined
python3 ${logcombiner_script} -b 0 aim_species.trees tmp.species.trees
rm -f mig_rates
rm -f arrows.data
Rscript tmp.script.r
rm -f tmp.species.trees
rm -f tmp.script.r
rm Rplots.pdf
rm -f mig_rates
rm -f arrows.data
for i in tmp.species*.{log,pdf}
do
    new_name=`echo ${i} | sed 's/tmp.species/aim_species/g'`
    mv ${i} ${new_name}
done
cd -

# Clean up.
rm -f tmp.script.r
