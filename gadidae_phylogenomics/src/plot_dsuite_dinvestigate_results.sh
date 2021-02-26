# m_matschiner Sat Nov 2 15:17:07 CET 2019

# Load modules.
module load R/3.4.4

# Plot the results in each dsuite dinvestigate result table as manhattan plot.
for table in ../res/tables/*__50_25.txt
do
    table_id=`basename ${table%.txt}`
    plot=../res/plots/${table_id}.png
    Rscript plot_dsuite_dinvestigate_results.r ${table} ${plot}
done