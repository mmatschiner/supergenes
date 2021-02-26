# m_matschiner Sun Feb 21 22:10:31 CET 2021

# Load modules.
module load Ruby/2.7.1-GCCcore-8.3.0

# Plot the dbbaa output of all dsuite analyses.
for dbbaa_res in ../res/dsuite/*/*_gadoga/*/*__BBAA.txt 
do
    region_id=`basename ${dbbaa_res%__BBAA.txt}`
    order=../data/tables/dsuite_order_${region_id}.txt
    combine_res=${dbbaa_res%__BBAA.txt}__combine.txt
    ruby plot_d.rb ${dbbaa_res} ${order} 0.35 ${dbbaa_res%.txt}.svg ${combine_res}
done

# Plot the dfix output of all dsuite analyses.
for dfix_res in ../res/dsuite/*/*_gadoga/*/*__tree.txt
do
    region_id=`basename ${dfix_res%__tree.txt}`
    order=../data/tables/dsuite_order_${region_id}.txt
    combine_res=${dfix_res%__tree.txt}__combine.txt
    ruby plot_d.rb ${dfix_res} ${order} 0.35 ${dfix_res%.txt}.svg ${combine_res}
done
