# m_matschiner Tue Sep 1 00:46:43 CEST 2020

# Load modules.
module load Ruby/2.7.1-GCCcore-8.3.0

# Set the dict.
dict=../data/assemblies/gadMor2.dict

# Set the repeat content table.
mkdir -p ../res/tables
table=../res/tables/repeat_content.txt
rm -f ${table}
touch ${table}

# Use a ruby script to determine the repeat content in regions.
window_size=1000000
for lg in LG0{1..9} LG{10..23}
do
    window_start=1
    window_end=$(( ${window_start} + ${window_size} - 1 ))
    lg_length=`cat ${dict} | grep "SN:${lg}" | cut -f 3 | cut -d ":" -f 2`
    while (( ${window_end} < ${lg_length} ))
    do
        ruby check_repetitive_content.rb ../data/masks/Atlantic_cod_repeats.tab ${lg} ${window_start} ${window_end} >> ${table}
        window_start=$(( ${window_end} + 1 ))
        window_end=$(( ${window_start} + ${window_size} - 1 ))
    done
done
