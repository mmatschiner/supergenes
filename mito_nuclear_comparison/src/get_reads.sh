# m_matschiner Wed Oct 23 14:49:45 CEST 2019

# Make the data directory if it doesn't exist yet.
mkdir -p ../data/reads

# Set the read table file.
read_table1=../data/tables/PRJNA301002.txt
read_table2=../data/tables/PRJEB12469.txt
read_table3=../data/tables/ERX1350534.txt

# Download the fastq files for the high-coverage specimens plus arctogadus from arnason and halldorsdottir (2019).
tail -n +2 ${read_table1} | cut -f 5,7,10 | grep -e SRR2906188 -e SRR2906185 -e SRR2906193 -e SRR2906345 -e SRR2906339 -e SRR2906208 -e SRR2906256 -e SRR2906368 -e SRR2906362 -e SRR2906197 -e SRR2906248 > tmp.txt
#  Download fastq files from malmstrøm et al. (2016).
tail -n +2 ${read_table2} | cut -f 5,7,10 | grep -e ERR1473878 -e ERR1473879 -e ERR1473880 -e ERR1473881 -e ERR1473882 -e ERR1473883 -e ERR1473884 -e ERR1473885 -e ERR1473886 -e ERR1473887 -e ERR1473888 -e ERR1473889 -e ERR1473890 -e ERR1473875 -e ERR1473876 -e ERR1473877 >> tmp.txt
# Download fastq files for gadoga from kirubakaran et al. (2016).
tail -n +2 ${read_table3} | cut -f 5,7,10 >> tmp.txt
while read line
do
    specimen_id=`echo ${line} | cut -d " " -f 1`
    species1=`echo ${line} | cut -d " " -f 2`
    species2=`echo ${line} | cut -d " " -f 3`
    link1=`echo ${line} | cut -d " " -f 4 | cut -d ";" -f 1`
    link1_file=`basename ${link1}`
    link2=`echo ${line} | cut -d " " -f 4 | cut -d ";" -f 2`
    link2_file=`basename ${link2}`
    species_id="${species1:0:3}${species2:0:3}"
    if [ ! -f ../data/reads/${species_id}_${specimen_id}_R1.fastq.gz ]
    then
        wget ${link1}
        mv ${link1_file} ../data/reads/${species_id}_${specimen_id}_R1.fastq.gz
    fi
    if [ ! -f ../data/reads/${species_id}_${specimen_id}_R2.fastq.gz ]
    then
        wget ${link2}
        mv ${link2_file} ../data/reads/${species_id}_${specimen_id}_R2.fastq.gz
    fi
done < tmp.txt

# Clean up.
rm -f tmp.txt