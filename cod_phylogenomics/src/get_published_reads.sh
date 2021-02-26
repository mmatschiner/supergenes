# m_matschiner Wed Oct 23 14:49:45 CEST 2019

# Make the data directory if it doesn't exist yet.
mkdir -p ../data/reads

# Try to load the aspera module.
module load Aspera-CLI/3.9.0.1326.6985b21

# Check that the aspera connect command-line utility is installed.
if ! [ -x "$(command -v ascp)" ]
then
    echo "ERROR: The aspera connect command-line utility must be installed to download from the European Nucleotide Archive."
    exit
fi

# Try to set the aspera key.
asp_key=/cluster/software/Aspera-CLI/3.9.0.1326.6985b21/etc/asperaweb_id_dsa.openssh
if [ ! -f ${asp_key} ]
then
    echo "ERROR: Aspera key not found!"
    echo "  Replace key path in this file."
    exit 1
fi

# Set the read table file.
read_table_ah19=../data/tables/PRJNA301002.txt
read_table_m16=../data/tables/PRJEB12469.txt
read_table_b19=../data/tables/PRJEB29231.txt
read_table_t17=../data/tables/ERX1622640_ERX1622647.txt
read_table_k16=../data/tables/ERX1350534.txt

# Select fastq files from arnason and halldorsdottir (2019).
tail -n +2 ${read_table_ah19} | cut -f 5,7,10 | grep -e SRR2906193 -e SRR2906345 -e SRR2906368 -e SRR2906362 -e SRR2906197 > tmp.txt
# Select fastq files from malmstrøm et al. (2016).
tail -n +2 ${read_table_m16} | cut -f 5,7,10 | grep -e ERR1473882 -e ERR1473883 -e ERR1473884 -e ERR1473885 -e ERR1473886 >> tmp.txt
# Select fastq files from Torresen et al. (2017).
tail -n +2 ${read_table_t17} | grep ERR1551885 | cut -f 4,6,7 >> tmp.txt
# Select fastq files from barth et al. (2019).
tail -n +2 ${read_table_b19} | grep -e AVE1409010 -e BOR1205007 -e KIE1103020 -e LOW1504007 -e TVE120514 | cut -f 4,6,7 >> tmp.txt
# Select fastq files from kirubakaran et al. (2016).
tail -n +2 ${read_table_k16} | grep ERR1278928 | cut -f 5,7,10 >> tmp.txt

# Download all selected fastq files.
while read line
do
    specimen_id=`echo ${line} | cut -d " " -f 1`
    species1=`echo ${line} | cut -d " " -f 2`
    species2=`echo ${line} | cut -d " " -f 3`
    link1=`echo ${line} | cut -d " " -f 4 | cut -d ";" -f 1 | sed 's/ftp.sra.ebi.ac.uk/era-fasp@fasp.sra.ebi.ac.uk:/g'`
    link1_file=`basename ${link1}`
    link2=`echo ${line} | cut -d " " -f 4 | cut -d ";" -f 2 | sed 's/ftp.sra.ebi.ac.uk/era-fasp@fasp.sra.ebi.ac.uk:/g'`
    link2_file=`basename ${link2}`
    species_id="${species1:0:3}${species2:0:3}"
    if [ ! -f ../data/reads/${species_id}_${specimen_id}_R1.fastq.gz ]
    then
        ascp -QT -l 300m -P33001 -i ${asp_key} ${link1} .
        ascp -QT -l 300m -P33001 -i ${asp_key} ${link2} .
        if [[ -f ${link1_file} && -f ${link2_file} ]] # Check that both files have been downloaded.
        then
            # Get number of lines in both files.
            link1_file_length=`zcat ${link1_file} | wc -l`
            link2_file_length=`zcat ${link2_file} | wc -l`
            if (( ${link1_file_length} == ${link2_file_length} )) # Check that both files have the same number of lines.
            then
                mv ${link1_file} ../data/reads/${species_id}_${specimen_id}_R1.fastq.gz
                mv ${link2_file} ../data/reads/${species_id}_${specimen_id}_R2.fastq.gz
            else
                echo "WARNING: Incomplete download for ${species_id}_${specimen_id}!"
                echo "  Files ${link1_file} and ${link2_file} have ${link1_file_length} and ${link2_file_length} lines."
                echo "  Rerun this script!"
                rm ${link1_file}
                rm ${link2_file}
            fi
        else
            echo "WARNING: Files could not be downloaded. Rerun this script!"
        fi
    fi
done < tmp.txt

# Clean up.
rm -f tmp.txt
