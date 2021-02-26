# m_matschiner Thu Dec 5 12:23:27 CET 2019

# Load modules.
module --quiet purge
module load BCFtools/1.9-foss-2018b
module load Ruby/2.6.1-GCCcore-7.3.0

# Make the snapp xml directory.
mkdir -p ../res/snapp/xmls

# Get the snapp_prep.rb script.
if [ ! -f snapp_prep.rb ]
then
    wget https://raw.githubusercontent.com/mmatschiner/snapp_prep/master/snapp_prep.rb
fi

# Make xml files for colinear and inverted regions.
for gzvcf in ../res/gatk/gadMor2/original_gadoga/colinear.vcf.gz ../res/gatk/gadMor2/original_gadoga/inversion_lg??.vcf.gz
do
    # Set the name of the xml file.
    xml_base=`basename ${gzvcf%.vcf.gz}.xml`
    xml=../res/snapp/xmls/${xml_base}

    # Only make an xml file if it doesn't already exist.
    if [ ! -f ${xml} ]
    then
	    # Uncompress the gzvcf.
        if [ ${xml_base} == "inversion_lg12.xml" ]
        then
            bcftools view -s ^Gadmor_lfc1,Gadmor_lfc2 ${gzvcf} > tmp.vcf
        else
            bcftools view ${gzvcf} > tmp.vcf
        fi
        
	    # Get the list of samples in the vcf.
	    samples=`head -n 10000 tmp.vcf | grep "#" | tail -n 1 | cut -f 10-`
        
        # Write a temporary samples file.
        echo -e "species\tspecimen" > tmp.spc.txt
        last_sample=""
        for sample in ${samples}
        do
            echo -e "${sample:0:10}_spc\t${sample}" >> tmp.spc.txt
            last_sample=${sample}
        done

        # Write a temporary constraints file.
        echo -ne "lognormal(0,3.83,0.093)\tcrown\t" > tmp.con.txt # The constraint is from the largest tree topology subset of the aim analysis with a prior distribution on the root age.
        samples_unique=`for sample in ${samples}; do echo ${sample:0:10}; done | sort | uniq`
        for sample in ${samples_unique}
        do
            echo -n "${sample:0:10}_spc" >> tmp.con.txt
            if [ ! ${sample} == ${last_sample} ]
            then
                echo -n "," >> tmp.con.txt
            fi
        done

	    # Generate the snapp xml files.
	    ruby snapp_prep.rb -v tmp.vcf -t tmp.spc.txt -c tmp.con.txt -l 100000 -m 1000 -x ${xml} -o ${xml_base%.xml}

	    # Clean up.
	    rm -f tmp.spc.txt
	    rm -f tmp.con.txt
	    rm -f tmp.vcf
    fi
done
