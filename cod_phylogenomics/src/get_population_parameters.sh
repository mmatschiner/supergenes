# m_matschiner Sun Jun 14 21:40:39 CEST 2020

# Make the result directory.
mkdir -p ../res/tables

# Load modules.
module load Ruby/2.7.1-GCCcore-8.3.0
module load BCFtools/1.10.2-GCC-8.3.0

# Set the vcf.
gzvcf=../res/gatk/gadMor2/original_gadoga/concatenated.masked.vcf.gz

# Get parameters for each region.
for region in colinear inversion_lg01 inversion_lg02 inversion_lg07 inversion_lg12
do
    # Extract the region.
    if [ ${region} == "colinear" ]
    then
        bcftools view -s ^Arcgla,Borsai --min-ac 1:minor -m2 -M2 -e 'F_MISSING > 0.2' -O z -o tmp.vcf.gz ${gzvcf} LG03,LG04,LG05,LG06,LG08,LG09,LG10,LG11,LG13,LG14,LG15,LG16,LG17,LG18,LG19,LG20,LG21,LG22,LG23
    elif [ ${region} == "inversion_lg01" ]
    then
        bcftools view -s ^Arcgla,Borsai --min-ac 1:minor -m2 -M2 -e 'F_MISSING > 0.2' -O z -o tmp.vcf.gz ${gzvcf} LG01:9114741-26192489
    elif [ ${region} == "inversion_lg02" ]
    then
        bcftools view -s ^Arcgla,Borsai --min-ac 1:minor -m2 -M2 -e 'F_MISSING > 0.2' -O z -o tmp.vcf.gz ${gzvcf} LG02:18489307-24050282
    elif [ ${region} == "inversion_lg07" ]
    then
        bcftools view -s ^Arcgla,Borsai --min-ac 1:minor -m2 -M2 -e 'F_MISSING > 0.2' -O z -o tmp.vcf.gz ${gzvcf} LG07:13606502-23016726
    elif [ ${region} == "inversion_lg12" ]
    then
        bcftools view -s ^Arcgla,Borsai,Gadmor_lfc1,Gadmor_lfc2 --min-ac 1:minor -m2 -M2 -e 'F_MISSING > 0.2' -O z -o tmp.vcf.gz ${gzvcf} LG12:589105-13631347
    fi

    # Make a table of all samples.
    samples=`bcftools view -h tmp.vcf.gz | tail -n 1 | cut -f 10-`
    rm -f tmp.samples.txt
    touch tmp.samples.txt
    for sample in ${samples}
    do
         echo -e "${sample:0:10}_spc\t${sample}" >> tmp.samples.txt
    done

    # Get the genome size.
    if [ ${region} == "colinear" ]
    then
        genome_size=`cat ../data/assemblies/gadMor2.dict | grep "SN:LG" | grep -v LG01 | grep -v LG02 | grep -v LG07 | grep -v LG12 | cut -f 3 | cut -d ":" -f 2 | awk '{ sum += $1 } END { print sum }'`
    elif [ ${region} == "inversion_lg01" ]
    then
        genome_size=`echo "26192489-9114741" | bc`
    elif [ ${region} == "inversion_lg02" ]
    then
        genome_size=`echo "24050282-18489307" | bc`
    elif [ ${region} == "inversion_lg07" ]
    then
        genome_size=`echo "23016726-13606502" | bc`
    elif [ ${region} == "inversion_lg12" ]
    then
        genome_size=`echo "13631347-589105"| bc`
    else
        echo "ERROR: Unexpected region!"
        exit 1
    fi

    # Get the unique population ids.
    pops=`cat tmp.samples.txt | cut -f 1 | sort | uniq`
    for pop in ${pops}
    do
        # Get the samples ids for this pop.
        cat tmp.samples.txt | grep ${pop} | cut -f 2 > tmp.samples_this_pop.txt

        # Make a vcf with only samples of this species.
        bcftools view -S tmp.samples_this_pop.txt -o tmp.pop.vcf tmp.vcf.gz

        # Calculate population parameters.
        ruby get_population_parameters_from_vcf.rb tmp.pop.vcf ${genome_size} ../res/tables/pop_stats_${region}_${pop}.txt

        # Clean up.
        rm -f tmp.pop.vcf
        rm -f tmp.samples_this_pop.txt
    done

    # Clean up.
    rm -f tmp.vcf.gz
    rm -f tmp.samples.txt
done
