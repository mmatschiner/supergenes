# m_matschiner Tue Dec 3 17:42:14 CET 2019

# Load modules.
module purge
module load BCFtools/1.9-foss-2018b

# Repeat gatk1 for each of the alternative references.
for ref in ../data/assemblies/gadMor2.fasta ../data/assemblies/melAeg_in_gadMor2_coords.fasta
do
    # Set the reference id.
    ref_id=`basename ${ref%.fasta}`

    # Repeat first with the original gadoga specimen and then with the alternative one.
    for gadoga_setting in original_gadoga alternative_gadoga
    do

        # Set the vcf file.
        gzvcf=../res/gatk/${ref_id}/${gadoga_setting}/concatenated.masked.vcf.gz

        # Remove all sites that are not complete and bi-allelic.
        echo -n "Extracting inversion on LG01..."
        bcftools view --min-ac 1:minor -m2 -M2 -e 'GT[*] = "mis"' -O z -o ../res/gatk/${ref_id}/${gadoga_setting}/inversion_lg01.vcf.gz ${gzvcf} LG01:9114741-26192386
        echo " done."
        echo -n "Extracting inversion on LG02..."
        bcftools view --min-ac 1:minor -m2 -M2 -e 'GT[*] = "mis"' -O z -o ../res/gatk/${ref_id}/${gadoga_setting}/inversion_lg02.vcf.gz ${gzvcf} LG02:18489307-24048607
        echo " done."
        echo -n "Extracting inversion on LG07..."
        bcftools view --min-ac 1:minor -m2 -M2 -e 'GT[*] = "mis"' -O z -o ../res/gatk/${ref_id}/${gadoga_setting}/inversion_lg07.vcf.gz ${gzvcf} LG07:13610591-23019113
        echo " done."
        echo -n "Extracting inversion on LG12..."
        bcftools view --min-ac 1:minor -m2 -M2 -e 'GT[*] = "mis"' -O z -o ../res/gatk/${ref_id}/${gadoga_setting}/inversion_lg12.vcf.gz ${gzvcf} LG12:638100-14327837
        echo " done."
        echo -n "Extracting colinear regions..."
        bcftools view --min-ac 1:minor -m2 -M2 -e 'GT[*] = "mis"' -O z -o ../res/gatk/${ref_id}/${gadoga_setting}/colinear.vcf.gz ${gzvcf} LG01:1-9114740,LG01:26192387-28298372,LG02:1-18489306,LG03,LG04,LG05,LG06,LG07:1-13610590,LG07:23019114-31229947,LG08,LG09,LG10,LG11,LG12:1-638099,LG12:14327838-27174696,LG13,LG14,LG15,LG16,LG17,LG18,LG19,LG20,LG21,LG22,LG23
        echo " done."

        # Write separate files with sites that are biallelic, allowing missing data.
        echo -n "Extracting inversion on LG01..."
        bcftools view --min-ac 1:minor -m2 -M2 -O z -o ../res/gatk/${ref_id}/${gadoga_setting}/inversion_lg01.full.vcf.gz ${gzvcf} LG01:9114741-26192386
        echo " done."
        echo -n "Extracting inversion on LG02..."
        bcftools view --min-ac 1:minor -m2 -M2 -O z -o ../res/gatk/${ref_id}/${gadoga_setting}/inversion_lg02.full.vcf.gz ${gzvcf} LG02:18489307-24048607
        echo " done."
        echo -n "Extracting inversion on LG07..."
        bcftools view --min-ac 1:minor -m2 -M2 -O z -o ../res/gatk/${ref_id}/${gadoga_setting}/inversion_lg07.full.vcf.gz ${gzvcf} LG07:13610591-23019113
        echo " done."
        echo -n "Extracting inversion on LG12..."
        bcftools view --min-ac 1:minor -m2 -M2 -O z -o ../res/gatk/${ref_id}/${gadoga_setting}/inversion_lg12.full.vcf.gz ${gzvcf} LG12:638100-14327837
        echo " done."
        echo -n "Extracting colinear regions..."
        bcftools view --min-ac 1:minor -m2 -M2 -O z -o ../res/gatk/${ref_id}/${gadoga_setting}/colinear.full.vcf.gz ${gzvcf} LG01:1-9114740,LG01:26192387-28298372,LG02:1-18489306,LG03,LG04,LG05,LG06,LG07:1-13610590,LG07:23019114-31229947,LG08,LG09,LG10,LG11,LG12:1-638099,LG12:14327838-27174696,LG13,LG14,LG15,LG16,LG17,LG18,LG19,LG20,LG21,LG22,LG23
        echo " done."
    done
done
