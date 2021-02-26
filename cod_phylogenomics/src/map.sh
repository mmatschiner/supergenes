# m_matschiner Wed Nov 20 13:58:50 CET 2019

# Make the output directory.
mkdir -p ../res/mapping/gadMor2
mkdir -p ../res/mapping/melAeg_in_gadMor2_coords

# Make the log directory.
mkdir -p ../log/misc

# Repeat mapping for each of the alternative references.
for ref in ../data/assemblies/gadMor2.fasta ../data/assemblies/melAeg_in_gadMor2_coords.fasta
do
    # Set the reference id.
    ref_id=`basename ${ref%.fasta}`

    # Map each pair of fastq files against the gadmor reference.
    for i in ../data/reads/*_R1.fastq.gz
    do
        species_id=`basename ${i} | cut -d "_" -f 1`
        library_id=`basename ${i} | cut -d "_" -f 2`
        if [ "${species_id}_${library_id}" == "Arcgla_ERR1473882" ]
        then
	        specimen_id=Arcgla
        elif [ "${species_id}_${library_id}" == "Arcgla_ERR1473883" ]
        then
	        specimen_id=Arcgla
        elif [ "${species_id}_${library_id}" == "Borsai_ERR1473884" ]
        then
	        specimen_id=Borsai
        elif [ "${species_id}_${library_id}" == "Borsai_ERR1473885" ]
        then
	        specimen_id=Borsai
        elif [ "${species_id}_${library_id}" == "Gadcha_ERR1473886" ]
        then
	        specimen_id=Gadcha
        elif [ "${species_id}_${library_id}" == "Gadmac_SRR2906345" ]
        then
	        specimen_id=Gadmac
        elif [ "${species_id}_${library_id}" == "Gadmor_BAT1309Z03" ]
        then
	        specimen_id=Gadmor_bat1
        elif [ "${species_id}_${library_id}" == "Gadmor_BAT1309Z10" ]
        then
	        specimen_id=Gadmor_bat2
        elif [ "${species_id}_${library_id}" == "Gadmor_TWI1307Z07" ]
        then
	        specimen_id=Gadmor_twc1
        elif [ "${species_id}_${library_id}" == "Gadmor_TWI1307Z15" ]
        then
	        specimen_id=Gadmor_twc2
        elif [ "${species_id}_${library_id}" == "Gadmor_TWI1307Z12" ]
        then
            specimen_id=Gadmor_two1
        elif [ "${species_id}_${library_id}" == "Gadmor_TWI1307Z17" ]
        then
            specimen_id=Gadmor_two2
        elif [ "${species_id}_${library_id}" == "Gadmor_ICO0304Z23" ]
        then
	        specimen_id=Gadmor_ico1
        elif [ "${species_id}_${library_id}" == "Gadmor_ICO0304Z20" ]
        then
	        specimen_id=Gadmor_ico2
        elif [ "${species_id}_${library_id}" == "Gadmor_ICC0304Z11" ]
        then
	        specimen_id=Gadmor_icc1
        elif [ "${species_id}_${library_id}" == "Gadmor_LOF1103Z11" ]
        then
	        specimen_id=Gadmor_lfo1
        elif [ "${species_id}_${library_id}" == "Gadmor_ERR1551885" ]
        then
	        specimen_id=Gadmor_lfo2
        elif [ "${species_id}_${library_id}" == "Gadmor_LOF1106Z04" ]
        then
	        specimen_id=Gadmor_lfc1
        elif [ "${species_id}_${library_id}" == "Gadmor_LOF1106Z24" ]
        then
            specimen_id=Gadmor_lfc2
        elif [ "${species_id}_${library_id}" == "Gadmor_LOF1106Z11" ]
        then
	        specimen_id=Gadmor_lfc3
        elif [ "${species_id}_${library_id}" == "Gadmor_AVE1409Z09" ]
        then
	        specimen_id=Gadmor_avc1
        elif [ "${species_id}_${library_id}" == "Gadmor_AVE1409Z10" ]
        then
	        specimen_id=Gadmor_avc2
        elif [ "${species_id}_${library_id}" == "Gadmor_AVE1403Z11" ]
        then
	        specimen_id=Gadmor_avo1
        elif [ "${species_id}_${library_id}" == "Gadmor_AVE1403Z10" ]
        then
	        specimen_id=Gadmor_avo2
        elif [ "${species_id}_${library_id}" == "Gadmor_LOW1503Z06" ]
        then
	        specimen_id=Gadmor_low1
        elif [ "${species_id}_${library_id}" == "Gadmor_LOW1504Z07" ]
        then
            specimen_id=Gadmor_low2
        elif [ "${species_id}_${library_id}" == "Gadmor_KIE1103Z20" ]
        then
            specimen_id=Gadmor_kie1
        elif [ "${species_id}_${library_id}" == "Gadmor_KIE1102Z06" ]
        then
            specimen_id=Gadmor_kie2
        elif [ "${species_id}_${library_id}" == "Gadmor_BOR1205Z07" ]
        then
            specimen_id=Gadmor_bor1
        elif [ "${species_id}_${library_id}" == "Gadmor_BOR1205Z03" ]
        then
            specimen_id=Gadmor_bor2
        elif [ "${species_id}_${library_id}" == "Gadoga_SRR2906193" ]
        then
	        specimen_id=Gadoga
        elif [ "${species_id}_${library_id}" == "Gadoga_ERR1278928" ]
        then
            specimen_id=Gadog2
        else
	        echo "ERROR: Unexpected combination of species id and library id: ${species_id}_${library_id}!"
	        exit 1
        fi
        fastq1=${i}
        fastq2=${i%_R1.fastq.gz}_R2.fastq.gz
        bam=../res/mapping/${ref_id}/${specimen_id}_${library_id}.bam
        out=../log/misc/map.${ref_id}_${specimen_id}_${library_id}.out
        log=../log/misc/map.${ref_id}_${specimen_id}_${library_id}.log
        rm -f ${out}
        if [ ${specimen_id} == "Gadmor_lfo2" ]
        then
            sbatch -o ${out} --partition=bigmem --mem-per-cpu=10G map.slurm ${specimen_id} ${ref} ${bam} ${log} ${fastq1} ${fastq2}
        else
            echo sbatch -o ${out} map.slurm ${specimen_id} ${ref} ${bam} ${log} ${fastq1} ${fastq2}
        fi
    done
done
