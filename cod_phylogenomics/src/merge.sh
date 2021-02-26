# m_matschiner Fri Nov 22 00:54:07 CET 2019

# Make the log directory.
mkdir -p ../log/misc

# Repeat merging for each of the alternative references.
for ref in ../data/assemblies/gadMor2.fasta ../data/assemblies/??????_in_gadMor2_coords.fasta
do
    # Set the reference id.
    ref_id=`basename ${ref%.fasta}`

    # Rename bam files for specimens that have only one bam file.
    if [ ! -f ../res/mapping/${ref_id}/Gadcha.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadcha_ERR1473886.bam ../res/mapping/${ref_id}/Gadcha.bam
        mv -f ../res/mapping/${ref_id}/Gadcha_ERR1473886.bam.bai ../res/mapping/${ref_id}/Gadcha.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmac.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmac_SRR2906345.bam ../res/mapping/${ref_id}/Gadmac.bam
        mv -f ../res/mapping/${ref_id}/Gadmac_SRR2906345.bam.bai ../res/mapping/${ref_id}/Gadmac.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_avc1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_avc1_AVE1409Z09.bam ../res/mapping/${ref_id}/Gadmor_avc1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_avc1_AVE1409Z09.bam.bai ../res/mapping/${ref_id}/Gadmor_avc1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_avc2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_avc2_AVE1409Z10.bam ../res/mapping/${ref_id}/Gadmor_avc2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_avc2_AVE1409Z10.bam.bai ../res/mapping/${ref_id}/Gadmor_avc2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_avo1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_avo1_AVE1403Z11.bam ../res/mapping/${ref_id}/Gadmor_avo1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_avo1_AVE1403Z11.bam.bai ../res/mapping/${ref_id}/Gadmor_avo1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_avo2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_avo2_AVE1403Z10.bam ../res/mapping/${ref_id}/Gadmor_avo2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_avo2_AVE1403Z10.bam.bai ../res/mapping/${ref_id}/Gadmor_avo2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_bat1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_bat1_BAT1309Z03.bam ../res/mapping/${ref_id}/Gadmor_bat1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_bat1_BAT1309Z03.bam.bai ../res/mapping/${ref_id}/Gadmor_bat1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_bat2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_bat2_BAT1309Z10.bam ../res/mapping/${ref_id}/Gadmor_bat2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_bat2_BAT1309Z10.bam.bai ../res/mapping/${ref_id}/Gadmor_bat2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_bor1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_bor1_BOR1205Z07.bam ../res/mapping/${ref_id}/Gadmor_bor1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_bor1_BOR1205Z07.bam.bai ../res/mapping/${ref_id}/Gadmor_bor1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_bor2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_bor2_BOR1205Z03.bam ../res/mapping/${ref_id}/Gadmor_bor2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_bor2_BOR1205Z03.bam.bai ../res/mapping/${ref_id}/Gadmor_bor2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_icc1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_icc1_ICC0304Z11.bam ../res/mapping/${ref_id}/Gadmor_icc1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_icc1_ICC0304Z11.bam.bai ../res/mapping/${ref_id}/Gadmor_icc1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_ico1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_ico1_ICO0304Z23.bam ../res/mapping/${ref_id}/Gadmor_ico1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_ico1_ICO0304Z23.bam.bai ../res/mapping/${ref_id}/Gadmor_ico1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_ico2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_ico2_ICO0304Z20.bam ../res/mapping/${ref_id}/Gadmor_ico2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_ico2_ICO0304Z20.bam.bai ../res/mapping/${ref_id}/Gadmor_ico2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_kie1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_kie1_KIE1103Z20.bam ../res/mapping/${ref_id}/Gadmor_kie1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_kie1_KIE1103Z20.bam.bai ../res/mapping/${ref_id}/Gadmor_kie1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_kie2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_kie2_KIE1102Z06.bam ../res/mapping/${ref_id}/Gadmor_kie2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_kie2_KIE1102Z06.bam.bai ../res/mapping/${ref_id}/Gadmor_kie2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_lfc1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_lfc1_LOF1106Z04.bam ../res/mapping/${ref_id}/Gadmor_lfc1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_lfc1_LOF1106Z04.bam.bai ../res/mapping/${ref_id}/Gadmor_lfc1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_lfc2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_lfc2_LOF1106Z24.bam ../res/mapping/${ref_id}/Gadmor_lfc2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_lfc2_LOF1106Z24.bam.bai ../res/mapping/${ref_id}/Gadmor_lfc2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_lfo1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_lfo1_LOF1103Z11.bam ../res/mapping/${ref_id}/Gadmor_lfo1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_lfo1_LOF1103Z11.bam.bai ../res/mapping/${ref_id}/Gadmor_lfo1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_lfo2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_lfo2_ERR1551885.bam ../res/mapping/${ref_id}/Gadmor_lfo2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_lfo2_ERR1551885.bam.bai ../res/mapping/${ref_id}/Gadmor_lfo2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_low1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_low1_LOW1503Z06.bam ../res/mapping/${ref_id}/Gadmor_low1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_low1_LOW1503Z06.bam.bai ../res/mapping/${ref_id}/Gadmor_low1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_low2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_low2_LOW1504Z07.bam ../res/mapping/${ref_id}/Gadmor_low2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_low2_LOW1504Z07.bam.bai ../res/mapping/${ref_id}/Gadmor_low2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_twc1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_twc1_TWI1307Z07.bam ../res/mapping/${ref_id}/Gadmor_twc1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_twc1_TWI1307Z07.bam.bai ../res/mapping/${ref_id}/Gadmor_twc1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_twc2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_twc2_TWI1307Z15.bam ../res/mapping/${ref_id}/Gadmor_twc2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_twc2_TWI1307Z15.bam.bai ../res/mapping/${ref_id}/Gadmor_twc2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_two1.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_two1_TWI1307Z12.bam ../res/mapping/${ref_id}/Gadmor_two1.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_two1_TWI1307Z12.bam.bai ../res/mapping/${ref_id}/Gadmor_two1.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadmor_two2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadmor_two2_TWI1307Z17.bam ../res/mapping/${ref_id}/Gadmor_two2.bam
        mv -f ../res/mapping/${ref_id}/Gadmor_two2_TWI1307Z17.bam.bai ../res/mapping/${ref_id}/Gadmor_two2.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadoga.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadoga_SRR2906193.bam ../res/mapping/${ref_id}/Gadoga.bam
        mv -f ../res/mapping/${ref_id}/Gadoga_SRR2906193.bam.bai ../res/mapping/${ref_id}/Gadoga.bam.bai
    fi
    if [ ! -f ../res/mapping/${ref_id}/Gadog2.bam ]
    then
        mv -f ../res/mapping/${ref_id}/Gadog2_ERR1278928.bam ../res/mapping/${ref_id}/Gadog2.bam
        mv -f ../res/mapping/${ref_id}/Gadog2_ERR1278928.bam.bai ../res/mapping/${ref_id}/Gadog2.bam.bai
    fi

    # Merge sets of bam files per specimen.
    if [ ! -f ../res/mapping/${ref_id}/Arcgla.bam ]
    then
        out=../log/misc/merge.${ref_id}_Arcgla.out
        log=../log/misc/merge.${ref_id}_Arcgla.log
        rm -f ${out}
        sbatch -o ${out} merge.slurm ../res/mapping/${ref_id}/Arcgla.bam ${ref} ${log} ../res/mapping/${ref_id}/Arcgla_ERR147388?.bam
    fi
    if [ ! -f ../res/mapping/${ref_id}/Borsai.bam ]
    then
        out=../log/misc/merge.${ref_id}_Borsai.out
        log=../log/misc/merge.${ref_id}_Borsai.log
        rm -f ${out}
        sbatch -o ${out} merge.slurm ../res/mapping/${ref_id}/Borsai.bam ${ref} ${log} ../res/mapping/${ref_id}/Borsai_ERR147388?.bam
    fi
done
