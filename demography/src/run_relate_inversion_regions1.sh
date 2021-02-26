#!/bin/sh

# JMIB [jmi.barth@gmail.com], 2020-06-10
# Wrapper for demography estimation using Relate.

## Copy input files for inversion regions to new names.
cp ../res/relate/input/gadmor.inversion_lg01.adjusted.vcf.gz ../res/relate/input/LG01_inversion.vcf.gz
cp ../res/relate/input/gadmor.inversion_lg02.adjusted.vcf.gz ../res/relate/input/LG02_inversion.vcf.gz
cp ../res/relate/input/gadmor.inversion_lg07.adjusted.vcf.gz ../res/relate/input/LG07_inversion.vcf.gz
cp ../res/relate/input/gadmor.inversion_lg12.adjusted.vcf.gz ../res/relate/input/LG12_inversion.vcf.gz

## Run Relate for each inversion region.
sbatch run_relate_prep_inv.slurm LG01_inversion.vcf.gz
sbatch run_relate_prep_inv.slurm LG02_inversion.vcf.gz
sbatch run_relate_prep_inv.slurm LG07_inversion.vcf.gz
sbatch run_relate_prep_inv.slurm LG12_inversion.vcf.gz