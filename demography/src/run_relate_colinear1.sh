#!/bin/sh

# JMIB [jmi.barth@gmail.com], 2020-06-10
# Wrapper for demography estimation using Relate.

## Load modules.
module load BCFtools/1.9-foss-2018b

##Index the input files.
for gzvcf in ../res/relate/input/gadmor.colinear.vcf.gz ../res/relate/input/gadmor.inversion_lg??.adjusted.vcf.gz
do
	tabix -p vcf ${gzvcf}
done

## Isolate data per lg.
for lg in LG{03..06} LG{08..11} LG{13..23}
do
	bcftools view ../res/relate/input/gadmor.colinear.vcf.gz ${lg} -h | bgzip > ../res/relate/input/${lg}_colinear.vcf.gz
done

## Run Relate for each colinear chromosome
for chr in ../res/relate/input/*_colinear.vcf.gz
do
	sbatch run_relate_prep.slurm ${chr//..\/res\/relate\/input\//}
done
