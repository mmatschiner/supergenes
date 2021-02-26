#!/bin/sh

# JMIB [jmi.barth@gmail.com], 2020-06-10
# Wrapper for demography estimation using Relate.

## Run Relate again to infer the demography.
sbatch run_relate_prep_inv_demography.slurm LG01
sbatch run_relate_prep_inv_demography.slurm LG02
sbatch run_relate_prep_inv_demography.slurm LG07
sbatch run_relate_prep_inv_demography.slurm LG12